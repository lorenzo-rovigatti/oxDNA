from typing import List
import numpy as np
import argparse
from os import path
from sys import stderr
from dataclasses import dataclass
from collections import namedtuple
import oxpy
from oxDNA_analysis_tools.UTILS import geom
from oxDNA_analysis_tools.UTILS.data_structures import Monomer, TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe, strand_describe, get_input_parameter

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "input_file",
                                              "monomers"])

@dataclass
class Duplex:
    """
        Defines a nucleic acid duplex structure

        Parameters:
            index (int): Unique identifier for the duplex
            start1 (int): Particle ID of the first particle in the first strand
            end1 (int): Particle ID of the last particle in the first strand
            start2 (int): Particle ID of the first particle in the complementary strand
            end2 (int): Particle ID of the last particle in the complementary strand
            axis (np.array): Normalized vector fit to the center of the duplex
            pos (np.array): start position in 3D space of the axis
    """
    time : int
    index : int
    start1 : int
    end1 : int
    start2 : int
    end2 : int
    axis : np.ndarray
    pos : np.ndarray


def find_duplex(monomers:List[Monomer]) -> List[Duplex]:
    """
        Finds the duplexes in a structure

        Parameters:
            monomers (List[Monomer]): List of monomers in the structure

        Returns:
            List[Duplex]: List of duplexes
    """
    def terminating_conditions(m):
        if m.id in assigned_monomers: # already assigned this nucleotide to another duplex
            return True
        if m.n3 is None: # this nucleotide is the first in a strand
            return True
        if m.pair is None: # this nucleotide is not part of a duplex
            return True
        if monomers[m.pair].n5 != monomers[m.n3].pair: # this nucleotide is at the start of a duplex
            return True
        if monomers[m.n3].pair is None: #special case to catch 3' overhangs
            return True
        return False

    duplex_index = 0
    duplex_length = 0
    duplex_list = []
    assigned_monomers = set({0}) # set of IDs because objects are not hashable.  Start with 0 filled to create an initial duplex if paired.
    d = Duplex(0, 0, 0, 0, 0, 0, np.zeros(3), np.zeros(3)) #initialize d so it is in scope

    for m in monomers:
        # the current duplex is over for some reason
        if terminating_conditions(m):
            if duplex_length > 3: # we only consider duplexes of length 4 or more
                d.end1 = monomers[m.id - 1].id
                d.start2 = monomers[m.id - 1].pair
                duplex_list.append(d)
                duplex_index += 1
            
            # reset for next round
            if m.pair is not None:
                d = Duplex(0, duplex_index, m.id, 0, 0, m.pair, np.zeros(3), np.zeros(3))
            duplex_length = 0

        duplex_length += 1
        assigned_monomers.add(m.id)
        if m.pair is not None:
            assigned_monomers.add(m.pair)

    return duplex_list

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):

    duplexes_at_step = []

    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id*chunk_size].offset)
        inp["confs_to_analyse"] = str(chunk_size)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = hb_list \n } \n }'

        if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
            print("WARNING: Sequence dependence not set for RNA model, wobble base pairs will be ignored", file=stderr)

        backend = oxpy.analysis.AnalysisBackend(inp)
        i = 0

        while backend.read_next_configuration():
            pairs = backend.config_info().get_observable_by_id("my_obs").get_output_string(backend.config_info().current_step).strip().split('\n')
            for p in pairs[1:]:
                p = p.split()
                a = int(p[0])
                b = int(p[1])
                ctx.monomers[a].pair = ctx.monomers[b].id
                ctx.monomers[b].pair = ctx.monomers[a].id

            duplex_list = find_duplex(ctx.monomers)

            if "RNA" in inp["interaction_type"]:
                for d in duplex_list:
                    d.axis, d.pos = geom.get_RNA_axis(backend.config_info().particles(), d)
            else:
                for d in duplex_list:
                    d.axis, d.pos = geom.get_DNA_axis(backend.config_info().particles(), d)

            duplexes_at_step.append(duplex_list)
            i +=1

        return duplexes_at_step

def duplex_finder(traj_info:TrajInfo, top_info:TopInfo, inputfile:str, monomers:List[Monomer], ncpus=1) -> List[List[Duplex]]:
    """
        Finds the duplexes in a trajectory

        Parameters:
            traj_info (TrajInfo): Information about the trajectory
            top_info (TopInfo): Information about the topology
            inputfile (str): Path to the input file
            monomers (List[Monomer]): List of monomers in the structure
            ncpus (int): Number of CPUs to use

        Returns:
            List[List[Duplex]]: List of lists of duplexes, one list for each step in the trajectory
    """
    ctx = ComputeContext(traj_info, top_info, inputfile, monomers)

    duplexes_at_step = []
    def callback(i, r):
        duplexes_at_step.extend(r)

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    return duplexes_at_step

def cli_parser(prog="duplex_finder.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Fit vectors to every duplex in the structure")
    parser.add_argument('input', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help="The trajectory file from the simulation")
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file',  type=str, nargs=1, help='name of the file to write the angle list to')
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    #Process command line arguments:
    inputfile = args.input[0]
    traj_file = args.trajectory[0]

        #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        outfile = "angles.txt"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    top_file = get_input_parameter(inputfile, "topology")
    top_info, traj_info = describe(top_file, traj_file)
    system, monomers = strand_describe(top_file)

    duplexes_at_step = duplex_finder(traj_info, top_info, inputfile, monomers, ncpus)

    #print duplexes to a file
    print("INFO: Writing duplex data to {}.  Use duplex_angle_plotter to graph data".format(outfile), file=stderr)
    with open(outfile, 'w') as f:
        f.write("time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n")
        for i in range (0, len(duplexes_at_step)):
            for j in range(0, len(duplexes_at_step[i])):
                line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t[{},{},{}]\n'.format(i,duplexes_at_step[i][j].index,duplexes_at_step[i][j].start1,duplexes_at_step[i][j].end1,duplexes_at_step[i][j].start2,duplexes_at_step[i][j].end2,duplexes_at_step[i][j].axis[0],duplexes_at_step[i][j].axis[1],duplexes_at_step[i][j].axis[2],duplexes_at_step[i][j].pos[0],duplexes_at_step[i][j].pos[1],duplexes_at_step[i][j].pos[2])
                f.write(line)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()