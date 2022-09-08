import numpy as np
import argparse
from os import path
from sys import stderr
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe
import oxpy

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "input_file",
                                              "visualize",
                                              "conversion_factor"])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id*chunk_size].offset)
        inp["confs_to_analyse"] = str(chunk_size)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = pair_energy \n } \n }'

        if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
            print("WARNING: Sequence dependence not set for RNA model, wobble base pairs will be ignored", file=stderr)

        backend = oxpy.analysis.AnalysisBackend(inp)

        # The 8 energies are:
        # 0 fene
        # 1 bexc
        # 2 stack
        # 3 nexc
        # 4 hb
        # 5 cr_stack
        # 6 cx_stack
        # 9 Debye-Huckel
        # 7 total
        if ctx.visualize:
            energies = np.zeros((ctx.top_info.nbases, 9))

        while backend.read_next_configuration():
            e_txt = backend.config_info().get_observable_by_id("my_obs").get_output_string(0).strip().split('\n')
            if ctx.visualize:
                for e in e_txt[1:]:
                    if not e[0] == '#':
                        e = e.split()
                        p = int(e[0])
                        q = int(e[1])
                        l = np.array([float(x) for x in e[2:]])*ctx.conversion_factor
                        energies[p] += l
                        energies[q] += l
            else:
                print(e_txt[0])
                for e in e_txt[1:]:
                    if not e[0] == '#':
                        e = e.split()
                        p = int(e[0])
                        q = int(e[1])
                        l = np.array([float(x) for x in e[2:]])*ctx.conversion_factor
                        print("{} {} {} {} {} {} {} {} {} {} {}".format(p, q, l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7], l[8]))
                    else: 
                        print(e)
        if ctx.visualize:
            return energies
        else:
            return
                
def output_bonds(traj_info:TrajInfo, top_info:TopInfo, inputfile:str, visualize:bool=False, conversion_factor:float=1, ncpus:int=1):
    """
        Computes the potentials in a trajectory

        Parameters:
            traj_info (TrajInfo): Information about the trajectory.
            top_info (TopInfo): Information about the topology.
            inputfile (str): Path to the input file.
            visualize (bool): (optional) If True, the energies are saved as a mean-per-particle oxView file.  If False, they are printed to the screen.
            conversion_factor (float): (optional) Conversion factor for the energies. 1 for oxDNA SU, 41.42 for pN nm.
            ncpus (int): (optional) Number of CPUs to use.

        Returns:
            energies (np.array): If visualize is True, the energies are saved as a mean-per-particle oxView file.  If False, they are printed to the screen and None is returned.
    """
    
    ctx = ComputeContext(traj_info, top_info, inputfile, visualize, conversion_factor)

    energies = np.zeros((ctx.top_info.nbases, 9))
    def callback(i, r):
        nonlocal visualize, energies
        if visualize:
            energies += r
        else:
            print(r)

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    if visualize:
        return energies
    else:
        return None

def cli_parser(prog="output_bonds.py"):
    parser = argparse.ArgumentParser(prog = prog, description="List all the interactions between nucleotides")
    parser.add_argument('inputfile', type=str, nargs=1, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-v', type=str, nargs=1, dest='outfile', help='if you want instead average per-particle energy as an oxView JSON')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-u', '--units', type=str, nargs=1, dest='units', help="(optional) The units of the energy (pNnm or oxDNA)")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    traj_file = args.trajectory[0]
    inputfile = args.inputfile[0]

    top_info, traj_info  = describe(None, traj_file)

    try:
        outfile = args.outfile[0]
        visualize = True
    except:
        visualize = False

    #if path.dirname(inputfile) != getcwd():
    #    sim_directory = path.dirname(inputfile)
    #else:
    #    sim_directory = ""

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    if args.units:
        if args.units[0] == "pNnm":
            units = "pN nm"
            conversion_factor = 41.42
        elif args.units[0] == "oxDNA":
            units = "oxDNA su"
            conversion_factor = 1
        else:
            print("Unrecognized units:", args.units[0], file=stderr)
            exit(1)
    else:
        units = "oxDNA su"
        conversion_factor = 1
        print("INFO: no units specified, assuming oxDNA su", file=stderr)

    energies = output_bonds(traj_info, top_info, inputfile, visualize, conversion_factor, ncpus)

    if visualize:
        energies /= traj_info.nconfs
        for i, potential in enumerate(["FENE","bexc", "stack", "nexc", "hb", "cr_stack", "cx_stack", "Debye-Huckel", "Total"]):
            if '.json' in outfile:
                fname = '.'.join(outfile.split('.')[:-1])+"_"+potential+'.json'
            else:
                fname = outfile+"_"+potential+'.json'
            with open(fname, 'w+') as f:
                f.write("{{\n\"{} ({})\" : [".format(potential, units))
                f.write(', '.join([str(x) for x in energies[:,i]]))
                f.write("]\n}")
            print("INFO: Wrote oxView overlay to:", fname)


if __name__ == "__main__":
    main()