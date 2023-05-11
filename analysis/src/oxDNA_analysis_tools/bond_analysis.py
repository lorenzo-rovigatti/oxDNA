import argparse
from sys import exit, stderr
from os import path
from collections import namedtuple
from typing import Tuple, Dict
import numpy as np
import matplotlib.pyplot as plt
import oxpy
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser, get_chunk_size
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_input_parameter

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "designed_pairs",
                                              "input_file"])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
        count_correct_bonds = []
        count_incorrect_bonds = []
        tot_bonds = []
        out_array = np.zeros(ctx.top_info.nbases, dtype=int)
        while backend.read_next_configuration():
            conf_corr_bonds = 0
            conf_incorr_bonds = 0
            conf_tot_bonds = 0
            pairs = backend.config_info().get_observable_by_id("my_obs").get_output_string(backend.config_info().current_step).strip().split('\n')
            for p in pairs[1:]:
                p = p.split()
                a = int(p[0])
                b = int(p[1])
                if a in ctx.designed_pairs.keys():
                    if ctx.designed_pairs[a] == b:
                        conf_corr_bonds += 1
                        conf_tot_bonds += 1
                        out_array[a] += 1
                        out_array[b] += 1
                    else:
                        conf_incorr_bonds += 1
                        conf_tot_bonds += 1

            count_correct_bonds.append(conf_corr_bonds)
            count_incorrect_bonds.append(conf_incorr_bonds)
            tot_bonds.append(conf_tot_bonds)

        return(np.array(tot_bonds), np.array(count_correct_bonds), np.array(count_incorrect_bonds), out_array)


def bond_analysis(traj_info:TrajInfo, top_info:TopInfo, pairs:Dict[int, int], inputfile:str, ncpus:int=1) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    '''
        Compare the bond occupancy of a trajectory with a designed structure

        Parameters: 
            traj_info (TrajInfo): Object containing the trajectory information
            top_info (TopInfo): Object containing the topology information
            pairs (dict): Designed pairs ({p1 : q1, p2 : q2})
            inputfile (str): The path to the input file used to run the simulation
            ncpus (int): (optional) number of cores to use

        Returns:
            total_bonds (np.ndarray): Number of formed bonds among the specified nucleotides at each step in the simulation
            incorrect_bonds (np.ndarray): Number of missbonds among specified nucleotides at each step in the simulation
            correct_bonds (np.ndarray): Number of correct bonds among specified nucleotides at each step in the simulation
            nt_array (np.ndarray): per-nucleotide correct bond occupancy
    '''
    ctx = ComputeContext(traj_info, top_info, pairs, inputfile)

    total_bonds = np.empty(traj_info.nconfs)
    correct_bonds = np.empty(traj_info.nconfs)
    incorrect_bonds = np.empty(traj_info.nconfs)
    nt_array = np.zeros(ctx.top_info.nbases, dtype=int)

    chunk_size = get_chunk_size()
    def callback(i, r):
        nonlocal total_bonds, correct_bonds, incorrect_bonds, nt_array
        total_bonds[i*chunk_size:i*chunk_size+len(r[0])] = r[0]
        correct_bonds[i*chunk_size:i*chunk_size+len(r[1])] = r[1]
        incorrect_bonds[i*chunk_size:i*chunk_size+len(r[2])] = r[2]
        nt_array += r[3]

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    nt_array  = nt_array / traj_info.nconfs

    return(total_bonds, correct_bonds, incorrect_bonds, nt_array)

def oxView_overlay(nt_array:np.ndarray, outfile:str):
    print("INFO: Writing bond occupancy data to {}".format(outfile))
    with open(outfile, "w+") as file:
        file.write("{\n\"occupancy\" : [")
        file.write(str(nt_array[0]))
        for n in nt_array[1:]:
            file.write(", {}".format(n))
        file.write("] \n}")
    return

def plot_trajectories(correct_bonds:np.ndarray, incorrect_bonds:np.ndarray, designed_bonds:int, plotname:str):
    fig, ax = plt.subplots()
    ax.plot(correct_bonds, alpha=0.5, label="Correct bonds")
    ax.plot(incorrect_bonds, alpha=0.5, label="Incorrect bonds")
    ax.axhline(designed_bonds, alpha=0.3, color='red', label="Designed bonds")
    plt.xlabel('Configuration')
    plt.ylabel('Number of Bonds')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plotname)
    return

def cli_parser(prog="bond_analysis.py"):
    #read data from files
    parser = argparse.ArgumentParser(prog = prog, description="Compare the bonds found at each trajectory with the intended design")
    parser.add_argument('inputfile', type=str, help="The inputfile used to run the simulation")
    parser.add_argument('trajectory', type=str, help="The trajecotry file to compare against the designed pairs")
    parser.add_argument('designed_pairs', type=str, help="The file containing the desired nucleotides pairings in the format `a b`")
    parser.add_argument('-o', metavar='output_file', type=str, dest='outfile', help="Name of the file to save the output oxView overlay to")
    parser.add_argument('-t', metavar='trajectory_plot', type=str, dest='traj_plot', help='Name of the file to save the trajecotry plot to')
    parser.add_argument('-p', metavar='num_cpus', type=int, dest='parallel', help="(optional) How many cores to use")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    #run system checks
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    # Parse CLI input
    inputfile = args.inputfile
    traj_file = args.trajectory
    designfile = args.designed_pairs
    if args.outfile:
        outfile = args.outfile
        outfile = outfile.strip(".json")+".json"
    else:
        outfile = 'bonds.json'
        print("INFO: No oxView name provided, defaulting to \"{}\"".format(outfile), file=stderr)
    if args.traj_plot:
        plotfile = args.traj_plot
        plotfile = plotfile.strip(".png")+".png"
    else:
        plotfile = 'bonds.png'
        print("INFO: No bond plot name provided, defaulting to \"{}\"".format(plotfile), file=stderr)

    # Get trajectory metadata
    top_file = get_input_parameter(inputfile, "topology")
    if not path.exists(top_file):
        raise RuntimeError("Topology '{}' not found, the topology specified in the input file must be present.".format(top_file))
    
    top_info, traj_info = describe(top_file, traj_file)

    # Convert the designed pairs into a dict
    with open(designfile, 'r') as file:
        pairs_txt = file.readlines()

    pairs = {int(p[0]) : int(p[1]) for p in [p.split() for p in pairs_txt]}

    if args.parallel:
        ncpus = args.parallel
    else:
        ncpus = 1

    # Compute bond information
    total_bonds, correct_bonds, incorrect_bonds, nt_array = bond_analysis(traj_info, top_info, pairs, inputfile, ncpus)

    # Summarize output and generate output filess
    print("\nSummary:\navg bonds: {}\navg correct bonds: {}/{}\navg missbonds: {}".format(np.mean(total_bonds),np.mean(correct_bonds), len(pairs), np.mean(incorrect_bonds)))

    oxView_overlay(nt_array, outfile)

    plot_trajectories(correct_bonds, incorrect_bonds, len(pairs), plotfile)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()