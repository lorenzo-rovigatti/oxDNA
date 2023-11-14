import sys
from os import path
from typing import Dict, Tuple
from collections import namedtuple
import argparse
import numpy as np
import matplotlib.pyplot as plt
import oxpy
from oxDNA_analysis_tools.UTILS.RyeReader import describe
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.data_structures import TrajInfo

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "input_file",
                                              "n1",
                                              "n2"])

def get_r(conf, nucid:int, pair_dict:Dict) -> np.ndarray:
    """
        Returns a normalized vector pointing from base midpoint of the `nucid`-th base pair to the midpoint of the next base pair.
        
        Parameters:
            conf (oxpy.config_info) : The current configuration
            nucid (int) : ID of the nucleotide in the first strand to compute the vector from

        Returns:
            (np.ndarray) : The vector pointing from the midpoint of nucid's base pair to the +1 base pair.
    """
    box = np.array(conf.box.box_sides)
    pair = pair_dict[nucid]
    next_pair = pair_dict[nucid+1]

    firstA = conf.particles()[nucid].base_site()
    firstB = conf.particles()[pair].base_site()
    secondA = conf.particles()[nucid+1].base_site()
    secondB = conf.particles()[next_pair].base_site()

    first_midpos = (firstA + firstB) / 2
    second_midpos = (secondA + secondB) / 2 

    r = second_midpos - first_midpos 		
    r -= box * np.rint(r / box)
    return r

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    l0 = 0
    correlations = np.zeros(ctx.n2 - ctx.n1)
    correlations_counter = np.zeros_like(correlations)
    # Read in a conf and identify h-bonds
    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(ctx.input_file)
        inp["list_type"] = "cells"
        inp["trajectory_file"] = ctx.traj_info.path
        inp["analysis_bytes_to_skip"] = str(ctx.traj_info.idxs[chunk_id*chunk_size].offset)
        inp["confs_to_analyse"] = str(chunk_size)
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = hb_list \n } \n }'

        if (not inp["use_average_seq"] or inp.get_bool("use_average_seq")) and "RNA" in inp["interaction_type"]:
            print("WARNING: Sequence dependence not set for RNA model, wobble base pairs will be ignored", file=sys.stderr)

        backend = oxpy.analysis.AnalysisBackend(inp)
        if backend.config_info().particles()[ctx.n1].strand_id != backend.config_info().particles()[ctx.n2].strand_id:
            raise RuntimeError("nucid_1 and nucid_2 must be on the same strand")

        while backend.read_next_configuration():
            pairs = backend.config_info().get_observable_by_id("my_obs").get_output_string(backend.config_info().current_step).strip().split('\n')[1:]

            # Extract paired nucleotides into a dict
            pair_dict = {}
            for p in pairs:
                p1 = int(p.split()[0])
                p2 = int(p.split()[1])
                if p1 < ctx.n1 or p1 > ctx.n2:
                    continue
                pair_dict[p1] = p2


            for j in range(ctx.n1, ctx.n2):
                # If there's no base pair, there's no midpoint
                if not j in pair_dict or not j+1 in pair_dict:
                    print("WARNING: Nucleotide {} or {} is unpaired.  Skipping...".format(j, j+1))
                    continue

                # Get the midpoint of the base pairs
                r0 = get_r(backend.config_info(), j, pair_dict)
                l0 += np.linalg.norm(r0)
                r0 = r0 / np.linalg.norm(r0)

                # Get the correlation between r0 and each subsequence base step
                for k in range(j, ctx.n2):
                    if not k in pair_dict or not k+1 in pair_dict:
                        continue

                    rk = get_r(backend.config_info(), k, pair_dict)
                    rk = rk / np.linalg.norm(rk)
                    correlations[k-j] += np.dot(r0, rk)
                    correlations_counter[k-j] += 1

    return(l0, correlations, correlations_counter)

def persistence_length(traj_info:TrajInfo, inp_file:str, n1:int, n2:int, ncpus:int=1) -> Tuple[float, np.ndarray]:
    """
        Computes the persistence length of a bonded sequence of nucleotides.

        Note that while n1 and n2 must be on the same strand, there can be multiple distinct duplexes within the strand.

        Only paired nucleotides will be considered in the persistence length calculation

        Parameters:
            traj_info (TrajInfo): TrajInfo object for the trajectory you want to analyze
            inp_file (str): The path to the input file used to run the simulation
            n1 (int): ID of the particle to start the analysis
            n2 (int): ID of the particle to end the analysis. 
        
        Returns:
            (Tuple[float, np.ndarray]) : Tuple containing the average contour length and correlations between each pair step
    """

    ctx = ComputeContext(traj_info, inp_file, n1, n2)
    l0 = 0
    correlations = np.zeros(n2 - n1)
    correlations_counter = np.zeros_like(correlations)

    def callback(i, r):
        nonlocal l0, correlations, correlations_counter
        l0 += r[0]
        correlations += r[1]
        correlations_counter += r[2]

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    l0 /= traj_info.nconfs
    correlations /= correlations_counter
    return(l0, correlations)

def fit_PL(correlations:np.ndarray, plt_name:str) -> float:
    """
        Fits persistence length from tangent vector correlations

        Parameters:
            correlations (np.ndarray) : Array of offsets vs correlation
            plt_name (str) : Name to save the resulting plot to

        Returns:
            (float) : Persistence length in nucleotides
    """
    # Fit the PL to the correlations
    x = np.arange(0, len(correlations))
    log_corr = np.log(correlations)
    A, B = np.polyfit(x, log_corr, 1)
    pl = -1/A

    # Make a plot
    fig, ax = plt.subplots()
    ax.scatter(x, log_corr, alpha=0.5)
    trend = np.poly1d([A, B])
    ax.plot(x, trend(x), label=f'Fit y = (-x/{pl:.1f})+{B:.3f}')
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_xlabel('Offset')
    ax.set_ylabel('ln(correlation)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(plt_name, dpi=300)
    print("INFO: Saving figure to", plt_name, file=sys.stderr)

    return pl

def cli_parser(prog="persistence_length.py"):
    parser = argparse.ArgumentParser(prog=prog, description="Calculates persistence length and contour length of a paired sequence of DNA.")
    parser.add_argument('traj_file', type=str, nargs=1, help="The trajectory file to analyze")
    parser.add_argument('input', type=str, nargs=1, help="The input file associated with the trajectory")
    parser.add_argument('nucid_1', type=int, nargs=1, help="Nucleotide ID to start the calculation at")
    parser.add_argument('nucid_2', type=int, nargs=1, help="Nucleotide ID to stop the calculation at")
    parser.add_argument('-p', '--parallel', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-d', '--data', metavar='data_file', nargs=1, help='If set, the correlations will be written to a txt file in the format `offset correlation`')
    parser.add_argument('-n', '--plot_name', nargs=1, help='Name to save the plot showing the fit of persistence length to correlations.  Defaults to persistence_length.png')
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python", "matplotlib", "numpy"])

    traj = args.traj_file[0]
    inp_file = args.input[0]
    n1 = args.nucid_1[0]
    n2 = args.nucid_2[0]

    # -p sets the number of cpus to use
    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    ti, di = describe(None, traj)

    l0, correlations = persistence_length(di, inp_file, n1, n2, ncpus)

    if args.data:
        print("INFO: Writing correlations to", args.data[0], file=sys.stderr)
        with open(args.data[0], 'w+') as f:
            for i, c in enumerate(correlations):
                f.write(f'{i} {c}\n')

    if args.plot_name:
        plot_name = args.plot_name[0]
    else:
        plot_name = 'persistence_length.png'

    pl = fit_PL(correlations, plot_name)
    print("Persistence length: {:.1f} nucleotides".format(pl))
    print("Overall bonded contour length between n1 and n2 is:", l0*0.8518, "nm")

    print("--- %s seconds ---" % (time.time() - start_time))
                
if __name__ == '__main__':
    main()
