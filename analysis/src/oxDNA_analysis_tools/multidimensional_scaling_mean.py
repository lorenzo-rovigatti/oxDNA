#!/usr/bin/env python3
#Multidimensional_scaling_mean
#Written by: Erik Poppleton
#Date: 3/28/2022
#Computes the RMSD of the contact map for a structure.  The average structure is determined
#by Scikit.learn's MDS algorithm, then subtracts the contact map of each individual structure from the man
#This is used to compute a per-nucleotide deviation in the contact map, which can be visualized with oxView

import numpy as np
import argparse
from os import path
from sys import exit, stderr
from json import dumps
from collections import namedtuple
from typing import Tuple
from sklearn.manifold import MDS
from oxDNA_analysis_tools.config import check
from oxDNA_analysis_tools.contact_map import contact_map
from oxDNA_analysis_tools.distance import vectorized_min_image
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, TopInfo, TrajInfo

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info"])

DevsContext = namedtuple("DevsContext",["traj_info",
                                        "top_info",
                                        "masked_mean_coords"])

#at 2.5 you start to see the hard edges caused by end-loops and see some loop interactions
CUTOFF = 2.5

def make_heatmap(contact_map:np.ndarray):
    """
    Convert a matrix of contact distances to a visual contact map.

    Parameters:
        contact_map (numpy.array): An array of all pairwise distances between nucleotides.
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    a = ax.imshow(contact_map, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
       ylabel="nucleotide id",
       xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance", rotation = 270)
    plt.show()

def devs_mds(ctx:DevsContext, chunk_size:int, chunk_id:int, ): 
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    
    np_poses = np.asarray([c.positions for c in confs])

    devs = np.zeros((ctx.top_info.nbases, ctx.top_info.nbases))

    for c in np_poses:
        c_map = vectorized_min_image(c, c, confs[0].box[0])
        masked_distances = np.ma.masked_array(c_map, ~(c_map < CUTOFF))

        # Fill the masked values with the cutoff.  Not sure if this is the best practice here.
        masked_distances = np.ma.filled(masked_distances, CUTOFF)
        masked_mean = np.ma.filled(ctx.masked_mean_coords, CUTOFF)

        diff = masked_distances - masked_mean
        diff = np.square(diff)
        devs += diff

    return devs

def multidimensional_scaling_mean(traj_info:TrajInfo, top_info:TopInfo, ncpus:int=1) -> Tuple[Configuration, np.ndarray]:
    """
        Compute the mean configuration of a trajectory using MDS.

        Parameters:
            traj_info (TrajInfo): Information about the trajectory.
            top_info (TopInfo): Information about the topology.
            ncpus (int): (optional) Number of CPUs to use.
    """
    example_conf = get_confs(top_info, traj_info, 1, 1)[0]
    
    distances = contact_map(traj_info, top_info, ncpus)

    mean_distances = distances / traj_info.nconfs
    masked_mean = np.ma.masked_array(mean_distances, ~(mean_distances < CUTOFF))

    print("INFO: fitting local distance data", file=stderr)
    mds = MDS(n_components=3, metric=True, max_iter=3000, eps=1e-12, dissimilarity="precomputed", n_jobs=1, n_init=1)
    out_coords = mds.fit_transform(masked_mean, init=example_conf.positions) #without the init you can get a left-handed structure.
    a1s = np.zeros((top_info.nbases, 3))
    a3s = np.zeros((top_info.nbases, 3))
    mean_conf = Configuration(0,example_conf.box, np.array([0,0,0]), out_coords, a1s , a3s)

    return mean_conf, masked_mean

def distance_deviations(traj_info:TrajInfo, top_info:TopInfo, masked_mean:np.ndarray, ncpus:int=1) -> np.ndarray:
    """
        Compute the deviations of the contact map from the mean.

        Parameters:
            traj_info (TrajInfo): Information about the trajectory.
            top_info (TopInfo): Information about the topology.
            masked_mean (np.ndarray): The mean contact map with values greater than a cutoff masked.
            ncpus (int): (optional) Number of CPUs to use.
    """
    # Compute the deviations from the mean
    ctx = DevsContext(traj_info, top_info, masked_mean)
    
    devs = np.zeros((top_info.nbases, top_info.nbases))
    def callback(i, r):
        nonlocal devs
        devs += r

    oat_multiprocesser(traj_info.nconfs, ncpus, devs_mds, callback, ctx)

    devs = np.ma.masked_array(devs, ~(devs != 0.0))
    devs = devs / traj_info.nconfs
    devs = np.mean(devs, axis=0)
    devs = np.sqrt(devs)

    return devs

def cli_parser(prog="multidimensional_scaling_mean.py"):
    #get commandline arguments
    parser = argparse.ArgumentParser(prog = prog, description="Calculate molecular contacts, and assembles an average set of contacts based on MDS")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-o', '--output', metavar='output', type=str, nargs=1, help='the name of the .dat file where the mean will be written')
    parser.add_argument('-d', '--dev_file', metavar='dev_file', type=str, nargs=1, help='the name of the .json file where the devs will be written')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()
    traj = args.trajectory[0]
    top_info, traj_info = describe(None, traj)

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    check(['python', 'numpy'])

    mean_conf, masked_mean = multidimensional_scaling_mean(traj_info, top_info, ncpus)

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else:
        outfile = "mean_mds.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    write_conf(outfile,mean_conf, include_vel=traj_info.incl_v)
    print("INFO: Wrote mean to {}".format(outfile), file=stderr)

    devs = distance_deviations(traj_info, top_info, masked_mean, ncpus)

    #-d names the deviations file
    if args.dev_file:
        devfile = args.dev_file[0].split(".")[0] + ".json"
    else:
        devfile = "devs_mds.json"
        print("INFO: No deviations file name provided, defaulting to \"{}\"".format(devfile), file=stderr)

    with open(devfile, "w") as file:
        file.write(
            dumps({"contact deviation" : list(devs)})
        )
    print("INFO: wrote file {}".format(devfile), file=stderr)

if __name__ == '__main__':
    main()