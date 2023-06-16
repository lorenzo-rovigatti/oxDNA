#!/usr/bin/env python

import numpy as np
from sys import exit, stderr
from collections import namedtuple
from typing import List
import argparse
import os
import matplotlib.pyplot as plt
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import get_chunk_size, oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs

import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info", 
                                              "p1s",
                                              "p2s"])

def min_image(p1:np.ndarray, p2:np.ndarray, box:float) -> float:
    """
    Calculates distance between two particles taking PBC into account

    Parameters:
        p1 (np.ndarray): The first particle's position
        p2 (np.ndarray): The second particle's position
        box (float): The size of the box (assumes a cubic box)

    Returns:
        (float): The distance between the two particles
    """
    p1 = p1 - (np.floor(p1/box) * box)
    p2 = p2 - (np.floor(p2/box) * box)
    diff = p1 - p2
    diff = diff - (np.round(diff/box)*box)
    return float(np.linalg.norm(diff))

def vectorized_min_image(p1s:np.ndarray, p2s:np.ndarray, box:float) -> np.ndarray:
    """
    Calculates all mutual distances between two sets of points taking PBC into account
    
    Paramters:
        p1s (np.ndarray) : the first set of points (Nx3 array)
        p2s (np.ndarray) : the second set of points (Mx3 array)
        box (float) : The size of the box (assumes a cubic box)

    returns:
        (np.array) : the distances between the points (NxM array)
    """

    p1s = p1s - (np.floor(p1s/box) * box)
    p2s = p2s - (np.floor(p2s/box) * box)
    diff = p1s[np.newaxis,:,:] - p2s[:,np.newaxis,:]
    diff = diff - (np.round(diff/box)*box)
    return np.linalg.norm(diff, axis=2)

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    box = confs[0].box
    distances = np.empty((len(ctx.p1s), len(confs)))

    for i, conf in enumerate(confs):
        distances[:,i] = [min_image(conf.positions[p1], conf.positions[p2], box)* 0.85 for p1, p2 in zip(ctx.p1s, ctx.p2s)]
    
    return distances

def distance(traj_infos:List[TrajInfo], top_infos:List[TopInfo], p1ss:List[List[int]], p2ss:List[List[int]], ncpus:int=1) -> List[List[float]]:
    """
        Compute the distance between two lists of particles

        Parameters:
            traj_infos (List[TrajInfo]): A list of TrajInfo objects
            top_infos (List[TopInfo]): A list of TopInfo objects
            p1ss (List[List[int]]): A list of particle indices for each trajectory
            p2ss (List[List[int]]): A list of particle indices for each trajectory

        Returns:
            distances (List[List[float]]): A list of distances for each trajectory
    """
    distances = [[] for _ in traj_infos]
    for i, (traj_info, top_info, p1s, p2s) in enumerate(zip(traj_infos, top_infos, p1ss, p2ss)):
        
        ctx = ComputeContext(traj_info, top_info, p1s, p2s)
        
        chunk_size = get_chunk_size()
        distances[i] = [[None]*traj_info.nconfs for _ in p1s]
        def callback(j, r):
            nonlocal distances
            for k, d in enumerate(r):
                distances[i][k][chunk_size*j:chunk_size*j+len(d)] = d
        print("INFO: Working on trajectory: {}".format(traj_info.path), file=stderr)

        oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    return distances

def cli_parser(prog="distance.py"):
    #handle commandline arguments
    #this program has no positional arguments, only flags
    parser = argparse.ArgumentParser(prog = prog, description="Finds the ensemble of distances between any two particles in the system")
    parser.add_argument('-i', '--input', metavar='input', nargs='+', action='append', help='A trajectory, and a list of particle pairs to compare.  Can call -i multiple times to plot multiple datasets.')
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The name to save the graph file to')
    parser.add_argument('-f', '--format', metavar='<histogram/trajectory/both>', nargs=1, help='Output format for the graphs.  Defaults to histogram.  Options are \"histogram\", \"trajectory\", and \"both\"')
    parser.add_argument('-d', '--data', metavar='data_file', nargs=1, help='If set, the output for the graphs will be dropped as a json to this filename for loading in oxView or your own scripts')
    parser.add_argument('-n', '--names', metavar='names', nargs='+', action='append', help='Names of the data series.  Will default to particle ids if not provided')
    parser.add_argument('-p', '--parallel', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-c', metavar='cluster', dest='cluster', action='store_const', const=True, default=False, help="Run the clusterer on each configuration's distance?")
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python", "matplotlib", "numpy"])

    #-i requires 4 or more arguments, the topology file of the structure, the trajectory to analyze, and any number of particle pairs to compute the distance between.
    try:
        trajectories = [i[0] for i in args.input]
        p1ss = [i[1::2] for i in args.input]
        p2ss = [i[2::2] for i in args.input]
        p1ss = [[int(j) for j in i] for i in p1ss]
        p2ss = [[int(j) for j in i] for i in p2ss]

    except Exception as e:
        print("ERROR:", e)
        parser.print_help()
        exit(1)
    
    #get number of distances to calculate
    n_dists = sum([len(l) for l in p1ss])

    #Make sure that the input is correctly formatted
    if len(p1ss) != len(p2ss):
        raise RuntimeError(" bad input arguments\nPlease supply an even number of particles")

    # Get metadata on the inputs
    top_infos = []
    traj_infos = []
    for traj in trajectories:
        top_info, traj_info = describe(None, traj)
        top_infos.append(top_info)
        traj_infos.append(traj_info)

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else: 
        print("INFO: No outfile name provided, defaulting to \"distance.png\"", file=stderr)
        outfile = "distance.png"

    #-f defines which type of graph to produce
    hist = False
    lineplt = False
    if args.format:
        if "histogram" in args.format:
            hist = True
        if "trajectory" in args.format:
            lineplt = True
        if "both" in args.format:
            hist = True
            lineplt = True
        if hist == lineplt == False:
            raise RuntimeError("Unrecognized graph format\nAccepted formats are \"histogram\", \"trajectory\", and \"both\"")
    else:
        print("INFO: No graph format specified, defaulting to histogram", file=stderr)
        hist = True

    #-c makes it run the clusterer on the output
    cluster = args.cluster

    # -p sets the number of cpus to use
    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    distances = distance(traj_infos, top_infos, p1ss, p2ss, ncpus)

    # -n sets the names of the data series
    if args.names:
        names = args.names[0]
        if len(names) < n_dists:
            print("WARNING: Names list too short.  There are {} items in names and {} distances were calculated.  Will pad with particle IDs".format(len(names), n_dists), file=stderr)
            for i in range(len(names), len(distances)):
                names.append("{}-{}".format([j for sl in p1ss for j in sl][i], [j for sl in p2ss for j in sl][i]))
        if len(names) > n_dists:
            print("WARNING: Names list too long. There are {} items in names and {} distances were calculated.  Truncating to be the same as distances".format(len(names), n_dists), file=stderr)
            names = names[:n_dists]

    else:
        print("INFO: Defaulting to particle IDs as data series names", file=stderr)
        names = ["{}-{}".format(p1, p2) for p1, p2 in zip([i for sl in p1ss for i in sl], [i for sl in p2ss for i in sl])]
    
    # -d will dump the distances as json files for loading with the trajectories in oxView
    if args.data:
        from json import dump
        if len(trajectories) > 1:
            print("INFO: distance lists from separate trajectories are printed to separate files for oxView compatibility.  Trajectory numbers will be appended to your provided data file name.", file=stderr)
            file_names = ["{}_{}.json".format(args.data[0].strip('.json'), i) for i,_ in enumerate(trajectories)]
        else:
            file_names = [args.data[0].strip('.json')+'.json']
        names_by_traj = [['{}-{}'.format(p1, p2) for p1, p2 in zip(p1l, p2l)] for p1l, p2l in zip(p1ss, p2ss)]
        
        for file_name, ns, dist_list in zip(file_names, names_by_traj, distances):
            obj = {}
            for n, d in zip(ns, dist_list):
                obj[n] = d        
            with open(file_name, 'w+') as f:
                print("INFO: writing data to {}.  This can be opened in oxView using the Order parameter selector".format(file_name), file=stderr)
                dump(obj, f)

    #convert the distance list into numpy arrays because they're easier to work with
    for i, l in enumerate(distances):
        distances[i] = np.array(l)
    
    means = [np.mean(i, axis=1) for i in distances]
    medians = [np.median(i, axis=1) for i in distances]
    stdevs = [np.std(i, axis=1) for i in distances]

    #get some min/max values to make the plots pretty
    lower = min((l.min() for l in distances))
    upper = max((l.max() for l in distances))

    #those horrific list comprehensions unpack lists of lists into a single list
    print("input:\t", end='')
    [print("{}-{}\t".format(p1, p2), end='') for p1, p2 in zip([i for sl in p1ss for i in sl], [i for sl in p2ss for i in sl])]
    print("")

    print("name:\t", end='')
    [print("{}\t".format(t), end='') for t in names[:n_dists]]
    print("")

    print("mean:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in [i for sl in means for i in sl]]
    print("")

    print("stdev:\t", end='')
    [print("{:.2f}\t".format(s), end='') for s in [i for sl in stdevs for i in sl]]
    print("")

    print("median:\t", end='')
    [print("{:.2f}\t".format(m), end='') for m in [i for sl in medians for i in sl]]
    print("")

    #make a histogram
    if hist == True:
        if lineplt == True:
            #if making two plots, automatically append the plot type to the output file name
            out = outfile[:outfile.find(".")]+"_hist"+outfile[outfile.find("."):]
        else:
            out = outfile
        bins = np.linspace(np.floor(lower-(lower*0.1)), np.ceil(upper+(upper*0.1)), 60)
        graph_count = 0
        for traj_set in distances:
            for dist_list in traj_set:
                a = plt.hist(dist_list, bins, weights=np.ones(len(dist_list)) / len(dist_list),  alpha=0.5, histtype=u'stepfilled', edgecolor='k', label=names[graph_count])
                graph_count += 1
        plt.xlabel("Distance (nm)")
        plt.ylabel("Normalized frequency")
        plt.legend()
        #plt.show()
        plt.tight_layout()
        print("INFO: Writing histogram to file {}".format(out), file=stderr)
        plt.savefig("{}".format(out))

    #make a trajectory plot
    if lineplt == True:
        if hist == True:
            #clear the histogram plot
            plt.clf()
            #if making two plots, automatically append the plot type to the output file name
            out = outfile[:outfile.find(".")]+"_traj"+outfile[outfile.find("."):]
        else:
            out = outfile
        graph_count = 0
        for traj_set in distances:
            for dist_list in traj_set:
                a = plt.plot(dist_list, alpha=0.5, label=names[graph_count])
                graph_count += 1
        plt.xlabel("Simulation Steps")
        plt.ylabel("Distance (nm)")
        plt.legend()
        #plt.show()
        plt.tight_layout()
        print("INFO: Writing trajectory plot to file {}".format(out), file=stderr)
        plt.savefig("{}".format(out))

    if cluster == True:
        if not all([x == trajectories[0] for x in trajectories]):
            raise RuntimeError("Clustering can only be run on a single trajectory")

        from oxDNA_analysis_tools.clustering import perform_DBSCAN

        labs = perform_DBSCAN(traj_infos[0], top_infos[0], distances[0].T, "euclidean", 12, 8)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()