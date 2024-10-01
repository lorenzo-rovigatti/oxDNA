import os
import numpy as np
import argparse
import matplotlib.pyplot as plt
from sys import stderr
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.config import check
from oxDNA_analysis_tools.distance import vectorized_min_image

from time import time
start_time = time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info"])

def compute(ctx:ComputeContext, chunk_size:int,  chunk_id:int) -> np.ndarray:
    """
    Computes the average distance between every pair of nucleotides and creates a matrix of these distances.

    Parameters:
        ctx (ComputeContext): A named tuple containing trajectory info, topology info, and the number of configurations to process.
        chunk_id (int): The id of the chunk to process.
    
    Returns:
        distances (numpy.array): A NxN matrix containing pairwise distances between every pair of nucleotides.
    """
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)

    np_poses = np.asarray([c.positions for c in confs])
    distances = np.zeros((ctx.top_info.nbases, ctx.top_info.nbases))
    for c in np_poses:
        distances += vectorized_min_image(c, c, confs[0].box[0])

    return distances

def contact_map(traj_info:TrajInfo, top_info:TopInfo, ncpus=1) -> np.ndarray: 
    """
        Calculates the average contact map for a trajectory.

        Parameters:
            traj_info (TrajInfo): Information about the trajectory.
            top_info (TopInfo): Information about the topology.
            ncpus (int): (optional) The number of cores to use.
    """
    ctx = ComputeContext(traj_info, top_info)

    distances = np.zeros((top_info.nbases, top_info.nbases))
    def callback(i, r):
        nonlocal distances
        distances += r

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    # Normalize the distances and convert to nm
    distances /= traj_info.nconfs
    distances *= 0.8518

    return distances

def cli_parser(prog="contact_map.py"):
    #get commandline arguments
    parser = argparse.ArgumentParser(prog = prog, description="Calculate and display the contact map for a structure")
    parser.add_argument('trajectory', type=str, nargs=1, help="The file containing the configurations of which the contact map is needed")
    parser.add_argument('-g', metavar='graph', dest='graph', nargs=1, type=str, help='Filename for the plot')
    parser.add_argument('-d', metavar='data', dest='data', nargs=1, help='The name of the file to save the contact map as a pickle.')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    return parser

def main():
    check(["python", "numpy", "matplotlib"])

    # Get arguments and file metadata
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()
    traj = args.trajectory[0]
    top_info, traj_info = describe(None, traj)

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    distances = contact_map(traj_info, top_info, ncpus)

    # Plot the contact map
    if args.graph:
        graph_name = args.graph[0]
    else:
        print("INFO: No graph name provided, defaulting to 'contact_map.png'", file=stderr)
        graph_name = "contact_map.png"

    fig, ax = plt.subplots()
    a = ax.imshow(distances, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
    ylabel="nucleotide id",
    xlabel="nucleotide id")
    b = fig.colorbar(a, ax=ax)
    b.set_label("distance (nm)", rotation = 270, labelpad=15)
    plt.tight_layout()
    print("INFO: Saving contact map to '{}'".format(graph_name), file=stderr)
    plt.savefig(graph_name)

    # Save the contact map as a pickle
    if args.data:
        data_name = args.data[0]
    else:
        print("INFO: No data name provided, defaulting to 'contact_map.pkl'", file=stderr)
        data_name = "contact_map.pkl"
    print("INFO: Saving contact map to '{}'".format(data_name), file=stderr)
    np.save(data_name, distances)

    print("--- %s seconds ---" % (time() - start_time))

if __name__ == "__main__":
    main()