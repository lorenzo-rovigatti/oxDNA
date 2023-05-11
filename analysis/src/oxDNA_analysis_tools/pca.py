import argparse
from typing import Tuple
import numpy as np 
import matplotlib.pyplot as plt
from sys import exit, stderr
from json import dumps
from warnings import catch_warnings, simplefilter
from os import path
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser, get_chunk_size
from oxDNA_analysis_tools.config import check
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, inbox

import time
start_time = time.time()

ComputeContext_cov = namedtuple("ComputeContext_cov",["traj_info",
                                                      "top_info",
                                                      "centered_ref_coords"])

ComputeContext_map = namedtuple("ComputeContext_map",["traj_info",
                                                      "top_info",
                                                      "centered_ref_coords",
                                                      "components"])

def align_positions(centered_ref_coords:np.ndarray, coords:np.ndarray) -> np.ndarray:
    """
    Single-value decomposition-based alignment of configurations

    This one only considers positions, unlike the one in align which also handles the a vectors

    Parameters
        centered_ref_coords (np.array): reference coordinates, centered on [0, 0, 0]
        coords (np.array): coordinates to be aligned

    Returns
        (np.array) : Aligned coordinates for the given conf
    """
    # center on centroid
    av1, reference_coords = np.zeros(3), centered_ref_coords.copy()
    av2 = np.mean(coords, axis=0)
    coords = coords - av2
    # correlation matrix
    a = np.dot(np.transpose(coords), reference_coords)
    u, _, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # check if we have found a reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    tran = av1 - np.dot(av2, rot)
    return  np.dot(coords, rot) + tran

def make_heatmap(covariance:np.ndarray):
    """
        Produces a heatmat plot of the covariance between every particle

        Parameters:
            covariance (numpy.array): The covariance matrix of particle positions

        Displays:
            A matplotlib imshow plot of the covariance
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    a = ax.imshow(covariance, cmap='viridis', origin='lower')
    ax.set(title = "interaction network",
       ylabel="nucleotide id*3",
       xlabel="nucleotide id*3")
    b = fig.colorbar(a, ax=ax)
    b.set_label("covariance", rotation = 270)
    plt.tight_layout()
    plt.savefig("heatmap.png")


def compute_cov(ctx:ComputeContext_cov, chunk_size:int, chunk_id:int):
    # get a chunk of confs and convert the positions to numpy arrays
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    covariation_matrix = np.zeros((ctx.top_info.nbases*3, ctx.top_info.nbases*3))
    for c in confs:
        c = inbox(c, center=True)
        c.positions = align_positions(ctx.centered_ref_coords, c.positions)
        difference_matrix = (c.positions - ctx.centered_ref_coords).flatten()
        covariation_matrix += np.einsum('i,j -> ij', difference_matrix, difference_matrix)

    return covariation_matrix

def map_confs_to_pcs(ctx:ComputeContext_map, chunk_size:int, chunk_id:int):
    """
    Transforms each configuration in a trajectory into a point in principal component space

    Parameters:
        ctx (ComputeContext_map) : A compute context which contains trajectory metadata, the align conf and the components
        cunk_id (int) : The id of the current chunk

    Returns:
        coordinates (numpy.array): The positions of each frame of the trajectory in principal component space.
    """

    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    coordinates = np.zeros((len(confs), ctx.top_info.nbases*3))
    for i, c in enumerate(confs):
        c = inbox(c, center=True)
        c.positions = align_positions(ctx.centered_ref_coords, c.positions)
        coordinates[i] = np.dot(ctx.components, c.positions.flatten())
    
    return coordinates

def pca(traj_info:TrajInfo, top_info:TopInfo, mean_conf:Configuration, ncpus:int=1) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
        Performs a PCA on a trajectory

        Parameters:
            traj_info (TrajInfo) : Information about the trajectory
            top_info (TopInfo) : Information about the topology
            mean_conf (Configuration) : The mean structure of the trajectory
            ncpus (int) : (optional) The number of CPUs to use for the computation

        Returns:
            (np.ndarray, np.ndarray, np.ndarray) : The structures mapped to coordinate space, the eigenvalues and the eigenvectors

    """
    
    # Create a ComputeContext which defines the problem to pass to the worker processes 
    ctx = ComputeContext_cov(
        traj_info, top_info, mean_conf.positions)

    covariation_matrix = np.zeros((ctx.top_info.nbases*3, ctx.top_info.nbases*3))
    def callback(i, r):
        nonlocal covariation_matrix
        covariation_matrix += r

    oat_multiprocesser(traj_info.nconfs, ncpus, compute_cov, callback, ctx)

    covariation_matrix /= (traj_info.nconfs-1)

    #now that we have the covatiation matrix we're going to use eigendecomposition to get the principal components.
    #make_heatmap(covariance)
    print("INFO: calculating eigenvectors", file=stderr)
    evalues, evectors = np.linalg.eig(covariation_matrix) #these eigenvalues are already sorted
    evectors = evectors.T #vectors come out as the columns of the array
    print("INFO: eigenvectors calculated", file=stderr)

    print("INFO: Saving scree plot to scree.png", file=stderr)
    plt.scatter(range(0, len(evalues)), evalues, s=25)
    plt.xlabel("component")
    plt.ylabel("eigenvalue")
    plt.tight_layout()
    plt.savefig("scree.png")

    total = sum(evalues)
    running = 0
    i = 0
    while running < 0.9:
        running += (evalues[i] / total)
        i += 1
    print("90% of the variance is found in the first {} components".format(i))

    ctx = ComputeContext_map(traj_info, top_info, mean_conf.positions, evectors)

    chunk_size = get_chunk_size()
    coordinates = np.zeros((traj_info.nconfs, ctx.top_info.nbases*3))
    def callback2(i, r):
        nonlocal coordinates, chunk_size
        coordinates[i*chunk_size:(i*chunk_size)+len(r)] = r

    oat_multiprocesser(traj_info.nconfs, ncpus, map_confs_to_pcs, callback2, ctx)

    return (coordinates, evalues, evectors)

def cli_parser(prog="pca.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Calculates a principal component analysis of nucleotide deviations over a trajectory")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('meanfile', type=str, nargs=1, help='The mean structure from oat mean')
    parser.add_argument('outfile', type=str, nargs=1, help='the name of the oxView .json overlay file where the PCs will be written')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")    
    parser.add_argument('-c', metavar='cluster', dest='cluster', action='store_const', const=True, default=False, help="Run the clusterer on each configuration's position in PCA space?")
    parser.add_argument('-n', metavar='num_components', nargs=1, type=int, dest='N', help="(optional) Print the first N components as oxView overlay files (defaults to 1)")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    check(["python", "numpy"])

    traj_file = args.trajectory[0]
    mean_file = args.meanfile[0]
    outfile = args.outfile[0]
    parallel = args.parallel

    # Get inputfile metadata
    top_info, traj_info = describe(None, traj_file)
    _, mean_info = describe(None, mean_file)

    # Get the mean structure and center it
    align_conf = get_confs(top_info, mean_info, 0, 1)[0]
    cms = np.mean(align_conf.positions, axis=0)
    align_conf.positions -= cms

    # -p sets the number of CPUs to use
    if parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    #-c makes it run the clusterer on the output
    cluster = args.cluster

    coordinates, evalues, evectors = pca(traj_info, top_info, align_conf, ncpus)

    #make a quick plot from the first three components
    print("INFO: Creating coordinate plot from first three eigenvectors.  Saving to coordinates.png", file=stderr)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(coordinates[:,0], coordinates[:,1], coordinates[:,2], c='g', s=25)
    plt.tight_layout()
    plt.savefig("coordinates2.png")

    #Create an oxView overlays for the first N components
    if args.N:
        N = args.N[0]
    else:
        N = 1
    prep_pos_for_json = lambda conf: list(
                        list(p) for p in conf
                        )
    print("INFO: Change the number of eigenvalues to sum and display by modifying the N variable in the script.  Current value: {}".format(N), file=stderr)
    for i in range(0, N): #how many eigenvalues do you want?
        f = outfile.strip(".json")+str(i)+".json"
        out = np.sqrt(evalues[i])*evectors[i]

        with catch_warnings(): #this produces an annoying warning about casting complex values to real values that is not relevant
            simplefilter("ignore")
            output_vectors = out.reshape(int(out.shape[0]/3), 3).astype(float)
        
        with open(f, "w+") as file:
            file.write(dumps({
                "pca" : prep_pos_for_json(output_vectors)
            }))

    print("--- %s seconds ---" % (time.time() - start_time))

    #If we're running clustering, feed the linear terms into the clusterer
    if cluster:
    #    print("INFO: Mapping configurations to component space...", file=stderr)
#
    #    #If you want to cluster on only some of the components, uncomment this
    #    #coordinates = coordinates[:,0:3]
#
        from oxDNA_analysis_tools.clustering import perform_DBSCAN
        labs = perform_DBSCAN(traj_info, top_info, coordinates, "euclidean", 12, 8)

if __name__ == '__main__':
    main()    