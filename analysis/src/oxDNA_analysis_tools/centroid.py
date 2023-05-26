from sys import stderr
from collections import namedtuple
from typing import Tuple, List
import numpy as np
import argparse
import os
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, inbox, write_conf, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.align import svd_align
import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "ref_coords",
                                              "indexes"])


def compute_centroid(ctx:ComputeContext, chunk_size, chunk_id:int) -> Tuple[np.ndarray, float, int]:
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    confs = [inbox(c) for c in confs]
    np_confs = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])
    centroid_candidate = np.zeros_like(np_confs[0])
    min_RMSD = np.inf
    centroid_id = -1

    for i, c in enumerate(np_confs):
        c[0] -= np.mean(c[0][ctx.indexes], axis=0) #didn't center earlier because you have to center on the indexed particles
        aligned_conf = svd_align(ctx.ref_coords.positions[ctx.indexes], c, ctx.indexes, ref_center=np.zeros(3))[0]
        RMSD = np.sqrt(np.mean(np.linalg.norm(aligned_conf[ctx.indexes] - ctx.ref_coords.positions[ctx.indexes], axis=1)**2))
        if RMSD < min_RMSD:
            min_RMSD = RMSD
            centroid_candidate = c
            centroid_id = i

    t = confs[centroid_id].time

    return (centroid_candidate, min_RMSD, t)

def centroid(traj_info:TrajInfo, top_info:TopInfo, ref_conf:Configuration, indexes:List[int]=[], ncpus=1) -> Tuple[Configuration, float]:
    '''
        Find the configuration in a trajectory closest to a provided reference configuration

        Parameters:
            traj_info (TrajInfo): Object containing information about the trajectory
            top_info (TopInfo): Object containing information about the topology
            ref_conf (Configuration): Object containing the reference configuration
            indexes (List[int]): (optional) Indexes of the particles to be used for alignment
            ncpus (int): (optional) Number of CPUs to use for alignment

        Returns:
            centroid_candidate (Configuration): The configuration with the lowest RMSD to the reference
            min_RMSD (float): The RMSD from the centroid to the reference
    '''
    if indexes == []:
        indexes = list(range(top_info.nbases))

    ref_conf = inbox(ref_conf)
    ref_cms = np.mean(ref_conf.positions[indexes], axis=0)
    ref_conf.positions -= ref_cms

    # create a ComputeContext which defines the problem to pass to the worker processes
    ctx = ComputeContext(
        traj_info, top_info, ref_conf, indexes
    )

    # What do we do with the output from the worker processes?
    min_RMSD = np.inf
    centroid_candidate = Configuration(0, ref_conf.box, np.zeros(3), np.zeros_like(ref_conf.positions), np.zeros_like(ref_conf.positions), np.zeros_like(ref_conf.positions))
    centroid_time = -1
    def callback(i, r):
        nonlocal min_RMSD, centroid_candidate, centroid_time
        centroid, RMSD, t = r
        if RMSD < min_RMSD:
            min_RMSD = RMSD
            centroid_time = t
            centroid_candidate.positions = centroid[0]
            centroid_candidate.a1s = centroid[1]
            centroid_candidate.a3s = centroid[2]

    oat_multiprocesser(traj_info.nconfs, ncpus, compute_centroid, callback, ctx)
    
    min_RMSD *= 0.8518
    centroid_candidate.time = centroid_time

    return centroid_candidate, min_RMSD

def cli_parser(prog = "centroid.py"):
    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = prog, description="Find the configuration in a trajectory closest to a provided reference configuration")
    parser.add_argument('reference_structure', type=str, nargs=1, help="The reference structure to search against")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the centroid to')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Alignment and RMSD based on a subset of particles given in a space-separated list in the provided file')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    #system check
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    #Get file paths
    ref = args.reference_structure[0].strip()
    traj = args.trajectory[0].strip()
    _, ref_info = describe(None, ref)
    top_info, traj_info = describe(None, traj)

    # -i comes with a list of particles indices representing a subset to compute the mean against.
    # Get the index list which is a space-separated list of particle ids.
    if args.index_file:
        index_file = args.index_file[0]
        with open(index_file, 'r') as f:
            indexes = f.readline().split()
            try:
                indexes = [int(i) for i in indexes]
            except:
                raise RuntimeError("The index file must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button")
    else:
        indexes = list(range(top_info.nbases))

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # get the mean structure from the file path
    ref_conf = get_confs(top_info, ref_info, 0, 1)[0]
    
    centroid_candidate, min_RMSD = centroid(traj_info, top_info, ref_conf, indexes, ncpus)

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
    else: 
        outfile = "centroid.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    write_conf(outfile, centroid_candidate, include_vel=traj_info.incl_v)
    print("INFO: Wrote centroid to {}".format(outfile), file=stderr)
    print("INFO: Min RMSD: {} nm".format(min_RMSD), file=stderr)
    print("INFO: Centroid time: {}".format(centroid_candidate.time), file=stderr)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()