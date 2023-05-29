import argparse
import os
import time
from typing import List, Tuple, Union
import numpy as np
from sys import stderr
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, inbox, conf_to_str

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "centered_ref_coords",
                                              "indexes",
                                              "center"])

def svd_align(ref_coords:np.ndarray, coords:np.ndarray, indexes:List[int], ref_center:np.ndarray=np.array([]), center:bool=True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Single-value decomposition-based alignment of configurations

    Parameters:
        ref_coords (numpy.ndarray): Reference coordinates.  Should be indexed before calling this function.
        coords (numpy.ndarray): Coordinates to be aligned
        indexes (List[int]): Indexes of the atoms to be aligned in the coords array
        ref_center (numpy.ndarray): (optional) The center of mass of the reference configuration. If not provided, it will be calculated (slightly slower for many confss).

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray] A tuple of the aligned coordinates (coords, a1s, a3s) for the given chunk
    """
    if len(ref_center) == 0:
        ref_center = np.mean(ref_coords, axis=0)
    av1 = ref_center
    av2 = np.mean(coords[0][indexes], axis=0)
    coords[0] = coords[0] - av2
    # correlation matrix
    a = np.dot(np.transpose(coords[0][indexes]), ref_coords - av1)
    u, _, vt = np.linalg.svd(a)
    rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    # check if we have found a reflection
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))
    tran = av1
    return  (np.dot(coords[0], rot) + tran,
             np.dot(coords[1], rot),
             np.dot(coords[2], rot))

def compute(ctx:ComputeContext, chunk_size, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    confs = [inbox(c, center=ctx.center) for c in confs]
    # convert to numpy repr
    np_coords = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])

    # align
    for i, c in enumerate(np_coords):
        c[0], c[1], c[2] = svd_align(ctx.centered_ref_coords, c, ctx.indexes, ref_center=np.zeros(3), center=ctx.center)
        confs[i].positions = c[0]
        confs[i].a1s = c[1]
        confs[i].a3s = c[2]
    #return confs
    out = ''.join([conf_to_str(c, include_vel=ctx.traj_info.incl_v) for c in confs])
    return out

def align(traj:str, outfile:str, ncpus:int=1, indexes:List[int]=[], ref_conf:Union[Configuration,None]=None, center:bool=True):
    """
        Align a trajectory to a ref_conf and print the result to a file.

        Parameters:
            traj (str) : The trajectory file name to align
            outfile (str) : The file name to write the aligned trajectory to
            ncpus (int) : (optional) How many cpus to parallelize the operation. default=1
            indexes (List[int]) : (optional) IDs of a subset of particles to consider for the alignment. default=all
            ref_conf (Configuration) : (optional) The configuration to align to. default=first conf

        Writes the aligned configuration to outfile
    """
    
    top_info, traj_info = describe(None, traj)

    if ref_conf == None:
        #read the first configuration and use it as the reference configuration for the rest
        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]

    if indexes == []:
        indexes = list(range(top_info.nbases))

    ref_conf = inbox(ref_conf) # Don't need to center now because we're going to after indexing anyway.

    # alignment requires the ref to be centered at 0.  Inboxing did not take the indexing into account.
    reference_coords = ref_conf.positions[indexes]
    ref_cms = np.mean(reference_coords, axis=0) # cms prior to centering
    reference_coords = reference_coords - ref_cms

    # Create a ComputeContext which defines the problem to pass to the worker processes 
    ctx = ComputeContext(
        traj_info, top_info, reference_coords, indexes, center
    )

    with open(outfile, 'w+') as f:
        def callback(i, r):
            nonlocal f
            f.write(r)

        oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)
    
    print(f"INFO: Wrote aligned trajectory to {outfile}", file=stderr)
    return

def cli_parser(prog="align.py"):
    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = prog, description="Aligns each frame in a trajectory to the first frame")    
    parser.add_argument('traj', type=str, nargs=1, help="The trajectory file to align")
    parser.add_argument('outfile', type=str, nargs=1, help='The name of the new trajectory file to write out')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-r', metavar='reference_structure', dest='reference_structure', nargs=1, help="Align to a provided configuration instead of the first frame.")
    parser.add_argument('-c', metavar='no_center', dest='no_center', action='store_const', const=True, default=False, help="Don't center the output.  Can avoid errors caused by small boxes.")
    return parser

def main():
    start_time = time.time()
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    #Parse command line arguments
    traj_file = args.traj[0]
    outfile = args.outfile[0]
    top_info, traj_info = describe(None, traj_file)

    #-r will make it align to a provided .dat file instead of the first configuration
    if args.reference_structure:
        #read reference configuration
        _, ref_info = describe(None, args.reference_structure[0])
        ref_conf = get_confs(top_info, ref_info, 0, 1)[0]
    else:
        ref_conf = None

    #-i will make it only run on a subset of nucleotides.
    #The index file is a space-separated list of particle IDs
    if args.index_file:
        index_file = args.index_file[0]
        with open(index_file, 'r') as f:
            indexes = f.readline().split()
            try:
                indexes = [int(i) for i in indexes]
            except:
                raise RuntimeError("The index file must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button")
    else: 
        indexes = []

    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    center = not args.no_center

    align(traj=traj_file, outfile=outfile, ncpus=ncpus, indexes=indexes, ref_conf=ref_conf, center=center)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()