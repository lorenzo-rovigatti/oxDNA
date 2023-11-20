import argparse
import os
import time
from typing import List, Union
import numpy as np
from sys import stderr
from collections import namedtuple
from random import randrange
from oxDNA_analysis_tools.align import svd_align
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, inbox, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, TopInfo, TrajInfo
start_time = time.time()

# A compute context is a single variable (in this case a namedtuple, but could also be a dict or a dataclass)
# which contains the arguments needed for a parallelized function.
# We use this style to avoid having to pass variable numbers of arguments to the parallelizer.
ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "centered_ref_coords",
                                              "indexes"])

# This is the function which computes the sum of particle positions for a chunk of a trajectory
# This function will be parallelized. 
# All parallelized functions MUST have these same three arguments                                             
def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    """
        Compute the sum of of positions and orientations for a chunk of configurations

        Parameters:
            ctx (ComputeContext): A namedtuple containing file information, the reference conf and the indexes to compute.
            chunk_size (int): The number of confs to compute in a chunk.
            chunk_id (int): The id of the chunk to compute.
    """
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    
    # Because of fix_diffusion, anything that performs alignment must be inboxed first.
    confs = (inbox(c, center=True) for c in confs)
    
    # convert to numpy repr for easier math
    np_coords = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])
    sub_mean = np.zeros(shape=[3,ctx.top_info.nbases,3])
    
    # sum over confs
    for c in np_coords:
        sub_mean += svd_align(ctx.centered_ref_coords, c, ctx.indexes, ref_center=np.zeros(3))
    
    return sub_mean

def mean(traj_info:TrajInfo, top_info:TopInfo, ref_conf:Union[Configuration,None]=None, indexes:List[int]=[], ncpus:int=1) -> Configuration:
    """
        Compute the mean structure of a trajectory.

        Parameters:
            traj_info (TrajInfo): Information about the trajectory
            top_info (TopInfo): Information about the topology
            ref_conf (Configuration): (optional) The reference configuration to align to. If None, a random configuraiton will be used.
            indexes (List[int]): (optional) The indexes of nucleotides included in the alignment. If None, all nucleotides will be aligned.
            ncpus (int): (optional) The number of CPUs to use. If None, 1 CPU will be used.
        
        Returns:
            Configuration (Configuration): The mean structure of the trajectory.
    """

    # Handle case where function was called from another script with incomplete arguments
    if indexes == []:
        indexes = list(range(top_info.nbases))
    if ref_conf == None:
        ref_conf_id = int(randrange(0, traj_info.nconfs))
        ref_conf = get_confs(top_info, traj_info, ref_conf_id, 1)[0]
    
    # alignment requires the ref to be centered at 0
    reference_coords = ref_conf.positions[indexes]
    ref_cms = np.mean(reference_coords, axis=0) # cms prior to centering
    reference_coords = reference_coords - ref_cms

    # The compute context is a single variable containin the arguments needed for the parallelized function. 
    ctx = ComputeContext(
        traj_info, top_info, reference_coords, indexes
    )

    # The callback function is called after each chunk is computed.
    # For computing the mean structure, we just add the sum of particle positions to a running total
    # we'll divide by the number of confs later.
    mean = np.zeros([3, top_info.nbases, 3])
    def callback(i, r):
        nonlocal mean
        mean += r

    # The parallelizer will call the "compute" function with ctx as an argument using "ncpus".
    # The "callback" function is called after each chunk is computed.
    # Chunk size is a global variable which can be set by calling `oat config -n <chunk_size>`
    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    # divide the sum by the number of confs to get the mean 
    mean /= traj_info.nconfs
    pos, a1s, a3s = mean

    # renormalize a1 and a3 because weird things happen if they aren't
    a1s = np.array([v/np.linalg.norm(v) for v in a1s])
    a3s = np.array([v/np.linalg.norm(v) for v in a3s])

    return Configuration(0,ref_conf.box,np.array([0,0,0]), pos, a1s , a3s)

# Sphinx can autogenerate cli documentation if you have a function which returns your parser.
# To get it to show up in the documentation you need to add the script to oxDNA/docs/oat/cli.md
def cli_parser(prog="mean.py"):
    # A standard way to create and parse command line arguments.
    parser = argparse.ArgumentParser(prog = prog, description="Computes the mean structure of a trajectory file")
    parser.add_argument('trajectory', type=str, nargs=1, help='The trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the mean structure to')
    parser.add_argument('-d', '--deviations', metavar='deviation_file', nargs=1, help='Immediately run oat deviations from the output')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-a', '--align', metavar='alignment_configuration', nargs=1, help='The id of the configuration to align to, otherwise random')
    return parser

# All scripts in oat must have a main method with no arguments to work with the command line interface.
def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    # Verify that dependencies are installed and a good version
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    # Get metadata about input files
    traj = args.trajectory[0]
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

    # -a specifies the id of the reference configuration, defaults to random
    # Get the reference configuration and bring it back in the box
    if args.align:
        ref_conf_id = int(args.align[0])
    else:
        ref_conf_id = int(randrange(0, traj_info.nconfs))
    ref_conf = get_confs(top_info, traj_info, ref_conf_id, 1)[0]
    ref_conf = inbox(ref_conf)

    # -p sets the number of cores to use.  Default is 1.
    if args.parallel:
        ncpus = args.parallel[0]
    else:
        ncpus = 1

    # Actually perform the mean computation
    mean_conf = mean(traj_info, top_info, ref_conf, indexes, ncpus)

    #-o names the output file
    if args.output:
        outfile = args.output[0]
    else:
        outfile = "mean.dat"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    # Create the mean configuration from the numpy arrays containing the positions and orientations
    # And write it to the outfile
    write_conf(outfile, mean_conf, include_vel=traj_info.incl_v)
    print("--- %s seconds ---" % (time.time() - start_time))

    # -d runs deviations.py after computing the mean
    if args.deviations:
        from oxDNA_analysis_tools import deviations
        dev_file = args.deviations[0]
        print("INFO: Launching compute_deviations", file=stderr)

        RMSDs, RMSFs = deviations.deviations(traj_info, top_info, mean_conf, indexes, ncpus)
        deviations.output(RMSDs, RMSFs, dev_file, dev_file.split('.')[0]+"_rmsd.png", dev_file.split('.')[0]+"_rmsd_data.json")

# This makes the script executable via the python interpreter since everything is inside the main() function.
if __name__ == '__main__':
    main()
