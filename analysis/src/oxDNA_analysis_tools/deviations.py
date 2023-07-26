import argparse
import os
from typing import List, Tuple
import numpy as np
from sys import stderr
from collections import namedtuple
from json import dumps
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser, get_chunk_size
from oxDNA_analysis_tools.UTILS.RyeReader import describe, inbox, get_confs
from oxDNA_analysis_tools.align import svd_align
import matplotlib.pyplot as plt
import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "mean_coords",
                                              "indexes"])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    confs = (inbox(c, center=True) for c in confs)
    confs = np.asarray([[c.positions, c.a1s, c.a3s] for c in confs])

    SFs = np.empty((len(confs), ctx.top_info.nbases))
    for i, c in enumerate(confs):
        aligned_conf = svd_align(ctx.mean_coords.positions[ctx.indexes], c, ctx.indexes, ref_center=np.zeros(3))[0]
        SFs[i] = np.power(np.linalg.norm(aligned_conf - ctx.mean_coords.positions, axis=1), 2)

    return SFs

def deviations(traj_info:TrajInfo, top_info:TopInfo, mean_conf:Configuration, indexes:List[int]=[], ncpus:int=1) -> Tuple[np.ndarray, np.ndarray]:
    """
        Find the deviations of a trajectory from a mean configuration

        Parameters:
            traj_info (TrajInfo): Information about the trajectory
            top_info (TopInfo): Information about the topology
            mean_conf (Configuration): The mean configuration
            indexes (List[int]): (optional) List of indexes of the particles to be used for alignment
            ncpus (int): (optional) Number of CPUs to use for alignment
        
        Returns:
            RMSDs (np.array): Root mean squared deviation for each configuration in the trajectory
            RMSFs (np.array): Average fluctuation for each particle in the structure
    """
    if indexes == []:
        indexes = list(range(top_info.nbases))

    mean_conf = inbox(mean_conf)
    ref_cms = np.mean(mean_conf.positions[indexes], axis=0)
    mean_conf.positions -= ref_cms

    # create a ComputeContext which defines the problem to pass to the worker processes
    ctx = ComputeContext(
        traj_info, top_info, mean_conf, indexes
    )

    chunk_size = get_chunk_size()
    MFs = np.empty((traj_info.nconfs, top_info.nbases))
    def callback(i, r):
        nonlocal MFs, chunk_size
        MFs[i*chunk_size:i*chunk_size+len(r)] = r

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    # Compute RMSDs and RMSF
    RMSDs = np.sqrt(np.mean(MFs, axis=1)) * 0.8518
    RMSFs = np.sqrt(np.mean(MFs, axis=0)) * 0.8518

    return (RMSDs, RMSFs)

def output(RMSDs:np.ndarray, RMSFs:np.ndarray, outfile:str='devs.json', plot_name:str='rmsd.png', data_file:str='rmsd_op.json'):
    """
        Create RMSF oxView overlay and RMSD plot

        Parameters:
            RMSDs (np.array): Root mean squared deviation for each configuration in the trajectory
            RMSFs (np.array): Average deviation for each particle in the structure
            outfile (str): (optional) Name of the oxView overlay file for the RMSF
            plot_name (str): (optional) Name of the RMSD plot
            data_file (str): (optional) Name of the oxView order parameter file for the RMSD
    """
    # Save the RMSDs and RMSFs to json files
    print("INFO: writing deviations to {}".format(outfile), file=stderr)
    with open(outfile, 'w') as f:
        f.write(
            dumps({
                "RMSF (nm)" : RMSFs.tolist()
            })
        )

    print("INFO: writing RMSDs to oxView order parameter file, {}".format(data_file), file=stderr)
    with open(data_file, 'w') as f:
        f.write(
            dumps({
                "RMSD (nm)" : RMSDs.tolist()
            })
        )

    print("INFO: writing RMSD plot to {}".format(plot_name), file=stderr)
    plt.plot(RMSDs)
    plt.axhline(np.mean(RMSDs), color='red')
    plt.xlabel('Configuration')
    plt.ylabel('RMSD (nm)')
    plt.tight_layout()
    plt.savefig(plot_name)

    return

def cli_parser(prog="deviations.py"):
    #handle commandline arguments
    #the positional arguments for this are: 
    # 1. the mean structure from compute_mean.py in json format
    # 2. the trajectory from which to compute the deviations
    # 3. the topology associated with the trajectory
    parser = argparse.ArgumentParser(prog = prog, description="Compute the RMSD of each nucleotide from the mean structure produced by compute_mean.py")
    parser.add_argument('mean_structure', type=str, nargs=1, help="The mean structure .json file from compute_mean.py")
    parser.add_argument('trajectory', type=str, nargs=1, help='the trajectory file you wish to analyze')
    parser.add_argument('-p', metavar='num_cpus', nargs=1, type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', nargs=1, help='The filename to save the RMSF json file to')
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Compute mean structure of a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-r', metavar='rmsd_plot', dest='rmsd_plot', nargs=1, help='The name of the file to save the RMSD plot to.')
    parser.add_argument('-d', metavar='rmsd_data', dest='rmsd_data', nargs=1, help='The name of the file to save the RMSD data in json format.')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    #system check
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy", "matplotlib"])

    # Get metadata about input files
    mean = args.mean_structure[0]
    traj = args.trajectory[0]
    _, mean_info = describe(None, mean)
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
    mean_conf = get_confs(top_info, mean_info, 0, 1)[0]
    RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, indexes, ncpus)

    #-o names the output file
    if args.output:
        outfile = args.output[0].strip()
        if not outfile.split(".")[-1] == 'json':
            outfile += ".json"
    else: 
        outfile = "devs.json"
        print("INFO: No outfile name provided, defaulting to \"{}\"".format(outfile), file=stderr)

    #-r names the file to print the RMSD plot to
    if args.rmsd_plot:
        plot_name = args.rmsd_plot[0]
    else:
        plot_name = 'rmsd.png'

    # -d names the file to print the RMSD data to
    if args.rmsd_data:
        data_file = args.rmsd_data[0]
    else:
        data_file = 'rmsd_op.json'

    output(RMSDs, RMSFs, outfile, plot_name, data_file)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()