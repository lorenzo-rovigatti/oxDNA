import time
start_time = time.time()
import argparse
import numpy as np
import matplotlib.pyplot as plt
from os import path
from json import dump
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.constants import FIG_DPI
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, inbox
from oxDNA_analysis_tools.UTILS.data_structures import TrajInfo, TopInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser, get_chunk_size

ComputeContext = namedtuple("ComputeContext", ["top_info", "traj_info", "indexes"])

# The parallelized function which actually does hard things
# In oxDNA, all particles have the same mass, so we can exclude mass from the calculation
def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    poses = np.array([inbox(c).positions for c in confs])
    poses = poses[:, ctx.indexes, :]
    n_particles = len(ctx.indexes)
    com = np.mean(poses, axis=1, keepdims=True)
    dists = np.linalg.norm(poses-com, axis=2)
    sq_dists = np.power(dists, 2)
    msd = np.sum(sq_dists, axis=1) / n_particles
    rg = np.sqrt(msd) * 0.8518  # Convert from simulation units to nm
    return(rg)

def rg(top_info:TopInfo, traj_info:TrajInfo, indexes:np.ndarray=np.array([]), ncpus:int=1) -> np.ndarray:
    """
    Compute radius of gyration for the provided trajectory

    Parameters:
        top_info (TopInfo): Information about the topology
        traj_info (TrajInfo): Information about the trajectory
        ncpus (int) : (optional) How many cpus to parallelize the operation. default=1

    Returns:
        np.ndarray: Radius of gyration at each configuration in nm
    """
    if indexes.size == 0:
        indexes = np.arange(top_info.nbases)

    # The ctx is the arguments for the parallelized function
    ctx = ComputeContext(top_info, traj_info, indexes)

    # Figure out how much stuff each process is working on
    chunk_size = get_chunk_size()

    # Take the output from each chunk and 
    output = np.zeros(traj_info.nconfs)
    def callback(i, r): # i is the number of the chunk, r is the data that comes back from processing the chunk
        nonlocal output
        output[i*chunk_size:i*chunk_size+len(r)] = r 

    # Call <compute> with args <ctx> <ncpus> times to process <nconfs> things then package the results with <callback>
    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    return output

def make_rg_plot(rg_values:np.ndarray, out:str):
    """
    Make a plot of radius of gyration over time

    Parameters:
        rg_values (np.ndarray): Radius of gyration values
        out (str): Path to save the figure
    """
    fig, ax = plt.subplots()
    ax.plot(rg_values, linewidth=0.5)
    ax.set_xlabel("Step")
    ax.set_ylabel("Radius of gyration (nm)")
    fig.tight_layout()
    fig.savefig(out, dpi=FIG_DPI)
    plt.close(fig)
    log(f"Wrote RG plot to {out}")

# This is what gets picked up by the cli documentation builder
def cli_parser(prog="program_name"):
    parser = argparse.ArgumentParser(prog = prog, description="Calculate radius of gyration over a trajectory")
    parser.add_argument('trajectory', type=str, help='The trajectory file you wish to analyze')
    parser.add_argument('-p', '--parallel', metavar='num_cpus', type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-i', '--index', metavar='index_file', dest='index_file', nargs=1, help='Compute the RG with only a subset of particles from a space-separated list in the provided file')
    parser.add_argument('-o', '--output', metavar='output_file', help='The filename to save the output plot to')
    parser.add_argument('-d', '--data', metavar='data_file', help='The filename to save the RG over time data to')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

# All main does is handle i/o
def main():
    # Get arguments from the CLI
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    #run system checks
    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    # Parse CLI input
    traj = args.trajectory
    top_info, traj_info = describe(None, traj)

    # -p sets the number of cores to use.  Default is 1.
    ncpus = args.parallel if args.parallel else 1

    # -i will use only a subset of nucleotides for alignment.
    # The index file is a space-separated list of particle IDs
    if args.index_file:
        index_file = args.index_file[0]
        with open(index_file, 'r') as f:
            indexes = f.readline().split()
            try:
                indexes = np.array([int(i) for i in indexes])
            except:
                raise RuntimeError("The index file must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button")
    else: 
        indexes = np.array([])

    # -d names the output file
    if args.data:
        datafile = args.data.strip()
        if not datafile.split(".")[-1] == 'json':
            datafile += ".json"
    else:
        datafile = "rg.json"
        log(f"No data file name provided, defaulting to {datafile}")

    # -o names the output plot file
    if args.output:
        plotfile = args.output.strip()
        if not plotfile.split(".")[-1] in ['png', 'jpg', 'jpeg', 'pdf','svg']:
            plotfile += ".png"
    else:
        plotfile = "rg.png"
        log(f"No plot file name provided, defaulting to {plotfile}")

    # Actually process data
    out = rg(top_info, traj_info, indexes=indexes, ncpus=ncpus)

    # Drop output as an oxView OP file
    with open(datafile, 'w+') as f:
        out_obj = {"rg" : [o for o in out]}
        dump(out_obj, f)
        log(f"Wrote oxView overlay file {datafile}")

    # Make a plot
    make_rg_plot(out, plotfile)

    print(f"Mean Radius of Gyration: {np.mean(out):.3f} nm")

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()