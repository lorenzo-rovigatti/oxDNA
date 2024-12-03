# This script is a skeleton for parallelized trajectory analysis
import time
start_time = time.time()
import argparse
import numpy as np
from os import path
from collections import namedtuple
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs
from oxDNA_analysis_tools.UTILS.data_structures import TrajInfo, TopInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser, get_chunk_size

# This could also be a dict or a class instance.  We just like namedtuple
ComputeContext = namedtuple("ComputeContext", ["top_info", "traj_info", "arrrrrg"])

# The parallelized function which actually does hard things
def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    pirate_stuff = ctx.arrrrrg
    return(np.arange(chunk_size)*pirate_stuff)

# All oat scripts have a function with the same name as the script for import into other scripts
def skeleton(top_info:TopInfo, traj_info:TrajInfo, optional_argument:int=0, ncpus:int=1) -> np.ndarray:
    """
    This is the docstring that gets picked up by the api documentation builder.  You need to add this script to oxDNA/docs/source/oat/{api.md,cli.md}

    Parameters:
        top_info (TopInfo): Information about the topology
        traj_info (TrajInfo): Information about the trajectory
        optional_argument (int): Tell the documentation about your great ideas (default=0)
        ncpus (int) : (optional) How many cpus to parallelize the operation. default=1

    Returns:
        np.ndarray: In this case just counting up to 2*chunk_size by 2s
    """
    # The ctx is the arguments for the parallelized function
    ctx = ComputeContext(top_info, traj_info, optional_argument)

    # Figure out how much stuff each process is working on (oat_multiprocessor also does this itself, so don't change it here)
    chunk_size = get_chunk_size()

    # The callback function takes data from each process as it finishes and puts it into some data structure
    # Here we assume that the parallelized function is returning one value per configuration
    output = np.zeros(traj_info.nconfs)
    def callback(i, r): # i is the number of the chunk, r is the data that comes back from processing the chunk
        nonlocal output
        output[i*chunk_size:i*chunk_size+len(r)] = r 

    # Call <compute> with args <ctx> <ncpus> times to process <nconfs> things then package the results with <callback>
    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    return output

# This is what gets picked up by the cli documentation builder
def cli_parser(prog="program_name"):
    parser = argparse.ArgumentParser(prog = prog, description="One sentence description of your program")
    parser.add_argument('trajectory', type=str, help='The trajectory file you wish to analyze')
    parser.add_argument('-f', '--flag1', type=int, help='An optional flag taking an int from the CLI')
    parser.add_argument('-p', '--parallel', metavar='num_cpus', type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', help='The filename to save the output to')
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
    if args.flag1:
        flag1 = args.flag1
    else:
        flag1 = 2

    # -p sets the number of cores to use.  Default is 1.
    ncpus = args.parallel if args.parallel else 1

    # -o names the output file
    if args.output:
        outfile = args.output
    else:
        outfile = "out.txt"
        log(f"No outfile name provided, defaulting to \"{outfile}\"")

    # Actually process data
    out = skeleton(top_info, traj_info, optional_argument=flag1, ncpus=ncpus)

    # Do something more complicated than this for the output
    with open(outfile, 'w+') as f:
        f.write(', '.join([str(o) for o in out]))
        log(f"Wrote output to file {outfile}")

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()