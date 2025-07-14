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
# This function needs to have exactly three arguments: an class/dict/tuple containing the needed arguments, the number of configurations per chunk, and the current chunk ID
# I name my example variables stupid things to make it clear to you that unlike some of the function names, variable names are arbitrary
# Make your variables descriptive of what they actually are, I see too many uses of boilerplate code which re-use example variable names.
def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    pirate_stuff = ctx.arrrrrg
    return(np.arange(chunk_size)*pirate_stuff)

# All oat scripts have a function with the same name as the script for import into other scripts
# This is convention, if you pull-request something into oat, I will ask for this.
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
    # You only need this if you're going to return a value-per-configuration.  
    # If you're aggregating in the callback, it's unnecessary.
    chunk_size = get_chunk_size()

    # The callback function takes data from each process as it finishes and puts it into some data structure
    # In this example, the parallelized function is returning one value per configuration
    # You can also return multiple values or an array or an object
    # You just need to modify the callback function to handle whatever data type you're returning.
    output = np.zeros(traj_info.nconfs) # pre-allocate the memory for your result
    def callback(i, r): # i is the number of the chunk, r is the return value of the parallelized compute function
        nonlocal output # You need this to access the variable outside the local scope
        output[i*chunk_size:i*chunk_size+len(r)] = r # This is how you can insert per-configuration values into a pre-allocated array

    # Call <compute> with args <ctx> <ncpus> times to process <nconfs> things then package the results with <callback>
    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    return output

# This is what gets picked up by the cli documentation builder
# The function name is not arbitrary, the documentation looks for 'cli_parser' specifically.
def cli_parser(prog="skeleton"):
    parser = argparse.ArgumentParser(prog = prog, description="One sentence description of your program")
    parser.add_argument('trajectory', type=str, help='The trajectory file you wish to analyze')
    parser.add_argument('-f', '--flag1', type=int, help='An optional flag taking an int from the CLI')
    parser.add_argument('-p', '--parallel', metavar='num_cpus', type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-o', '--output', metavar='output_file', help='The filename to save the output to')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

# All main does is handle i/o
# This function name is also not optional, the command line interface looks for a 'main' function with no arguments.
def main():
    # Get arguments from the CLI
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    # Run system checks
    # I check all imported libraries on all scripts
    # set_quiet(True) will turn off the 'INFO:...' prints
    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    # Parse the CLI input
    # Here you can set default values and read files
    # You can also set defaults in cli_parser, but if you want to do anything more complicated, it needs to b ehere.
    # After parsing the input, you should have all the variables ready to call the core function
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
    # We do this in a separate function so that it's also importable into Python scripts
    # This gives users flexibility to use OAT as either a Python library or a set of CLI tools.
    out = skeleton(top_info, traj_info, optional_argument=flag1, ncpus=ncpus)

    # Do something with the output
    # You might want to make this a separate function
    with open(outfile, 'w+') as f:
        f.write(', '.join([str(o) for o in out]))
        log(f"Wrote output to file {outfile}")

    # We generally end scripts with this
    print("--- %s seconds ---" % (time.time() - start_time))

# This lets the Python file be both importable and callable
# Don't add extra stuff to this block unless you know what you're doing
if __name__ == '__main__':
    main()