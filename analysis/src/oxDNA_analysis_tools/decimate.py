import argparse
import time
import os
from collections import namedtuple
from copy import deepcopy
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from typing import Union
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, conf_to_str
import mmap
ComputeContext = namedtuple("ComputeContext", [
    "traj_info",
    "top_info"
])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
        
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_size*chunk_id, chunk_size)

    out = ''.join([conf_to_str(c,  include_vel=ctx.traj_info.incl_v) for c in confs])
    return out


def decimate(traj: str, outfile: str, start: int, stop: int, stride: int):
    """
        Reduce the number of configurations in a trajectory.

        Parameters:
            traj (str) : The trajectory file name to decimate
            outfile (str) : The file name to write the decimates trajectory to
            start (int) : (optional) Starting configuration for the new trajectory.  Accepts negative indexes.  default=0
            stop (int) : (optional) Process up to this conf (exclusive).  Accepts negative indexes. 
            stride (int) : (optional) Include only every stride-th conf. (default=10)

        Writes the trajectory directly to outfile.
    """
    # Describe the trajectory and extract indices
    _, traj_info = describe(None, traj)
    selected_idxs = traj_info.idxs[start:stop:stride]
    
    # Pre-extract offsets and sizes to minimize attribute lookups
    selected_data = [(idx.offset, idx.size) for idx in selected_idxs]
    
    with open(traj, 'rb') as infile, open(outfile, 'wb') as out_file:
        # Memory-map the input file for faster random access
        with mmap.mmap(infile.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            for offset, size in selected_data:
                out_file.write(mm[offset : offset + size])
    log(f"Wrote decimated trajectory to {outfile}")

def cli_parser(prog="decimate.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Creates a smaller trajectory only including start/stop/stride frames from the input.")    
    parser.add_argument('traj', type=str, help="The trajectory file to decimate")
    parser.add_argument('outfile', type=str, help='The name of the new trajectory file to write out')
    parser.add_argument('-s', '--start', dest='start', default=0, type=int, help='First conf to write to the output file.')
    parser.add_argument('-e', '--stop', dest='stop', default=-1, type=int, help='Process up to this conf (exclusive).  Accepts negative indexes.')
    parser.add_argument('-d', '--stride', dest='stride', default=10, type=int, help='Write out every this many confs (default=10)')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

def main():
    start_time  = time.time()
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python"])

    #Parse command line arguments
    traj = args.traj
    outfile = args.outfile
    start = args.start
    stop = args.stop
    stride = args.stride

    decimate(traj=traj, outfile=outfile, start=start, stop=stop, stride=stride)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()