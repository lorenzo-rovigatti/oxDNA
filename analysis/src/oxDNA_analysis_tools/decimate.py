import argparse
import time
import os
from collections import namedtuple
from copy import deepcopy
from sys import stderr
from typing import Union
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, conf_to_str

ComputeContext = namedtuple("ComputeContext", [
    "traj_info",
    "top_info"
])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
        
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_size*chunk_id, chunk_size)

    out = ''.join([conf_to_str(c,  include_vel=ctx.traj_info.incl_v) for c in confs])
    return out


def decimate(traj:str, outfile:str, ncpus:int=1, start:int=0, stop:Union[int,None]=None, stride:int=10):
    """
        Reduce the number of configurations in a trajectory.

        Parameters:
            traj (str) : The trajectory file name to decimate
            outfile (str) : The file name to write the decimates trajectory to
            ncpus (int) : (optional) How many cpus to parallelize the operation. default=1
            start (int) : (optional) Starting configuration for the new trajectory.  Accepts negative indexes.  default=0
            stop (int) : (optional) Process up to this conf (exclusive).  Accepts negative indexes. 
            stride (int) : (optional) Include only every stride-th conf. (default=10)
    """
    top_info, traj_info = describe(None, traj)
    
    #things outside the stop/start/stride bounds just don't exist anymore.
    my_di = deepcopy(traj_info)
    my_di.idxs = traj_info.idxs[start:stop:stride]
    my_di.nconfs = len(my_di.idxs)

    ctx = ComputeContext(
        my_di,
        top_info
    )

    with open(outfile, 'w+') as f:
        def callback(i, r):
            nonlocal f
            f.write(r)

        oat_multiprocesser(my_di.nconfs, ncpus, compute, callback, ctx)

    print(f"INFO: Wrote decimated trajectory to {outfile}", file=stderr)
    return

def cli_parser(prog="decimate.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Creates a smaller trajectory only including start/stop/stride frames from the input.")    
    parser.add_argument('traj', type=str, help="The trajectory file to decimate")
    parser.add_argument('outfile', type=str, help='The name of the new trajectory file to write out')
    parser.add_argument('-p', dest='parallel', default=1, type=int, help="(optional) How many cores to use")
    parser.add_argument('-s', dest='start', default=0, type=int, help='First conf to write to the output file.')
    parser.add_argument('-e', dest='stop', default=None, type=int, help='Process up to this conf (exclusive).  Accepts negative indexes.')
    parser.add_argument('-d', dest='stride', default=10, type=int, help='Write out every this many confs (default=10)')
    return parser

def main():
    start_time  = time.time()
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python"])

    #Parse command line arguments
    traj = args.traj
    outfile = args.outfile
    ncpus = args.parallel
    start = args.start
    stop = args.stop
    stride = args.stride

    decimate(traj=traj, outfile=outfile, ncpus=ncpus, start=start, stop=stop, stride=stride)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()