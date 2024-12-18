import argparse
from typing import Union
from os import remove, path
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from collections import namedtuple
from numpy import round, zeros_like
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, conf_to_str
import time
start_time = time.time()

ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "d",
                                              "a"])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    
    for conf in confs:
        if ctx.d is not None: #round positions
            conf.positions = round(conf.positions, ctx.d)
            conf.a1s = round(conf.a1s, ctx.d)
            conf.a3s = round(conf.a3s, ctx.d)
        if ctx.a: #discard a vectors
            conf.a1s = zeros_like(conf.a1s, dtype=int)
            conf.a3s = zeros_like(conf.a1s, dtype=int)

    out = ''.join([conf_to_str(c,  include_vel=False) for c in confs])
    return out

def minify(traj_info:TrajInfo, top_info:TopInfo, out:str, d:Union[int,None]=None, a:bool=False, ncpus=1):
    """
        Make a trajectory smaller by discarding some precision.

        Parameters:
            traj_info (TrajInfo): Information about the trajectory
            top_info (TopInfo): Information about the topology
            out (str): Path to the output file
            d (int): Number of digits to round to
            a (bool): Discard the a vectors
        
        The output will be written to out.
    """
    try:
        remove(out)
    except:
        pass

    ctx = ComputeContext(traj_info, top_info, d, a)

    with open(out, 'w+') as f:
        def callback(i, r):
            nonlocal f
            f.write(r)

        oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    log(f"Wrote aligned trajectory to {out}")

    return

def cli_parser(prog="minify.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Compress given configuration.")
    parser.add_argument('trajectory', type=str, help='the trajectory file you wish to analyze')
    parser.add_argument('outfile',    type=str, help='minified file')
    parser.add_argument('-p', '--parallel', metavar='num_cpus', type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-a', '--no_a', action = 'store_true', help='Discard a vectors.')
    parser.add_argument('-d', '--decimals', type=int, help='Round positions and orientations to the specified number of digits.')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    logger_settings.set_quiet(args.quiet)
    traj_file = args.trajectory
    out = args.outfile

    top_info, traj_info = describe(None, traj_file)

    # -p sets the number of parallel processes
    if args.parallel:
        ncpus = args.parallel
    else:
        ncpus = 1

    # -d sets the decimals of the output
    if args.decimals:
        d = args.decimals
    else:
        d = None

    # -a sets the a vectors to 0
    if args.no_a:
        a = True
    else:
        a = False

    minify(traj_info, top_info, out, d, a, ncpus)
    
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()