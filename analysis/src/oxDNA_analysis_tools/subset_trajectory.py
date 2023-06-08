import os
import argparse
from sys import stderr
from collections import namedtuple
from copy import deepcopy
from typing import List
from oxDNA_analysis_tools.UTILS.oat_multiprocesser import oat_multiprocesser
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, strand_describe, conf_to_str, get_top_string
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, System, TopInfo, TrajInfo



ComputeContext = namedtuple("ComputeContext",["traj_info",
                                              "top_info",
                                              "indexes"])

def compute(ctx:ComputeContext, chunk_size:int, chunk_id:int):
    confs = get_confs(ctx.top_info, ctx.traj_info, chunk_id*chunk_size, chunk_size)
    outstr = [[] for _ in range(len(ctx.indexes))]
    for conf in confs:
        sub_confs = [Configuration(conf.time, conf.box, conf.energy, conf.positions[i], conf.a1s[i], conf.a3s[i]) for i in ctx.indexes]
        for i, sub_conf in enumerate(sub_confs):
            outstr[i].append(conf_to_str(sub_conf,  include_vel=ctx.traj_info.incl_v))

    return [''.join(out) for out in outstr]

def write_topologies(system:System, indexes:List[List[int]], outfiles:List[str], old_format=False):
    top_names = [o+ ".top" for o in outfiles]
    for idx, top_name in zip(indexes, top_names):
        idx = set(idx)
        new_sys = deepcopy(system)
        for s in new_sys:
            if s[0].n3 != None and s[0].id != s[-1].n5:
                print(f"WARNING: Strand {s.id} is circular. Subsetting the trajectory will cut circular strands", file=stderr)
            s.monomers = [n for n in s if n.id in idx]
            
        new_sys.strands = [s for s in new_sys if len(s.monomers) > 0]

        with open(top_name, 'w+') as f:
            f.write(get_top_string(new_sys, old_format))

    return top_names

def subset(traj_info:TrajInfo, top_info:TopInfo, system:System, indexes:List[List[int]], outfiles:List[str], ncpus=1):
    """
        Splits a trajectory into multiple trajectories, each containing a subset of the particles in the original configuration.

        Parameters:
            traj_info (TrajInfo): Information about the trajectory
            top_info (TopInfo): Information about the topology
            system (System): The system object describing the topology of the original structure.
            indexes (List[List[int]]): A list of lists of indexes.  Each list corresponds to one set of particles to include in one of the output trajectories.
            outfiles (List[str]): A list of output file names.  The number of elements in this list must match the number of elements in the indexes list.
            ncpus (int): (optional) The number of CPUs to use

        The output trajectories will be written to the files specified in outfiles.
    """
    # Create a ComputeContext which defines the problem to pass to the worker processes 
    ctx = ComputeContext(traj_info, top_info, indexes)

    dat_names = [o+ ".dat" for o in outfiles]
    files = [open(f, 'w+') for f in dat_names]
    def callback(i, r):
        for f, subset in zip(files, r):
            f.write(subset)

    oat_multiprocesser(traj_info.nconfs, ncpus, compute, callback, ctx)

    for f in files:
        f.close()

    # Write topology files
    top_names = write_topologies(system, indexes, outfiles, system.strands[0].is_old())

    print("INFO: Wrote trajectories: {}".format(dat_names), file=stderr)
    print("INFO: Wrote topologies: {}".format(top_names), file=stderr)

def cli_parser(prog="subset_trajectory.py"):
    #command line arguments
    parser = argparse.ArgumentParser(prog = prog, description="Extracts parts of a structure into separate trajectories")
    parser.add_argument('trajectory', type=str, help="The trajectory file to subset")
    parser.add_argument('topology', type=str, help="The topology file corresponding to the trajectory")
    parser.add_argument('-i', '--index', metavar='index', action='append', nargs=2, help='A space separated index file and the associated output file name.  This can be called multiple times')
    parser.add_argument('-p', metavar='num_cpus', type=int, dest='parallel', help="(optional) How many cores to use")
    parser.add_argument('-f', action='store_true', dest='old_format', help="Use the old 3'-5' topology format?")
    return(parser)

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    top_file  = args.topology
    traj_file = args.trajectory
    index_files = [i[0] for i in args.index]
    output_files = [i[1] for i in args.index]
    top_info, traj_info = describe(top_file, traj_file)
    system, _ = strand_describe(top_file)
    indexes = []
    outfiles = []
    for i, o in zip(index_files, output_files):
        with open(i) as f:
            data = f.readline().split()
            try:
                data = sorted([int(i) for i in data])
            except:
                raise RuntimeError("The index file {} must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button".format(i))
        indexes.append(data)
        outfiles.append(o)

    if args.parallel:
        ncpus = args.parallel
    else:
        ncpus = 1

    subset(traj_info, top_info, system, indexes, outfiles, ncpus)

if __name__ == '__main__':
    main()