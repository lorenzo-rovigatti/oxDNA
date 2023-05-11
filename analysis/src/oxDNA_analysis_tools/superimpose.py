import argparse
import os
from typing import List
import numpy as np
from sys import stderr
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, inbox, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration
from oxDNA_analysis_tools.align import svd_align


def superimpose(ref:Configuration, victims:List[str], indexes:List[int]=[]):
    """
    Superimposes one or more structures sharing a topology to a reference structure

    Parameters:
        ref (Configuration) : the reference configuration to superimpose to
        victims (List[Configuration]) : the configurations to superimpose on the reference
        indexes (List[int]) : the indexes of the particles to superimpose on the reference (default: all)

    Returns:
        Aligned configurations (List[Configuration])
    """

    if indexes == []:
        indexes = list(range(len(ref.positions)))

    ref = inbox(ref)
    # alignment requires the ref to be centered at 0.  Inboxing did not take the indexing into account.
    reference_coords = ref.positions[indexes]
    aligned = []

    for i, f in enumerate(victims):
        top_info, traj_info = describe(None, f)
        conf = get_confs(top_info, traj_info, 0, 1)[0]
        conf = inbox(conf)

        np_coords = np.asarray([conf.positions, conf.a1s, conf.a3s])

        conf.positions, conf.a1s, conf.a3s = svd_align(reference_coords, np_coords, indexes)
        aligned.append(conf)

    return aligned

def cli_parser(prog="superimpose.py"):
    parser = argparse.ArgumentParser(prog = prog, description="superimposes one or more structures sharing a topology to a reference structure")
    parser.add_argument('reference', type=str, nargs=1, help="The reference configuration to superimpose to")
    parser.add_argument('victims', type=str, nargs='+', help="The configurations to superimpose on the reference")
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))    
    args = parser.parse_args()

    #run system checks
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    #Get the reference configuration
    ref_file = args.reference[0]
    top_info, ref_info = describe(None, ref_file)
    ref_conf = get_confs(top_info, ref_info, 0, 1)[0]

    ref_conf = inbox(ref_conf)

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
        indexes = list(range(top_info.nbases))

    aligned = superimpose(ref_conf, args.victims, indexes)

    for i, conf in enumerate(aligned):
        write_conf("aligned{}.dat".format(i), conf, include_vel=ref_info.incl_v)
        print("INFO: Wrote file aligned{}.dat".format(i), file=stderr)

if __name__ == '__main__':
    main()