import argparse
import os
import numpy as np
from sys import stderr
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, inbox, write_conf
from oxDNA_analysis_tools.align import svd_align


def main():
    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="superimposes one or more structures sharing a topology to a reference structure")
    parser.add_argument('reference', type=str, nargs=1, help="The reference configuration to superimpose to")
    parser.add_argument('victims', type=str, nargs='+', help="The configuraitons to superimpose on the reference")
    parser.add_argument('-i', metavar='index_file', dest='index_file', nargs=1, help='Align to only a subset of particles from a space-separated list in the provided file')
    args = parser.parse_args()

    #run system checks
    from oxDNA_analysis_tools.config import check_dependencies
    check_dependencies(["python", "numpy"])

    #Get the reference configuration
    ref_file = args.reference[0]
    top_info, ref_info = describe(None, ref_file)
    ref_conf = get_confs(ref_info.idxs, ref_info.path, 0, 1, top_info.nbases)[0]

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
                print("ERROR: The index file must be a space-seperated list of particles.  These can be generated using oxView by clicking the \"Download Selected Base List\" button")
    else: 
        indexes = list(range(top_info.nbases))

    # alignment requires the ref to be centered at 0.  Inboxing did not take the indexing into account.
    reference_coords = ref_conf.positions[indexes]
    ref_cms = np.mean(reference_coords, axis=0) # cms prior to centering
    reference_coords = reference_coords - ref_cms

    for i, f in enumerate(args.victims):
        top_info, traj_info = describe(None, f)
        conf = get_confs(traj_info.idxs, traj_info.path, 0, 1, top_info.nbases)[0]
        conf = inbox(conf)

        np_coords = np.asarray([conf.positions, conf.a1s, conf.a3s])

        conf.positions, conf.a1s, conf.a3s = svd_align(reference_coords, np_coords, indexes)

        write_conf("aligned{}.dat".format(i), conf)
        print("INFO: Wrote file aligned{}.dat".format(i), file=stderr)

if __name__ == '__main__':
    main()