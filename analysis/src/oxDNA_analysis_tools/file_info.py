import os
from typing import List, Dict
import argparse
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs

def file_info(trajectories:List) -> Dict:
    """
        Extract metadata about trajectory files

        Parameters:
            trajectories (List[str]) : Filepaths to the trajectories to analyze

        Returns:
            info (Dict) : Dictionary with name, number of particles, number of confs, filezie, starting step and ending step.
    """
    
    info = {
        'name' : [],
        'particles' : [],
        'n_confs' : [],
        'traj_size' : [],
        't_start' : [],
        't_end' : [],
    }

    # Get info from each trajectory
    for t in trajectories:
        info['name'] = t
        ti, di = describe(None, t)
        info['particles'].append(ti.nbases)
        info['n_confs'].append(di.nconfs)
        info['traj_size'].append(os.stat(t).st_size / 1000000)

        first_conf = get_confs(ti, di, 0, 1)[0]
        info['t_start'].append(first_conf.time)

        last_conf = get_confs(ti, di, di.nconfs-1, 1)[0]
        info['t_end'].append(last_conf.time)

    return (info)


def print_info(info, labels):
    # prints each column nicely justified/padded
    def print_value(v, pad):
        print(str(v).ljust(pad+2), end='')

    # Munge data for printing
    info['name'] = labels
    info['traj_size'] = ['{:.2f} MB'.format(s) for s in info['traj_size']]
    info['t_start'] = ['{:.3g}'.format(t) for t in info['t_start']]
    info['t_end'] = ['{:.3g}'.format(t) for t in info['t_end']]

    # Get maximum value in each column
    w = [max([len(k), max([len(str(v)) for v in info[k]])]) for k in info.keys()]

    # Print column headers
    [print(k.ljust(w[i]+2), end='') for i, k in enumerate(info.keys())]
    print()
    
    # Print rows
    for i, _ in enumerate(labels):
        for j, k in enumerate(info.keys()):
            print_value(info[k][i], w[j])
        print()

def cli_parser(prog="file_info.py"):
    #handle commandline arguments
    parser = argparse.ArgumentParser(prog = prog, description="Prints metadata about trajectories")
    parser.add_argument('trajectories', type=str, nargs='+', help='One or more trajectories to get information on.')
    parser.add_argument('-l', '--labels', type=str, nargs='+', help='Labels for the files if not the filename')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    # Verify that dependencies are installed and a good version
    from oxDNA_analysis_tools.config import check
    check(["python"])

    # Command line arguments
    trajectories = args.trajectories
    if args.labels:
        labels = args.labels
        if len(labels) != len(trajectories):
            raise RuntimeError("Number of trajectories does not match the number of labels.")
    else:
        labels = [os.path.basename(t) for t in trajectories]

    info = file_info(trajectories)

    print_info(info, labels)

if __name__ == '__main__':
    main()
