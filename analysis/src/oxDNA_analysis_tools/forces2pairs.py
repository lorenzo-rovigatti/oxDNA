#!/usr/bin/env python

#Created by: Erik Poppleton
#Date: 6/29/18
#Python2
#Converts the forces file printed out by tiamat2oxdna to a pairs file containing all designed H-bonds

import sys
import argparse
import os
from typing import List, Tuple

def cli_parser(prog="forces2pairs"):
    parser = argparse.ArgumentParser(prog = prog, description="Convert an external force file to a list of particle pairs")
    parser.add_argument('force_file', type=str, nargs=1, help="The force file to generate pairs from")
    parser.add_argument('-o', '--output', type=str, nargs=1, help='name of the file to write the pair list to')
    return parser

def forces2pairs(force_file:str) -> List[Tuple]:
    """
    Returns a list of tuples containig the pairing information for a structure

    Parameters:
        force_file (str): path to the force file

    Returns:
        pairs (List[Tuple]): A list of tuples where each tuple corresponds to a pair found in the force file.
    """
    pairs = []
    a = b = -1
    with open(force_file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith("particle"):
                a = int(float(line.split("=")[1].strip()))
            if "ref_particle" in line:
                b = int(float(line.split("=")[1].strip()))
            if "}" in line:
                if a < b:
                    pairs.append((a, b))
                a = b = -1

    return pairs


def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    infile = args.force_file[0]

    try:
        out = args.output[0]
    except:
        print("INFO: No outfile provided, defaulting to pairs.txt", file=sys.stderr)
        out = 'pairs.txt'


    pairs = forces2pairs(infile)

    with open(out, 'w+') as f:
        for p in pairs:
            f.write('{} {}\n'.format(p[0], p[1]))

    print("INFO: pairing information written to {}".format(out))

if __name__ == '__main__':
    main()
    