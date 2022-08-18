#!/usr/bin/env python3
import os
from sys import stderr
from re import finditer

import argparse
from typing import List
import numpy as np

from oxDNA_analysis_tools.external_force_utils.force_reader import write_force_file
from oxDNA_analysis_tools.external_force_utils.forces import mutual_trap

# Copied from ExpertRNA
def parse_dot_bracket(input:str) -> List[int]:
    """
    Converts a dot-bracket string to a list of paired nucleotides

    Parameters:
        input (str): A dot-bracket string

    Returns:
        output (list): A list where each index corresponds to a nucleotide.  Value is -1 is unpaired or another index if paired.
    """
    output = np.full(len(input), -1)
    #I'm not sure this is the most efficent way to do this, but I'm lazy.
    more = True
    while more:
        more = False

        #finds matched parenthesis
        for x in finditer(r"\([^()]*\)", input):
            more = True
            output[x.start()] = x.end()-1
            output[x.end()-1] = x.start()

            #its recursive...
            input=input[0:x.start()] + "." + input[x.start()+1:x.end()-1] + "." + input[x.end():]

    return output

def cli_parser(prog="db_to_force.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Create an external forces file enforcing the current base-pairing arrangement")
    parser.add_argument('db_file', type=str, nargs=1, help="A text file containing dot-bracket notation of the base-pairing arrangement")
    parser.add_argument('-o', '--output', type=str, nargs=1, help='Name of the file to write the force list to')
    parser.add_argument('-s', '--strength', type=float, nargs=1, help='Strength of the forces')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    # Get input
    with open(args.db_file[0]) as f:
        db_str = f.read()

    # Check for strength input, otherwise default to 0.09 which won't explode most simulations.
    if args.strength:
        strength = args.strength[0]
        print("INFO: Using strength {}".format(strength), file=stderr)
    else:
        strength = 0.09
        print("INFO: No strength provided, defaulting to {}".format(strength), file=stderr)

    # convert the db string to an index list
    db_idx = parse_dot_bracket(db_str)

    force_list = []

    #p is particle id, q is paired particle id
    for p, q in enumerate(db_idx):
        if q != -1:
            force_list.append(mutual_trap(p, q, strength, 1.2, 1))

    # write the force file
    if args.output:
        outfile = args.output[0]
        print("INFO: Writing forces to {}".format(outfile), file=stderr)
    else:
        outfile = "external_forces.txt"
        print("INFO: No output filename found.  Defaulting to {}".format(outfile), file=stderr)

    write_force_file(force_list, outfile)

if __name__ == '__main__':
    main()