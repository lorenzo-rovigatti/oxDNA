#!/usr/bin/env python3
import os
from sys import stderr
import argparse
from typing import List, Dict
import numpy as np
from oxDNA_analysis_tools.external_force_utils.force_reader import write_force_file
from oxDNA_analysis_tools.external_force_utils.forces import mutual_trap

def parse_dot_bracket(input:str) -> np.ndarray:
    """
    Converts a dot-bracket string to a list of paired nucleotides.

    Accepts (), [], {}, and . characters, otherwise throws an error

    Parameters:
        input (str): A dot-bracket string

    Returns:
        (list): A list where each index corresponds to a nucleotide.  Value is -1 is unpaired or another index if paired.
    """
    output = np.full(len(input), -1)
    paren_queue = []
    square_queue = []
    curly_queue = []

    for i, c in enumerate(input.strip()):
        if c == '.':
            continue
        elif c == '(':
            paren_queue.append(i)
        elif c == '[':
            square_queue.append(i)
        elif c == '{':
            curly_queue.append(i)
        elif c == ')':
            pair = paren_queue.pop()
            output[i] = pair
            output[pair] = i
        elif c == ']':
            pair = square_queue.pop()
            output[i] = pair
            output[pair] = i
        elif c == '}':
            pair = curly_queue.pop()
            output[i] = pair
            output[pair] = i
        else:
            raise RuntimeError("Encountered invalid character '{}' in dot bracket".format(c))

    return output

def db_to_forcelist(db_str:str, stiff:float, reverse:bool, r0:float=1.2, PBC:bool=True, rate:float=0, stiff_rate:float=0) -> List[Dict]:
    """
        Convert a dot-bracket string to oxDNA mutual traps

        Parameters:
            db_str (str): The dot-bracket string to convert.  Currently ignores pseudoknots.
            stiff (float): stiff of the mutual trap to create
            reverse (bool): Reverse the dot-bracket string before creating the forces?

        Returns:
            List[Dict]: A list of force dictionaries
    """
    if reverse:
        db_str = db_str[::-1]
        db_str = db_str.replace('(', '&')
        db_str = db_str.replace(')', '(')
        db_str = db_str.replace('&', ')')
        db_str = db_str.replace('[', '&')
        db_str = db_str.replace(']', '[')
        db_str = db_str.replace('&', ']')
        db_str = db_str.replace('{', '&')
        db_str = db_str.replace('}', '{')
        db_str = db_str.replace('&', '}')

    # convert the db string to an index list
    db_idx = parse_dot_bracket(db_str)

    force_list = []

    #p is particle id, q is paired particle id
    for p, q in enumerate(db_idx):
        if q != -1:
            force_list.append(mutual_trap(p, q, stiff, r0, True, rate=rate, stiff_rate=stiff_rate))

    return force_list

def cli_parser(prog="db_to_force.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Create an external forces file enforcing the current base-pairing arrangement")
    parser.add_argument('db_file', type=str, nargs=1, help="A text file containing dot-bracket notation of the base-pairing arrangement")
    parser.add_argument('-o', '--output', type=str, nargs=1, help='Name of the file to write the force list to')
    parser.add_argument('-s', '--stiff', type=float, nargs=1, help='Stiffness of the mutual trap')
    parser.add_argument('-r', '--reverse', action='store_true', dest='reverse', default=False, help='Reverse the dot-bracket before writing forces')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    with open(args.db_file[0]) as f:
        db_str = f.read()

    # Default trap stiffness is 0.09 which won't explode most simulations.
    if args.stiff:
        stiff = args.stiff[0]
        print("INFO: Using stiffness {}".format(stiff), file=stderr)
    else:
        stiff = 0.09
        print("INFO: No stiffness provided, defaulting to {}".format(stiff), file=stderr)

    reverse = args.reverse

    force_list = db_to_forcelist(db_str, stiff, reverse)

    if args.output:
        outfile = args.output[0]
        print("INFO: Writing forces to {}".format(outfile), file=stderr)
    else:
        outfile = "external_forces.txt"
        print("INFO: No output filename found.  Defaulting to {}".format(outfile), file=stderr)

    write_force_file(force_list, outfile)

if __name__ == '__main__':
    main()