#!/usr/bin/env python3
import os
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
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
        np.ndarray: A list where each index corresponds to a nucleotide.  Value is -1 is unpaired or another index if paired.
    """
    open_symbols = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    close_symbols = ")]}>abcdefghijklmnopqrstuvwxyz"

    output = np.full(len(input), -1)
    queue = []

    for i, c in enumerate(input.strip()):
        if c == '.':
            continue

        depth = open_symbols.find(c)
        is_open = True
        if depth == -1:
            depth = close_symbols.find(c)
            is_open = False
            if depth == -1:
                raise RuntimeError("Encountered invalid character '{}' in dot bracket".format(c))
                                
        if depth >= len(queue):
            queue.append([])
        
        if is_open:
            queue[depth].append(i)
        else:
            pair = queue[depth].pop()
            output[i] = pair
            output[pair] = i

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
    open_symbols = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    close_symbols = ")]}>abcdefghijklmnopqrstuvwxyz"
    
    if reverse:
        db_str = db_str[::-1]
        for o, c in zip(open_symbols, close_symbols):
            db_str = db_str.replace(o, '&')
            db_str = db_str.replace(c, o)
            db_str = db_str.replace('&', c)

    # convert the db string to an index list
    db_idx = parse_dot_bracket(db_str)

    force_list = []

    #p is particle id, q is paired particle id
    for p, q in enumerate(db_idx):
        if q != -1:
            force_list.append(mutual_trap(p, q, stiff, r0, PBC=PBC, rate=rate, stiff_rate=stiff_rate))

    return force_list

def cli_parser(prog="db2forces.py"):
    parser = argparse.ArgumentParser(prog = prog, description="Create an external forces file enforcing the current base-pairing arrangement")
    parser.add_argument('db_file', type=str, nargs=1, help="A text file containing dot-bracket notation of the base-pairing arrangement")
    parser.add_argument('-o', '--output', type=str, nargs=1, help='Name of the file to write the force list to')
    parser.add_argument('-s', '--stiff', type=float, nargs=1, help='Stiffness of the mutual trap')
    parser.add_argument('-r', '--reverse', action='store_true', dest='reverse', default=False, help='Reverse the dot-bracket before writing forces')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    with open(args.db_file[0]) as f:
        db_str = f.read()

    # The file isn't just a db
    db_lines = db_str.split('\n')
    if len(db_lines) > 1:
        db_str = ''
        for l in db_lines:
            if l[0] == '(' or l[0] == '.':
                db_str = l
                break
    if db_str == '':
        raise RuntimeError("Dot-bracket strings must start with '.' or '('!  No dot-bracket string found in file")

    # Default trap stiffness is 0.09 which won't explode most simulations.
    if args.stiff:
        stiff = args.stiff[0]
        log("Using stiffness {}".format(stiff))
    else:
        stiff = 0.09
        log("No stiffness provided, defaulting to {}".format(stiff))

    reverse = args.reverse

    force_list = db_to_forcelist(db_str, stiff, reverse)

    if args.output:
        outfile = args.output[0]
        log("Writing forces to {}".format(outfile))
    else:
        outfile = "external_forces.txt"
        log("No output filename found.  Defaulting to {}".format(outfile))

    write_force_file(force_list, outfile)

if __name__ == '__main__':
    main()