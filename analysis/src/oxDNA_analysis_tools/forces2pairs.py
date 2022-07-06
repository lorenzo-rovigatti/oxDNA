#!/usr/bin/env python

#Created by: Erik Poppleton
#Date: 6/29/18
#Python2
#Converts the forces file printed out by tiamat2oxdna to a pairs file containing all designed H-bonds

import sys
import argparse
import os

def main():

    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Convert an external force file to a list of particle pairs")
    parser.add_argument('force_file', type=str, nargs=1, help="The force file to generate pairs from")
    parser.add_argument('-o', '--output', type=str, nargs=1, help='name of the file to write the pair list to')

    args = parser.parse_args()

    infile = args.force_file[0]

    #if there's no outfile argument, just print to stdout.
    try:
        out = args.output[0]
        outfile = open(out, 'w+')

    except:
        outfile = sys.stdout

    #Process the forces file
    with open(infile) as f:
        for i, line in enumerate(f.readlines()):
            line = line.strip()
            if line.startswith("particle"):
                a = line.split("=")[1].strip()
            if "ref_particle" in line:
                b = line.split("=")[1].strip()
            if "}" in line:
                if int(a) < int(b):
                    print(a, b, sep = ' ', file=outfile)
                a = -1
                b = -1

if __name__ == '__main__':
    main()
    