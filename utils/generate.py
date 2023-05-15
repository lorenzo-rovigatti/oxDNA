#!/usr/bin/env python

import sys, os
try:
    import numpy as np
except:
    import mynumpy as np
import base
import generators as gen

def double_strands():
    g = gen.StrandGenerator()
    s = base.System([30, 30, 30])
    s.add_strands(g.generate(10))
    s.add_strands(g.generate(10))
    s.add_strands(g.generate(25, dir=[0, 1, 0], start_pos=[0, 0, 7], double=False))
    s.add_strands(g.generate(25, dir=[1, 1, 1], start_pos=[-10, -10, -10]))
    s.print_crepy_output ("prova.mgl")
    s.print_lorenzo_output ("generated.dat", "generated.top")

def tetramers():
    g = gen.TetramerGenerator()
    s = base.System ([20, 20, 20])
    s.add_strands (g.generate())
    s.print_crepy_output ("tetramer.mgl")
    s.print_lorenzo_output ("tetramer.dat", "tetramer.top")

def read_strands(filename='caca.sqs', box_side=50):
    """
    The main() function for this script
    Reads a text file with the following format:
    - Each line contains the sequence for a single strand (A,C,T,G)
	      the nucleotides can be specified by (case insensitive) letter (A, C, G, T), 
			random (X), strong (S) and weak (W).
    - Options:
        - DOUBLE
        - CIRCULAR
        - DOUBLE CIRCULAR

    Ex: Two ssDNA (single stranded DNA)
    ATATATA
    GCGCGCG

    Ex: Two strands, one double stranded, the other single stranded.
    DOUBLE AGGGCT
    CCTGTA

	 Ex: One dsDNA that frays only on one side:
	 DOUBLE SSXXXXXXWW

	 Ex: One dsDNA that frays very little but has TATA in the middle:
	 DOUBLE SSXXXXTATAXXXXSS

    Ex: Two strands, one double stranded circular, the other single stranded circular.
    DOUBLE CIRCULAR AGAGGTAGGAGGATTTGCTTGAGCTTCGAGAGCTTCGAGATTCGATCAGGGCT
    CIRCULAR CCTGTAAGGAGATCGGAGAGATTCGAGAGGATCTGAGAGCTTAGCT
    """
    # we read the sequences from a file; each line is a strand; 
    # prepending DOUBLE tells us to generate a double strand
    # prepending CIRCULAR tells us to generate a (topologically) closed strand 

    # Check filename
    try:
        infile = open (filename, "r")
    except IOError:
        base.Logger.die("%s not found" % filename)

    # Check box_side
    if box_side > 1:
        pass
    else:
        base.Logger.log("Invalid box size (%f). Using default box size." % box_side, level=base.Logger.WARNING)
        box_side = 50

    double = gen.StrandGenerator()
    lines = infile.readlines()
    nlines = len(lines)
    print >> sys.stdout, "Found %i lines" % (nlines)

    box = np.array ([float(box_side), float(box_side), float(box_side)], np.float64)
    s = base.System(box)
    i = 1
    for line in lines:
        line = line.upper().strip()
        # skip empty lines
        if len(line) == 0: continue

        # Set flags
        bool_double = False
        bool_circular = False
        bool_rw = False
        type_strand = 'single strand'
        type_circular = 'linear'
        if line.find('DOUBLE') >= 0:
            bool_double = True
            type_strand = 'duplex'
        if line.find('CIRCULAR') >= 0:
            bool_circular = True
            type_circular = 'circular'
        if line.find('RW') >= 0:
            bool_rw = True
         
        # Parse input and output structures in lorenzo format
        raw_seq = line.split()[-1]

        import random
        raw_seq2 = "" 
        for x in raw_seq:
            if x in ['x', 'X']:
                raw_seq2 += random.choice(['A','C','G','T'])
            elif x in ['s', 'S']:
                raw_seq2 += random.choice(['C','G'])
            elif x in ['w', 'W']:
                raw_seq2 += random.choice(['A','T'])
            else:
                raw_seq2 += x

        seq = [base.base_to_number[x] for x in raw_seq2]
        length = len(raw_seq)
        print >> sys.stdout, "Adding %s %s of %i bases" % (type_circular, type_strand, length)
        cdm = np.random.random_sample(3) * s._box
        axis = np.random.random_sample(3)
        axis /= np.sqrt(np.dot(axis, axis))
        success = False
        while not success:
            if not bool_rw:
                success = s.add_strands(double.generate(len(seq), sequence=seq, dir=axis, \
                   start_pos=cdm, double=bool_double, circular=bool_circular), check_overlap=True)
            else:
                success = s.add_strands(double.generate_rw(sequence=seq, start_pos=cdm))
            
            cdm = np.random.random_sample(3) * s._box
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))
            print >> sys.stdout, "  try again with %i" % (i)
        print >> sys.stdout, "  done %i" % (i)
        i += 1

    s.print_lorenzo_output("generated.dat", "generated.top")


def main():
    """Parse command line options"""

    try:
        box_side = float(sys.argv[1])
        filename = sys.argv[2]
    except:
        base.Logger.die("Usage: %s <%s> <%s>" % (sys.argv[0], "box size", "file with sequences"))
        sys.exit(1)

    # Check filename before input
    try:
        inp = open (filename, 'r')
        inp.close()
    except:
        base.Logger.die("Could not open file '%s' for reading. Aborting" % filename)
        sys.exit(2)
    read_strands(filename, box_side)

if __name__ == "__main__":
        main()
