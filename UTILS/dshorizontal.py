#!/usr/bin/env python

import sys
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
    s.print_lorenzo_output ("prova.conf", "prova.top")
    
def tetramers():
    g = gen.TetramerGenerator()
    s = base.System ([20, 20, 20])
    s.add_strands (g.generate())
    s.print_crepy_output ("tetramer.mgl")
    s.print_lorenzo_output ("tetramer.conf", "tetramer.top")

def read_strands():
    # prendiamo le sequenze da un file; ogni riga uno strand; se e' un
    # double strand, ci deve essere DOUBLE davanti
    double = gen.StrandGenerator()
    try:
        infile = open ("caca.sqs", "r")
    except IOError:
        base.Logger.die("caca.sqs not found")

    if len(sys.argv) > 1:
        side = float(sys.argv[1])
    else: side = 50

    lines = infile.readlines()
    nlines = len(lines)
    print >> sys.stderr, "Found %i lines" % (nlines)
    s = base.System(np.array ([float(side), float(side), float(side)], np.float64))
    i = 1 
    for line in lines:
        line = line.upper().strip()
        # skip empty lines
        if len(line) == 0: continue
        if line[:6] == 'DOUBLE':
            line = line[6:]
            line = line.split()[-1]
            seq = [base.base_to_number[x] for x in line]
            length = len(line)
            print >> sys.stderr, "Adding duplex of %i bases" % (length)
            cdm = np.random.random_sample(3) * s._box
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))
            axis = np.array ([1.,0.,0])
            cdm = np.array ([0.,0.,0.])
            while not s.add_strands(double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm)):
                cdm = np.random.random_sample(3) * s._box
                axis = np.random.random_sample(3)
                axis /= np.sqrt(np.dot(axis, axis))
                print >> sys.stderr, "  riprovo %i" % (i)
            print >> sys.stderr, "  done %i" % (i)
        else:
            seq = [base.base_to_number[x] for x in line]
            cdm = np.random.random_sample(3) * s._box
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))
            print >> sys.stderr, "Adding single strand of %i bases" % (len(line))
            while not s.add_strand(double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=False)):
                cdm = np.random.random_sample(3) * s._box
                axis = np.random.random_sample(3)
                axis /= np.sqrt(np.dot(axis, axis))
                print >> sys.stderr, "  riprovo %i" % (i)
            print >> sys.stderr, "  done %i" % (i)
        #print >> sys.stderr, "Added strand..."
        i += 1

    s.print_crepy_output("prova.mgl")
    s.print_lorenzo_output("prova.conf", "prova.top")

read_strands()

