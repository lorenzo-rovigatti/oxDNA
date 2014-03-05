#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import subprocess

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

if len(sys.argv) < 4:
    base.Logger.log("Usage is %s configuration topology input [output] [cdmto0] [strand]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(sys.argv) > 4:
    output = sys.argv[4]
else: output = sys.argv[1] + ".tcl"

cdm20 = False
if len (sys.argv) > 5:
    if sys.argv[5] == "cdmto0":
        cdm20 = True

centre_strand_true = False
if len (sys.argv) > 6:
    centre_strand_true = True
    centre_strand = int(sys.argv[6])

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()

infile = sys.argv[3]
conffile = sys.argv[1]


s.map_nucleotides_to_strands()
try:
    open(infile)
except:
    base.Logger.log("unable to find file %s, exit" % infile, base.Logger.CRITICAL)
    sys.exit()
launchargs = [PROCESSDIR + 'output_bonds',infile,conffile,str(0)]
myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
s.read_H_bonds(myinput.stdout.readlines())
base.Logger.log("Working on conf 1", base.Logger.INFO)

# bring in box
for strand in s._strands:
            diff = np.rint (strand.cm_pos / s._box ) * s._box
            s.translate (-diff)

# center of mass to 0
if cdm20:
    # centre of mass of centre_strand to zero
    if centre_strand_true:
        base.Logger.log("setting centre of mass of strand %s to zero" %centre_strand, base.Logger.INFO)
        centre_of_mass = np.array([0.,0.,0.])
        n_nucleotides = 0
        for nucleotide in s._strands[centre_strand]._nucleotides:
            centre_of_mass += nucleotide.cm_pos
            n_nucleotides += 1
        centre_of_mass /= float(n_nucleotides)
        for strand in s._strands:
            strand.translate (-centre_of_mass)
    # centre of mass of system to zero
    else:
        base.Logger.log ("setting cdm to 0", base.Logger.INFO)
        cdm = np.array ([0.,0.,0.])
        for strand in s._strands:
            for n in strand._nucleotides:
                cdm += n.cm_pos
        cdm = cdm / float (s.get_N_Nucleotides())
        for strand in s._strands:
            strand.translate (-cdm)

s.print_tcl_cylinder_output(output, visibility="caca.vis", show_labels = False, show_box = False)

base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)

