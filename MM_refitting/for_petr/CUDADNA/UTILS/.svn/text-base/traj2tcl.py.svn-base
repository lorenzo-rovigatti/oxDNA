#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys

if len(sys.argv) < 3:
    base.Logger.log("Usage is %s configuration topology [output]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(sys.argv) > 3:
    output = sys.argv[3]
else: output = sys.argv[1] + ".tcl"
    
l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()
append = False
niter = 1
while s:
    base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
	
    for strand in s._strands:
		diff = np.rint(strand.cm_pos / s._box ) * s._box
		s.translate (-diff)

    cdm = np.array ([0.,0.,0.])
    for strand in s._strands:
        for n in strand._nucleotides:
            cdm += n.cm_pos
    cdm = cdm / float (s.get_N_Nucleotides())
    
    s.print_tcl_output(output, visibility="caca.vis")
    s = l.get_system()
    append = True
    niter += 1

base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)

