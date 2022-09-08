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
    base.Logger.log("Usage is %s configuration topology [output] [fixcm]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(sys.argv) > 3:
    output = sys.argv[3]
else: output = sys.argv[1] + ".xyz"

fixcm = False
if len(sys.argv) > 4:
    if sys.argv[4] == "fixcm":
        fixcm = True
    
l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()
append = False
niter = 1
while s:
    base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)

    # set center of mass to strand 0
    if fixcm:
        base.Logger.log("setting centre of mass to strand 0", base.Logger.INFO)
        centre_of_mass = np.array([0.,0.,0.])
        n_nucleotides = 0
        for nucleotide in s._strands[0]._nucleotides:
            centre_of_mass += nucleotide.cm_pos
            n_nucleotides += 1
        centre_of_mass /= float(n_nucleotides)
        for strand in s._strands:
            strand.translate (-centre_of_mass)

    s.print_vmd_xyz_output(output, append=append, same_colors=True, visibility = "caca.vis")
    s = l.get_system()
    append = True
    niter += 1

            
if niter < 2:
    base.Logger.log ("Dind't do anything...")
else:
    base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)

