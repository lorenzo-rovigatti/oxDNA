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
else: output = sys.argv[1] + ".xyz"

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()

visibility = [False, False]

append = False
while s:
    s.translate (-s._nucleotides[0].cm_pos)
    r = s.get_reduced ("caca.vis")
    r.print_tcl_output("isolated.dat.tcl")
    r.print_lorenzo_output("isolated.dat", "isolated.top")
    s = l.get_system()
    append = True


