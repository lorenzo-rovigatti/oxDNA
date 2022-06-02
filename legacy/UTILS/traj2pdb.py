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
else: output = sys.argv[1] + ".pdb"
    
l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()
append = False
while s:
    #s.bring_in_box_nucleotides()	
    s.print_pdb_output(output, append=append)
    s = l.get_system()
    append = True

base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)
