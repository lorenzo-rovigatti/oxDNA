#!/usr/bin/env python

import base
import readers
import os.path
import sys

if len(sys.argv) < 3:
    base.Logger.log("Usage is %s configuration topology [output]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(sys.argv) > 3:
    output = sys.argv[3]
else: output = sys.argv[1] + ".mgl"
    
l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()
s.print_crepy_output(output, same_colors=True)

base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)

