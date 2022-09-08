#!/usr/bin/env python

import base
import readers
import os.path
import sys
import subprocess

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

if len(sys.argv) < 4:
    base.Logger.log("Usage is %s input configuration topology [output]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

inputfile = sys.argv[1]

if len(sys.argv) > 4:
    output = sys.argv[4]
else: output = sys.argv[2] + ".dbk"
    
l = readers.LorenzoReader(sys.argv[2], sys.argv[3])
s = l.get_system()
s.map_nucleotides_to_strands()
launchargs = [PROCESSDIR + 'output_bonds', inputfile, sys.argv[2], "0"]
myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
s.read_H_bonds(myinput.stdout.readlines())

for line in myinput.stderr.readlines():
      if "CRITICAL" in line:
        base.Logger.log("Error running output_bonds'", base.Logger.INFO)

s.print_dot_bracket_output(output)

base.Logger.log("Output printed on '%s'" % output, base.Logger.INFO)

