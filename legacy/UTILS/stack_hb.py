#!/usr/bin/env python

#A utility that prints out the number of hydrogen bonds between different strands in the system 

import base
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import readers 
import subprocess

def get_total_energies (lines):
    estack, ehb = 0., 0.
    for line in lines:
        words = line.split()
        if len(words) == 10 and words[0][0] != '#':
            estack += float (words[4])
            ehb += float (words[6])
    return estack, ehb

def nhbst (lines):
    nst, nhb = 0, 0
    for line in lines:
        words = line.split ()
        if len(words) == 10 and words[0][0] != '#':
            estack = float (words[4])
            ehb = float (words[6])
            if estack < -0.1:
                nst += 1
            if ehb < -0.1:
                nhb += 1
    return nhb, nst

#PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")
PROCESSPROGRAM = os.path.join(os.path.dirname(__file__), "output_bonds.py")

if (len(sys.argv) < 3):
  print 'Usage %s input_file trajectory_file ' % sys.argv[0]
  sys.exit()
  
#now get topology file name:
inputfile = sys.argv[1]
conffile = sys.argv[2]
topologyfile = ""
fin = open(inputfile)
for line in fin:
  if "topology" in line:
    topologyfile = line.split('=')[1].replace(' ','').replace('\n','')

myreader = readers.LorenzoReader(conffile,topologyfile)
mysystem = myreader.get_system()

if not os.path.isfile(PROCESSPROGRAM):
	print "Cannot execute output_bonds program. Please make sure to go to process_data/ directory and type make"
	sys.exit(1)

counter = 0
while mysystem != False:
    mysystem.map_nucleotides_to_strands()
    launchargs = [PROCESSPROGRAM,inputfile,conffile,str(counter)]
    myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    #mysystem.read_H_bonds(myinput.stdout.readlines())
    a, b = nhbst (myinput.stdout.readlines())
    for line in myinput.stderr.readlines():
        if "CRITICAL" in line:
            print line
    print mysystem._time, a, b
    counter += 1
    mysystem = myreader.get_system()

