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
        if len(words) == 10:
            estack += float (words[3])
            ehb += float (words[2])
    return estack, ehb

def nhbst (lines):
    nst, nhb = 0, 0
    for line in lines:
        words = line.split ()
        if len(words) == 10:
            estack = float (words[3])
            ehb = float (words[2])
            if estack < -0.1:
                nst += 1
            if ehb < -0.1:
                nhb += 1
    return nhb, nst

def stackarray (lines, nbases):
    ret = [False for i in xrange(nbases - 1)]
    for line in lines:
        words = line.split ()
        if len(words) == 10:
            i = int (words[0])
            j = int (words[1])
            if j - i == 1:
                estack = float (words[4])
                if estack < -0.1:
                    ret[i] = True
    return ret

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

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

if not os.path.isfile(PROCESSDIR + "output_bonds"):
	print "Cannot execute output_bonds program. Please make sure to go to process_data/ directory and type make"
	sys.exit(1)

total = [0 for i in xrange(mysystem._N - 1)]
lengths = [0 for i in xrange(mysystem._N)]


counter = 0
while mysystem != False:
    mysystem.map_nucleotides_to_strands()
    launchargs = [PROCESSDIR + 'output_bonds',inputfile,conffile,str(counter)]
    myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    #mysystem.read_H_bonds(myinput.stdout.readlines())
    a = stackarray (myinput.stdout.readlines(), mysystem._N)
    for line in myinput.stderr.readlines():
        if "CRITICAL" in line:
            print line
    for i in xrange(len(a)):
        if a[i]:
            total[i] += 1

    i = 0
    while i < len(a):
        l = 0
        while a[i] == True:
            l += 1
            i += 1
            if i == len(a):
                break
        lengths[l] += 1
        i += 1

    print >> sys.stderr, mysystem._time, a
    print >> sys.stderr, mysystem._time, lengths
    print >> sys.stderr
    counter += 1
    mysystem = myreader.get_system()

print '# stacking probabilities'
for i in xrange(len(total)):
    print total[i]/float(counter)
print

print '# expected number of regions of n stacks'
for i in xrange (len (lengths)):
    print i, lengths[i]/float(counter)
print

