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

PROCESSDIR = os.path.join(os.path.dirname(__file__), "process_data/")

if (len(sys.argv) < 4):
  print 'Usage %s <input> <trajectory> <group>' % sys.argv[0]
  sys.exit()

group = int(sys.argv[3])

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

counter = 0

bonds = []
while mysystem != False:
    print >> sys.stderr, "# Mapping nucleotides..."
    mysystem.map_nucleotides_to_strands()
    print >> sys.stderr, "#         Done."
    launchargs = [PROCESSDIR + 'output_bonds',inputfile,conffile,str(counter)]
    print >> sys.stderr, "# Running output_bonds..."
    myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    print >> sys.stderr, "#         Done."
    print >> sys.stderr, "# Computing CSD..."

    mysystem.read_H_bonds(myinput.stdout.readlines())
    
    nmax = mysystem.N_strands
    bonds = [s.H_interactions.keys() for s in mysystem._strands]
    # C.S.D.
    clust = [i/group for i in range (nmax)]
    inclust = [0 for i in range (nmax)]
    
    check = True
    while check:
        check = False
        for i in xrange (nmax):
            for j in bonds[i]:
                if clust[i] != clust[j]:
                    check = True
                    #print clust[i], clust[j], "-->",
                    csn = min (clust[i], clust[j])
                    clust[i], clust[j] = csn, csn
                    #print clust[i], clust[j]
    
    for i in xrange (nmax):
        inclust[clust[i]] += 1
     
    inclust.sort()
    inclust.reverse()
    sum = 0
    for i in inclust:
        if i > group:
            print i/group,
    print

    counter += 1
    mysystem = myreader.get_system()


