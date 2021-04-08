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
#print PROCESSDIR

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

counter = 0

while mysystem != False:
	mysystem.map_nucleotides_to_strands()
	launchargs = [PROCESSDIR + 'output_bonds',inputfile,conffile,str(counter)]
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
 	stdout,stderr = myinput.communicate()
	linewise = stdout.split('\n')
	mysystem.read_H_bonds(linewise[:-1])

	for line in stderr.split('\n'):
      	  if "CRITICAL" in line:
              	  print line
	print '# configuration number: ',counter
	mysystem.show_H_interactions()
	counter += 1
	mysystem = myreader.get_system()


