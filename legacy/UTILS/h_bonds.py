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
import tempfile

command_for_data =  'analysis_data_output_1 = { \n name = stdout \n print_every = 1 \n col_1 = { \n type=pair_energy \n} \n}'
PROCESSPROGRAM = os.path.join(os.path.dirname(__file__), "../bin/DNAnalysis")
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
    line = line.lstrip()
    if not line.startswith('#'):
        if "topology" in line:
            topologyfile = line.split('=')[1].replace(' ','').replace('\n','')

myreader = readers.LorenzoReader(conffile,topologyfile)
mysystem = myreader.get_system()

if not os.path.isfile(PROCESSPROGRAM):
	print "Cannot execute output_bonds program. Please make sure to compile DNAnalysis in ../bin/ directory"
	sys.exit(1)

counter = 0



tempfile_obj = tempfile.NamedTemporaryFile()
launchcommand = inputfile + ' trajectory_file='+tempfile_obj.name+' '+command_for_data

launchargs = [PROCESSPROGRAM,inputfile ,'trajectory_file='+tempfile_obj.name,command_for_data]
#print command_for_data
#launchargs = [PROCESSPROGRAM,inputfile ,'trajectory_file='+conffile,command_for_data]

while mysystem != False:
	mysystem.map_nucleotides_to_strands()
	mysystem.print_lorenzo_output(tempfile_obj.name,'/dev/null')
	tempfile_obj.flush()
	#launchcommand = 'trajectory_file='+tempfile_obj.name+' '+command_for_data
	#os.system(PROCESSPROGRAM+' '+launchcommand)
	#launchargs = [PROCESSPROGRAM,launchcommand]
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout,stderr = myinput.communicate()
	linewise = stdout.split('\n')
	#print 'FEEDING:', linewise
	mysystem.read_H_bonds(linewise[:-1])

	for line in stderr.split('\n'):
      	  if "CRITICAL" in line:
              	  print line
	print '# configuration number: ',counter
	mysystem.show_H_interactions()
	counter += 1
	mysystem = myreader.get_system()


