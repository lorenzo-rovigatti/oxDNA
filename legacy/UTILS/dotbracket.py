#!/usr/bin/env python

#A utility that prints out the number of hydrogen bonds between different strands in the system 
# specify the secondary structure of which strands you want to display. If you set both strands equal, it shows the bonds of the strand with itself
#Still in beta version, tested only on few simple systems. Currently does not support pseudoknots!
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

def conf2dot(system,sid1,sid2):
	if sid1 > sid2:
		return 'Error, sid1 > sid2'
	offset = 0
	if sid1 != sid2:
		offset = 1+len(system._strands[sid1]._nucleotides)

	if sid1 == sid2:
		result = ['.'] * len(system._strands[sid1]._nucleotides)
	elif sid1 != sid2:
		result = ['.'] * len(system._strands[sid1]._nucleotides)  + ['+'] +  ['.'] * len(system._strands[sid2]._nucleotides)
				
	for nucleotide_id in range(len(system._strands[sid1]._nucleotides)):
		nucleotide = system._strands[sid1]._nucleotides[nucleotide_id]
		if(len(nucleotide.interactions) > 1):
			print >> sys.stderr, 'One nucleotide has more than 1 h-bonding interaction, skipping...'
			return ''
		elif len(nucleotide.interactions) == 1:
			second_id =  system._nucleotides[nucleotide.interactions[0]].index 
			sid = system._nucleotide_to_strand[second_id]
			if (sid == sid2 ) and result[nucleotide_id] == '.':
				local_id = offset + second_id - system._strands[sid]._nucleotides[0].index
				result[nucleotide_id] = '(' 
				result[local_id] = ')'

			elif (sid == sid1 ) and result[nucleotide_id] == '.':
				local_id =  second_id - system._strands[sid]._nucleotides[0].index
				result[nucleotide_id] = '(' 
				result[local_id] = ')'
	if sid1 != sid2:
		
		for nucleotide_id in range(len(system._strands[sid2]._nucleotides)):
			nucleotide = system._strands[sid1]._nucleotides[nucleotide_id]
			if(len(nucleotide.interactions) > 1):
				print >> sys.stderr, 'One nucleotide has more than 1 h-bonding interaction, skipping...'
				return ''
			elif len(nucleotide.interactions) == 1:
				second_id =  system._nucleotides[nucleotide.interactions[0]].index 
				sid = system._nucleotide_to_strand[second_id]
				if (sid == sid2) and result[nucleotide_id] == '.':
					local_id = offset + second_id - system._strands[sid]._nucleotides[0].index
					result[nucleotide_id] = '(' 
					result[local_id] = ')'

	return ''.join(result)					
					

if (len(sys.argv) < 5):
  print 'Usage %s input_file trajectory_file strand_id1 strand_id2' % sys.argv[0]
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

#strand ids
sid1 = int(sys.argv[3])
sid2 = int(sys.argv[4])

if(sid1 > sid2):
	pom = sid2
	sid2 = sid1
	sid1 = pom

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
	#mysystem.show_H_interactions()
	notation = conf2dot(mysystem,sid1,sid2)
	print notation
	counter += 1
	mysystem = myreader.get_system()


