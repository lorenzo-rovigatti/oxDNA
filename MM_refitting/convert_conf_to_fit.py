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

def get_nuc_in_strand_index(strand,global_index):
	for i in range(len(strand._nucleotides)):
		if strand._nucleotides[i].index == global_index:
			return i
	raise IOError('cannot convert nucid to strandid for '+ str(global_index) + ' in strand ' )  #if we got hwew, it means the nucleotiae is not in the strand 

def get_system_to_interaction_line(system,energy_fraction,weight=1):
	outline = str(weight) + ' ' 
	total_stack = 0
	H_bond_line = ''
	ST_line = ''
	
        for strand in system._strands:
	  for i in range(len(strand._nucleotides)):
		id1 = strand._nucleotides[i].index
                nucl = strand._nucleotides[i]
		for id2 in strand._nucleotides[i].all_interactions[base.INT_STACK].keys():
 			energy = nucl.all_interactions[base.INT_STACK][id2]
			if(id1 < id2 and energy < 0):
				sid = system._nucleotide_to_strand[id2]
				j = get_nuc_in_strand_index(system._strands[sid],id2)
				total_stack += nucl.all_interactions[base.INT_STACK][id2]
				ST_line += ' ST %d %d %d %d %f' % (system._nucleotide_to_strand[id1],i,system._nucleotide_to_strand[id2],j,nucl.all_interactions[base.INT_STACK][id2]) 	
		
        for strand in system._strands:
	  for i in range(len(strand._nucleotides)):
                nucl = strand._nucleotides[i]
		id1 = strand._nucleotides[i].index
		for id2 in nucl.all_interactions[base.INT_HYDR].keys():
 			energy = nucl.all_interactions[base.INT_HYDR][id2]
			if(id1 < id2 and energy < 0):
				sid = system._nucleotide_to_strand[id2]	
				j = get_nuc_in_strand_index(system._strands[sid],id2)
				H_bond_line += ' H %d %d %d %d %f' % (system._nucleotide_to_strand[id1],i,system._nucleotide_to_strand[id2],j,nucl.all_interactions[base.INT_HYDR][id2]) 	
        #print 'System pot en is ', system.E_pot, total_stack
	Etotal = system._N*system.E_pot - (1.0 - energy_fraction) * total_stack	
	outline =str(weight) + ' ' + str(Etotal) + ' '   + ST_line + ' ' + H_bond_line 
	return outline 


#the format of the fit file is:
#DUPLEXMELT <NO_OF_NUCLEOTIDES> <SIMULATION_BETA> <MIN_FIT_BETA> <FIT_BETA_STEP> <FIT_BETA_ARRAY_LENGTH>
#ORIGINAL PARAMETERS
#weight_of_state Energy_of_state(without entropic contribution in stacking) H strand1 nucid strandid2 nucid2 energy  ST ....



PROCESSDIR = os.path.join("../process_data/")

if (len(sys.argv) < 3):
  print 'Usage %s input_file trajectory_file ' % sys.argv[0]
  sys.exit()


  
#now get topology file name:
inputfile = sys.argv[1]
conffile = sys.argv[2]
topologyfile = ""
fin = open(inputfile)
T = 0
for line in fin:
  if "topology" in line:
    topologyfile = line.split('=')[1].replace(' ','').replace('\n','')
  if line.split('=')[0].replace(' ','') == 'T':
	TP = line.split('=')[1].replace(' ','').replace('\n','')
	if TP[len(TP)-1] == 'K':
		T = float(TP[:-1])
	elif TP[len(TP)-1] == 'c' or TP[len(TP)-1] == 'C':
		T = float(TP[:-1]) + 273.15
        else:
		T = float(TP[:-1])*3000.0
if T == 0:
  print >> sys.stderr, 'Did not manage to detect temperature'
  raise IOError('No tempereature detected')

beta = 3000.0/T

print >> sys.stderr, 'Reading T=',T,' beta ',beta 
beta_min = beta - 1.2
beta_step = 0.1
array_length = 22




myreader = readers.LorenzoReader(conffile,topologyfile)
mysystem = myreader.get_system()

if not os.path.isfile(PROCESSDIR + "output_bonds"):
	print "Cannot execute output_bonds program. Please make sure to go to process_data/ directory and type make"
	sys.exit(1)

counter = 0

print 'DUPLEXMELT',str(20),beta,beta_min,beta_step,array_length
print '1.077 1.077 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64'


FRACTION_DS =  0.18
energy_fraction = (1.0-FRACTION_DS)/(1.0-FRACTION_DS+(9.0/beta)*FRACTION_DS);

while mysystem != False:
	mysystem.map_nucleotides_to_strands()
	launchargs = [PROCESSDIR + 'output_bonds',inputfile,conffile,str(counter)]
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	mysystem.read_all_interactions(myinput.stdout.readlines())

	for line in myinput.stderr.readlines():
      	  if "CRITICAL" in line:
              	  print line
	counter += 1
        print get_system_to_interaction_line(mysystem,energy_fraction)
	mysystem = myreader.get_system()
   	
