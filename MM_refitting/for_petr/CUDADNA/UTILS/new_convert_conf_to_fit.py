#!/usr/bin/env python

#A utility that prints out the number of hydrogen bonds between different strands in the system 
import sys

sys.path.append('/usersVol2/sulc/MM_refitting/for_petr/CUDADNA/UTILS/')
sys.path.append('/usersVol2/sulc/MM_refitting/for_petr/CUDADNA/UTILS/process_data/')
#PROCESSDIR = '/usersVol2/sulc/SR_PMP_stuff/SR_PMP/UTILS/process_data/' 
PROCESSDIR = '/usersVol2/sulc/MM_refitting/for_petr/CUDADNA/UTILS/process_data/' 

import base
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import readers 
import subprocess

#STCK_BASE_NORM = 1.3448

def get_parameters_line(mod_file):
	#so far only assumes average parameters
	interesting_params = {'HYDR_EPS' : 1.077, 'STCK_BASE_EPS' : 1.3448 ,'STCK_FACT_EPS': 2.6568 , 'CRST_K' : 47.5}	
	input_vals = open(mod_file,'r')
	for line in input_vals.readlines():
		line = line.strip()
		if len(line) > 1 and line[0] != '#':
			key, val = line.split('=')
                        key = key.strip()
                        val = val.strip()
			if key in interesting_params.keys():
				interesting_params[key] = float(val)

	#now, a bit dirty solution:
	parline= ''
	for i in range(2):
		parline += str(interesting_params['HYDR_EPS']) + ' '
	for i in range(10):
		parline += str(interesting_params['STCK_BASE_EPS']) + ' ' 
	for i in range(10):
		parline += str(interesting_params['CRST_K']) + ' '

	parline += str(interesting_params['STCK_FACT_EPS'] / interesting_params['STCK_BASE_EPS'])


	return parline
	
def get_OPW_pairs(opfile,wfile):
	new_pairs = [] 
	input = open(opfile,'r')
	for line in input.readlines():
		if 'pair' in line:
			valA,valB = line.split('=')[1].split(',')
			valA = int(valA)
			valB = int(valB)
			new_pairs.append([valA,valB])

	input_w = open(wfile,'r')
	new_weights = [0] * (len(new_pairs) + 1)

	for i in range(len( new_weights)):
		line = input_w.readline()
		op, weight = line.split()
		op = int(op)
		weight = float(weight)
		new_weights[op] = weight
 		
	print new_pairs,new_weights
	return new_pairs, new_weights



def get_weight(mysystem,op_pairs, weights):
        op_value = 0
        for pair in op_pairs:
                nucA = mysystem._nucleotides[pair[0]]
                nucB = pair[1]
                if nucB in nucA.all_interactions[base.INT_HYDR].keys() and  nucA.all_interactions[base.INT_HYDR][nucB] < base.H_CUTOFF :
                        op_value += 1
       	if op_value >= len(weights):
		print 'Calculated op=',op_value, 'which is not compatible with ',weights 
        return weights[op_value]        

	
def get_nuc_in_strand_index(strand,global_index):
	for i in range(len(strand._nucleotides)):
		if strand._nucleotides[i].index == global_index:
			return i
	raise IOError('cannot convert nucid to strandid for '+ str(global_index) + ' in strand ' )  #if we got hwew, it means the nucleotiae is not in the strand 

def get_system_to_interaction_line(system,weight=1.):
	
	outline = str(weight) + ' ' 
	total_stack = 0
	H_bond_line = ''
	ST_line = ''
	CST_line = ''
	total_cross = 0
	
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
	Etotal = system._N*system.E_pot  # - (1.0 - energy_fraction) * total_stack

	if(  ST_line + H_bond_line == ''):
		ST_line = 'ST 0 0 0 1 0.00'
	
	outline =str(weight) + ' ' + str(Etotal) + ' ' + CST_line + ' '  + ST_line + ' ' + H_bond_line 
	return outline 


#the format of the fit file is:
#DUPLEXMELT <NO_OF_NUCLEOTIDES> <SIMULATION_BETA> <MIN_FIT_BETA> <FIT_BETA_STEP> <FIT_BETA_ARRAY_LENGTH>
#ORIGINAL PARAMETERS
#weight_of_state Energy_of_state(without entropic contribution in stacking) H strand1 nucid strandid2 nucid2 energy  ST ....



#PROCESSDIR = os.path.join("../process_data/")

if (len(sys.argv) < 4):
  print 'Usage %s input_file trajectory_file modoptions_file [order_parameter_file weight_file]' % sys.argv[0]
  sys.exit()


  
#now get topology file name:
inputfile = sys.argv[1]
conffile = sys.argv[2]
modfile = sys.argv[3]

if len(sys.argv) >= 5: 
	opfile = sys.argv[4]
	weightfile = sys.argv[5]
else:
	weightfile = None

outputfile = conffile + '.FIT'
print >> sys.stderr, ' Saving to ',outputfile

outfile = open(outputfile,'w')

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
beta_min = beta - 1.7
beta_step = 0.05
array_length = 60




myreader = readers.LorenzoReader(conffile,topologyfile)
mysystem = myreader.get_system()

if not os.path.isfile(PROCESSDIR + "output_bonds"):
	print "Cannot execute output_bonds program. Please make sure to go to process_data/ directory and type make"
	sys.exit(1)

counter = 0

print >> outfile, 'DUPLEXMELT',str(20),beta,beta_min,beta_step,array_length
#print '1.077 1.077 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64'
line = get_parameters_line(modfile)
print >> outfile, line
#print >> outfile, '1.077 1.077 1.077 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64 1.64 47.5 47.5 47.5 47.5 47.5 47.5 47.5 47.5 47.5 47.5 0'

FRACTION_DS =  0.18
#energy_fraction = (1.0-FRACTION_DS)/(1.0-FRACTION_DS+(9.0/beta)*FRACTION_DS);
#energy_fraction = 1.

if weightfile != None:
	ops, weights = get_OPW_pairs(opfile,weightfile)
	print >> sys.stderr, ops,weights
while mysystem != False:
	mysystem.map_nucleotides_to_strands()
	launchargs = [PROCESSDIR + 'output_bonds',inputfile,conffile,str(counter)]
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	mysystem.read_all_interactions(myinput.stdout.readlines())

	for line in myinput.stderr.readlines():
      	  if "CRITICAL" in line:
              	  print line
	counter += 1

	if weightfile != None:
		weight = get_weight(mysystem, ops,weights )
	else:
		weight = 1.
        #print get_system_to_interaction_line(mysystem,energy_fraction)
        print >> outfile, get_system_to_interaction_line(mysystem,weight)
	outfile.flush()
	
	mysystem = myreader.get_system()
   	
