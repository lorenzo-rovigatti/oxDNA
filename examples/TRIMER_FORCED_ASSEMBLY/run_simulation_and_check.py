#!/usr/bin/env python
import sys

PROCESSDIR = '../../UTILS/'

import subprocess
import base
import readers

import sys

def check_domain_status(conffile,topologyfile,a1,b1,a2,b2):
	r = readers.LorenzoReader(conffile,topologyfile)
	mysystem = r.get_system()
	mysystem.map_nucleotides_to_strands()
	print PROCESSDIR+'output_bonds.py'
	launchargs = [PROCESSDIR+'output_bonds.py',inputfile,conffile]
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	mysystem.read_H_bonds(myinput.stdout.readlines())
	#a1 paried with b1, a2 with b2
	print 'Launched my inpit and got' , inputfile, conffile
	bindex = b1
	totbonds = 0
	for i in range(a1,a2+1):
		if bindex in mysystem._nucleotides[i].interactions:
			totbonds += 1	
		bindex -= 1
	return totbonds, a2-a1+1




def gen_force_file(conffile,topology,bond_pairs,no_of_w,out_file_name):
	r = readers.LorenzoReader(conffile,topology)
	s = r.get_system()
	counter = 0
	outfile = open(out_file_name,'w')
	for pair in bond_pairs:
			a = pair[0]
			b = pair[1]
			if(s._nucleotides[a]._btype + s._nucleotides[b]._btype != 3):
				print 'Error, bases ',a,b,' are not complementary'
			ous = '{ \ntype = mutual_trap\nparticle = %d\nstiff = 0.9\nr0 = 1.2\nref_particle = %d\nPBC=1\n}\n' % (a,b)
			ous += '{ \ntype = mutual_trap\nparticle = %d\nstiff = 0.9\nr0 = 1.2\nref_particle = %d\nPBC=1\n}\n' % (b,a)
			outfile.write(ous)
			counter = counter + 1
			print 'Using mutual trap ', a,b
			if(counter > no_of_w):
				outfile.close()
				return
	outfile.close()


if len(sys.argv) != 4:
	print 'Usage: %s conf topology pairs_file' % (sys.argv[0])
	sys.exit(1)

conffile = sys.argv[1]
topologyfile = sys.argv[2]

import subprocess

bond_pairs = []
infile = open(sys.argv[3])
for line in infile.readlines():
		if(len(line.split()) > 1):
			a = int(line.split()[0])
			b = int(line.split()[1])
			bond_pairs.append( [a,b])

#print bond_pairs
counter = 1
finished = False
while not finished:
	print 'Generating forces with adding  traps  ',counter+1, ' for domain ', (counter-1)/2 + 1
	sys.stdout.flush()
	gen_force_file(sys.argv[1],sys.argv[2],bond_pairs,counter,'partial_forces.txt')
	inputfile = 'input'
	launchargs = ['./oxDNA',inputfile]
	myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
	myinput.wait()
	
	
	totbonds,maxbonds = check_domain_status(conffile,topologyfile,bond_pairs[counter - 1][0],bond_pairs[counter - 1][1], bond_pairs[counter][0],bond_pairs[counter ][1])
	print 'One itetation finished, for domain %d,  we had %d bonds out of %d  ' % ((counter-1)/2  + 1,totbonds,maxbonds)
	if(totbonds / float(maxbonds) > 0.5):	
		counter += 2
		print 'Making step for the next domain'
	else:
		print 'The bonds still not fully formed for the domain, continuing...'

	if(counter >= len(bond_pairs)):
		print 'Done, all iterations finished'	
		finished = True	
	sys.stdout.flush()
