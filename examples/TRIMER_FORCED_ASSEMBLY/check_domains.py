#!/usr/bin/env python

import base
import readers

import sys

PROCESSDIR = '../../UTILS/' 

def check_domain_status(mysystem,a1,b1,a2,b2):
	#a1 paried with b1, a2 with b2
	bindex = b1
	totbonds = 0
	for i in range(a1,a2+1):
		if bindex in mysystem._nucleotides[i].interactions:
			totbonds += 1	
		bindex -= 1
	return totbonds, a2-a1+1




	
if len(sys.argv) != 5:
	print 'Usage: %s input conf topology pairs_file' % (sys.argv[0])
	sys.exit(1)

inputfile = sys.argv[1]
conffile = sys.argv[2]
topology_file = sys.argv[3]
r = readers.LorenzoReader(sys.argv[2],sys.argv[3])
mysystem = s = r.get_system()

infile = open(sys.argv[4],'r')

index = 16

new_topology = range(s._N)
import subprocess
mysystem.map_nucleotides_to_strands()
launchargs = [PROCESSDIR + 'output_bonds.py',inputfile,conffile]
myinput = subprocess.Popen(launchargs,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
mysystem.read_H_bonds(myinput.stdout.readlines())

for line in myinput.stderr.readlines():
  if "CRITICAL" in line:
	  print line

domains = 0

while True:
	line1 = infile.readline()
	line2 = infile.readline()
	if not line2: break

	if(len(line1.split()) > 1 and len(line2.split()) > 1):
		a = int(line1.split()[0])
		b = int(line1.split()[1])
		A = int(line2.split()[0])
		B = int(line2.split()[1])
		tot,out_of = check_domain_status(s,a,b,A,B)
		print domains, ': ',tot,'/',out_of
		domains +=1	
					
#for i in range(len(new_topology)):
#	print i,new_topology[i]

