#!/usr/bin/env python

import base
import readers

import sys




if len(sys.argv) != 4:
	print 'Usage: %s conf topology pairs_file' % (sys.argv[0])
	sys.exit(1)

r = readers.LorenzoReader(sys.argv[1],sys.argv[2])
s = r.get_system()

infile = open(sys.argv[3],'r')

index = 16

new_topology = range(s._N+26,26+s._N+s._N)

while True:
	line1 = infile.readline()
	line2 = infile.readline()
	if not line2: break

	if(len(line1.split()) > 1 and len(line2.split()) > 1):
		a = int(line1.split()[0])
		b = int(line1.split()[1])
		A = int(line2.split()[0])
		B = int(line2.split()[1])
		reindexed = b
		for i in range(a,A+1):
			otype = s._nucleotides[i]._btype
			ntype = index + otype
			index += 4
			comptype = 3 - ntype
			new_topology[i] = ntype
			if(otype + s._nucleotides[reindexed]._btype != 3):
				print 'Error, base ',i,'not complementary to base ',reindexed
			new_topology[reindexed] = comptype

			reindexed -= 1	
								
#for i in range(len(new_topology)):
#	print i,new_topology[i]


final_topology = "%d %d\n" % (len(s._nucleotides), len(s._strands))
topology = ''
for s in s._strands:
    sc, st = s.get_output(base.OUT_LORENZO)
    topology += st

lines = topology.split('\n')

for i in range(len(lines)):
	vals = lines[i].split()
	if(i >= len(new_topology)):
		pass
		#print >> sys.stderr, 'impossible index ',i
	else:
		vals[1] = str(new_topology[i])
		for info in vals:
			final_topology += info + ' '
		final_topology += '\n'


print final_topology,
