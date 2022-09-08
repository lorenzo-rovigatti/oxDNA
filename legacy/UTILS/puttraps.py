#!/usr/bin/env python

import sys
import random as rnd
import base, readers

n_in_tetramer = 132
lambda_cdm = 0.1
lambda_patch = 0.1

if not len(sys.argv) >= 3:
	print >> sys.stderr, "Usage:", sys.argv[0], "configuration", "topology", "lambda_cmd", "lambda_patch"
	sys.exit(2)

try:
	lambda_cdm = float(sys.argv[3])
	lambda_patch = float(sys.argv[4])
except:
	print >> sys.stderr, "Could not read in the lambdas..."
	sys.exit(3)

try:
	print >> sys.stderr, "Reading system", sys.argv[1], sys.argv[2]
	r = readers.LorenzoReader (sys.argv[1], sys.argv[2])
	s = r.get_system()
except:
	print >> sys.stderr, "Could not read in ", sys.argv[1], sys.argv[2]
	sys.exit(-2)


print >> sys.stderr, "Assuming n_in_tetramer = 132"

forces = open ("ext_einst.dat", 'w')

# now for each tetramer we set the traps where they are supposed to be
# cdm: 12th particle of first strand
for i, n in enumerate(s.get_nucleotide_list()):
	if i%n_in_tetramer == 12:
		#print n.index, n.cm_pos
		print >> forces, "{\ntype=trap\nparticle=%d\nstiff=%g\npos0=%g, %g, %g\nrate=0.\ndir=0., 0., 0.\ngroup_name=coms\n}\n\n" % (n.index, lambda_cdm, n.cm_pos[0], n.cm_pos[1], n.cm_pos[2])

# for each "patch":
patches = [30, 30+33, 30 + 66, 30 + 99]
print >> sys.stderr, "patch indexes: ", patches
for i, n in enumerate(s.get_nucleotide_list()):
	if i%n_in_tetramer in patches:
		#print n.index, n.cm_pos
		print >> forces, "{\ntype=trap\nparticle=%d\nstiff=%g\npos0=%g, %g, %g\nrate=0.\ndir=0., 0., 0.\ngroup_name=patches\n}\n\n" % (n.index, lambda_patch, n.cm_pos[0], n.cm_pos[1], n.cm_pos[2])

forces.close()

print >> sys.stderr, "External forces file ext_einst.dat printed"

