#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import utils

if len(sys.argv) < 5:
	base.Logger.log("Usage is %s configuration1 topology1 configuration2 topology2 [configuration3 topology3] ... " % sys.argv[0], base.Logger.CRITICAL)
	sys.exit()

try:
	confs = sys.argv[1::2]
	tops = sys.argv[2::2]
except:
	base.Logger.log("Usage is %s configuration1 topology1 configuration2 topology2 [configuration3 topology3] ... " % sys.argv[0], base.Logger.CRITICAL)
	sys.exit(-1)

if len(tops) != len(confs):
	base.Logger.log ("There must be the same numbers of configurations and topologies", base.Logger.CRITICAL)
	sys.exit(-1)

systems = []
for i, conf in enumerate(confs):
	top = tops[i]
	base.Logger.log ("Working on %s %s" % (conf, top))
	systems.append (readers.LorenzoReader (conf, top).get_system())

assert (len(systems) == len(confs))

base.Logger.log ("# Bringing nucleotides in box")
# bring the nucleotides back in the box
for S in systems:
    base.Logger.log ("   # Bringing nucleotides in box or another system...")
    for s in S._strands:
        diff = np.rint(s.cm_pos / S._box ) * S._box
        s.translate (-diff)

# the biggest box has to be the first system, so we reshuffle them
has_max_box = 0
for S in systems[1:]:
    for k in [0, 1, 2]:
        if S._box[k] > systems[0]._box[k]:
            print >> sys.stderr, "the biggest box has to be the first system. Aborting"
            sys.exit (-2)

final = systems[0].copy()
njoined = 1
for S in systems[1:]:
    base.Logger.log ("# Trying to join")
    has_overlaps = True
    while has_overlaps:
        translated = S.copy()
        nold = final.N_strands
        translated.rotate(utils.get_random_rotation_matrix(), np.array([0.,0.,0.]))
        translated.translate(np.random.rand(3) * systems[0]._box)
        has_overlaps = False
        
        for s1 in translated._strands:
            if final.is_overlapping_better(s1):
                has_overlaps = True
                base.Logger.log ("   # Overlap, retrying")
                break

	#final = prova.copy()

    final = final.join(translated, systems[0]._box)
    njoined += 1
    base.Logger.log("   joined %i" %  njoined)

final.print_lorenzo_output ("joined.dat", "joined.top")
final.print_crepy_output ('joined.mgl', same_colors=True)

