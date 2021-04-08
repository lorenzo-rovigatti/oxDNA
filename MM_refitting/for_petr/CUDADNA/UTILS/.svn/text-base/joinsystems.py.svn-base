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

'''
if len(sys.argv) < 4:
    base.Logger.log("Usage is %s configuration topology [howmany] [box_side]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()
'''

s2 = readers.LorenzoReader('prova.conf','prova.top').get_system()
s1 = readers.LorenzoReader('joined.dat', 'joined.top').get_system()

systems = [s1, s2]
base.Logger.log ("# Bringing nucleotides in box")
# bring the nucleotides back in the box
for S in systems:
    base.Logger.log ("   # Bringing nucleotides in box or another system...")
    for s in S._strands:
        diff = np.rint(s.cm_pos / S._box ) * S._box
        s.translate (-diff)

'''
base.Logger.log ("# Setting cdm to 0.")
for S in systems:
    base.Logger.log ("   # Setting cdm to 0. for a system")
    redo = True
    while redo:
        cdm = np.array([0.,0.,0.])
        for s in S._strands:
            for n in s._nucleotides:
                cdm += n.cm_pos
        cdm = (1. / float(S.get_N_Nucleotides())) * cdm
        S.translate (-cdm) 
        #S._prepare(visibility=None)
        cdm = np.array([0.,0.,0.])
        for s in S._strands:
            for n in s._nucleotides:
                cdm += n.cm_pos
        cdm = (1. / float(S.get_N_Nucleotides())) * cdm
        if np.dot (cdm, cdm) < 1.:
            redo = False
'''

# vthe biggest box has to be the first system
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

