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

if len(sys.argv) < 4:
    base.Logger.log("Usage is %s configuration topology [howmany] [box_side]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

conf = sys.argv[1]
top = sys.argv[2]
howmany = 2
box_side = 100.
if len (sys.argv) > 3:
    howmany = int(sys.argv[3])
if len (sys.argv) > 4:
    box_side = float(sys.argv[4])

systems = []
for i in range(0, howmany):
    base.Logger.log ("# Reading system %i" % i)
    l1 = readers.LorenzoReader(conf, top)
    systems.append (l1.get_system())

base.Logger.log ("# Bringing nucleotides in box")
# bring the nucleotides back in the box
for S in systems:

    myr = S._strands[0].cm_pos
    for s in S._strands:
        s.translate (-myr)

    base.Logger.log ("   # Bringing nucleotides in box or another system...")
    for s in S._strands:
        diff = np.rint(s.cm_pos / S._box ) * S._box
        s.translate (-diff)

base.Logger.log ("# Setting cdm to 0.")
for S in systems:
    base.Logger.log ("   # Setting cdm to 0. for system another system")

    redo = True
    while redo:
        cdm = np.array([0.,0.,0.])
        for s in S._strands:
            for n in s._nucleotides:
                cdm += n.cm_pos
        cdm = (1. / float(S.get_N_Nucleotides())) * cdm
        S.translate (-cdm) 
        #S._prepare(visibility=None)
        redo = False

        cdm = np.array([0.,0.,0.])
        for s in S._strands:
            for n in s._nucleotides:
                cdm += n.cm_pos
        cdm = (1. / float(S.get_N_Nucleotides())) * cdm
        if np.dot (cdm, cdm) < 1.:
            redo = False

final_box = np.array([box_side, box_side, box_side])
final = systems[0].copy()

final.print_lorenzo_output ("justone.dat", "justone.top")

njoined = 1
for S in systems[1:]:
    base.Logger.log ("# Trying to join")
    has_overlaps = True
    while has_overlaps:
        translated = S.copy()
        nold = final.N_strands
        translated.rotate(utils.get_random_rotation_matrix(), np.array([0.,0.,0.]))
        translated.translate(np.random.rand(3) * final._box)
        has_overlaps = False
        
        for s1 in translated._strands:
            if final.is_overlapping_better(s1):
                has_overlaps = True
                base.Logger.log ("   # Overlap, retrying")
                break

        '''prova = final.join(translated, final_box)
        for s1 in prova._strands[nold:]:
            for s2 in prova._strands[:nold]:
                if s1.overlaps_with(s2, final_box):
                    has_overlaps = True
                    base.Logger.log ("   # Overlap, retrying")
                    break'''
		
	#final = prova.copy()

    final = final.join(translated, final_box)
    njoined += 1
    base.Logger.log("   joined %i" %  njoined)

final.print_lorenzo_output ("joined.dat", "joined.top")
final.print_crepy_output ('joined.mgl', same_colors=True)
