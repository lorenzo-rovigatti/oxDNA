#!/usr/bin/env python

import sys
try:
    import numpy as np
except:
    import mynumpy as np
from base import *
from generators import *

side = 57

g = StrandGenerator ()

s = System ([side, side, side])
sa = "CATAACACAATCACATCTCATCCATTCCACTCA"
for i in range (54):
    cdm = np.random.random_sample(3) * s._box
    axis = np.random.random_sample(3)
    axis /= np.sqrt(np.dot(axis, axis))
    print >> sys.stderr, "Adding single strand"
    seq = [base.base_to_number[x] for x in sa]
    while not s.add_strand(g.generate(len(sa), sequence=seq, dir=axis, start_pos=cdm, double=False)):
        cdm = np.random.random_sample(3) * s._box
        axis = np.random.random_sample(3)
        axis /= np.sqrt(np.dot(axis, axis))
        print >> sys.stderr, "  riprovo %i" % (i)
        print >> sys.stderr, "  done %i" % (i)

sb = "CATAACACAATCACATCTCACCACCAAACTTCA"
for i in range (6):
    cdm = np.random.random_sample(3) * s._box
    axis = np.random.random_sample(3)
    axis /= np.sqrt(np.dot(axis, axis))
    print >> sys.stderr, "Adding single strand"
    seq = [base.base_to_number[x] for x in sb]
    while not s.add_strand(g.generate(len(sb), sequence=seq, dir=axis, start_pos=cdm, double=False)):
        cdm = np.random.random_sample(3) * s._box
        axis = np.random.random_sample(3)
        axis /= np.sqrt(np.dot(axis, axis))
        print >> sys.stderr, "  riprovo %i" % (i)
        print >> sys.stderr, "  done %i" % (i)

sa = "CACTAACATACAACATCTCATAACACAATCACA"
sb =         "TGAGATGTGATTGTGTTATGAGATG"
for i in range (120):
    cdm = np.random.random_sample(3) * s._box
    axis = np.random.random_sample(3)
    axis /= np.sqrt(np.dot(axis, axis))
    print >> sys.stderr, "Adding double"
    while not s.add_strands (g.generate_double_offset(seqA=sa, seqB=sb, offset=13, dir=axis, start_pos=cdm)):
        cdm = np.random.random_sample(3) * s._box
        axis = np.random.random_sample(3)
        axis /= np.sqrt(np.dot(axis, axis))
        print >> sys.stderr, "  riprovo %i" % (i)
        print >> sys.stderr, "  done %i" % (i)

#s.renumber ()
s.print_vmd_xyz_output ("pallinz.xyz", same_colors=True)
s.print_lorenzo_output ("riprova.conf", "riprova.top")

