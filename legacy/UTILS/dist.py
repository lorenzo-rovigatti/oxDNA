#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import external_forces as forces 

if len(sys.argv) < 4:
    base.Logger.log("Usage is %s configuration id1 id2" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()

try:
    i1 = int(sys.argv[3])
    i2 = int(sys.argv[4])
except:
    print >> sys.stderr, "Supply nucleotides... Aborting"
    sys.exit (-1)

try:
    site = int (sys.argv[5])
except:
    site = None

niter = 1
while s:
    #base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
    n1 = s._nucleotides[i1]
    n2 = s._nucleotides[i2]
    dr = n1.distance (n2, PBC=True, box=s._box)
    if site is not None:
        if site == -1:
            dr = n1.cm_pos - n2.cm_pos
        elif site == 1:
            dr = n1.pos_stack - n2.pos_stack
        elif site == 2:
            dr = n1.pos_base - n2.pos_base
        elif site == 0:
            pass
        else:
            raise KeyError
        dr -= s._box * np.rint (dr / s._box)

    print s._time, np.sqrt(np.dot(dr, dr)), dr[0], dr[1], dr[2]
    s = l.get_system()
    niter += 1


