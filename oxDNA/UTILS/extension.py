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
    base.Logger.log("Usage is %s configuration topology external_forces_file" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

traps = forces.parse_traps (sys.argv[3])
l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()
'''
niter = 1
while s:
    #base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
    r1 = s._nucleotides[n1].cm_pos
    r2 = s._nucleotides[n2].cm_pos
    #dr = s._nucleotides[n1].distance (s._nucleotides[n2], PBC=True, box=s._box)
    dr = s._nucleotides[n1].distance (s._nucleotides[n2], PBC=False)
    print s._time, np.sqrt(np.dot(dr, dr)), dr[0], dr[1], dr[2]
    s = l.get_system()
    niter += 1
'''

if len(traps) < 2:
    print >> sys.stderr, "there are less than 2 traps. Aborting..."
    sys.exit (-1)

if traps[0].type != traps[1].type:
    print >> sys.stderr, "trap types don't match. Aborting..."
    sys.exit (-2)

i1 = traps[0].particle
i2 = traps[1].particle

niter = 1
while s:
    #base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
    n1 = s._nucleotides[i1]
    n2 = s._nucleotides[i2]
    dr = n1.distance (n2, PBC=False)
    #dr = s._nucleotides[n1].distance (s._nucleotides[n2], PBC=False)
    if traps[0].type == 'string':
        print s._time, np.sqrt(np.dot(dr, dr)), dr[0], dr[1], dr[2]
    elif traps[0].type == 'trap':
        f1 = traps[0].get_force (n1.cm_pos)
        f2 = traps[1].get_force (n2.cm_pos)
        print s._time, np.sqrt(np.dot(dr, dr)), f1[0], f1[1], f1[2], f2[0], f2[1], f2[2] 
    else:
        pass
    s = l.get_system()
    niter += 1


