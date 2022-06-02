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

if len(sys.argv) < 3:
    base.Logger.log("Usage is %s configuration topology" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()

try:
    i1 = 0
    i2 = len(s._strands[0]._nucleotides) - 1
    print >> sys.stderr, "Nucleotides", i1, i2
except:
    print >> sys.stderr, "Supply nucleotides... Aborting"
    sys.exit (-1)

L2 = 0.
l0 = 0.
Ll0 = 0.
Lmax = 1.025 * (i2 - i1 + 1)
niter = 1
while s:
    base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
    first = s._strands[0]._nucleotides[0]
    second = s._strands[0]._nucleotides[1]
    last = s._strands[0]._nucleotides[-1]
    
    r0N = first.distance (last, PBC=False)
    r01 = first.distance (second, PBC=False)
    
    l0 += np.sqrt (np.dot (r01, r01))
    L2 += np.dot (r0N, r0N)
    Ll0 += np.dot (r01, r0N)
    
    s = l.get_system()
    niter += 1

Ll0 /= float (niter)
L2 /= float (niter)
l0 /= float (niter)

Pl = Ll0 / l0
Kl = L2 / Lmax

print Pl, Kl, l0, L2, Ll0, "Pl, Kl, <l0>, <L2>, <L * l0>"


