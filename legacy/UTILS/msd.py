#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys

block_size = 10

if len(sys.argv) < 3:
    base.Logger.log("Usage is %s trajectory topology [output]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
stmp = l.get_system()
nconf, nblock = 0, 0
msd = np.zeros(block_size, dtype = np.float64)
nmsd = np.zeros(block_size, dtype = np.int64)
while stmp:
    if (nconf%block_size == 0):
        nblock += 1
        print >> sys.stderr, " -- Starting block %i " % (nblock)
        s = []

    nconf += 1
    print >> sys.stderr, " --- Working on conf %d, block %i" % (nconf, nblock)

    for i in range(len(s)):
        for j in range(i):
            #print "###", i, j
            strands1 = s[i]._strands
            strands2 = s[j]._strands
            
            #if len(strands1) != len(strands2):
            #    print >> sys.stdout, "## Disaster: different number of strands between confs. Aborting"
            #    sys.exit(-2)

            dt = i - j
            for k in range(len(strands1)):
                # we now do the stuff
                r1 = strands1[k].cm_pos
                r2 = strands2[k].cm_pos
                dr = r1 - r2 # newer minus older
                msd[dt] += np.dot(dr,dr)
                nmsd[dt] += 1

    ''' 
    # Now we work on the msd...
    for strand1 in stmp._strands:
        myindex = strand1.index
        print stmp._N, len(stmp._strands)
        print "# index: %i" % (myindex)
        r1 = strand1.cm_pos
        for oldsystem in s:
            print "## index: %i" % (myindex)
            strand2 = oldsystem._strands[myindex-1]
            r2 = strand2.cm_pos
            dr = r2 - r1
            '''

    # we append the configuration we just worked on
    s.append(stmp)
    
    # try to get the next one from the trajectory.
    # if l.get_system() fails, the while cycle finishes
    stmp = l.get_system()

#msd = msd / nmsd # numpy works element-wise; (convenient)

out = open ("msd.dat", "w")
for i in range(1, block_size - 1):
    print >> out, "%7i %10.5g" % (i, msd[i] / nmsd[i])
out.close ()

print >> sys.stderr, "All DONE"

