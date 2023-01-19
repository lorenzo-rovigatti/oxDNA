#!/usr/bin/env python
#MAKE SURE THAT THE FULL PATH TO UTILS DIRECTORY IS IN PYTHONPATH

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import math

prune = 1

def get_lj(s,nucid):
    #this function returns normalized vector pointing from base midpoint of nucid-th base pair to (nudic+1)th midpoint of the next bp
    box = s._box
    i1A = nucid
    i1B = len(s._strands[1]._nucleotides)-1-nucid
    firstA = s._strands[0]._nucleotides[i1A]
    firstB = s._strands[1]._nucleotides[i1B]     
    secondA = s._strands[0]._nucleotides[i1A+1]
    secondB = s._strands[1]._nucleotides[i1B-1]
    first_midpos = (firstA.get_pos_base() + firstB.get_pos_base()) / 2.0  
    second_midpos = (secondA.get_pos_base() + secondB.get_pos_base()) / 2.0  
    lj = second_midpos - first_midpos         
    lj -= box * np.rint(lj / box)
    return lj

if len(sys.argv) < 5:
    base.Logger.log("Usage is %s trajectory_file topology nucid_1  nucid_2 [prune]" % sys.argv[0], base.Logger.CRITICAL)
    print('The program calculates Kuhn length, Persistence length, end mean of nucid_1 to nucid_2 distance; the optional argument prune is a number of configurations to skip in the trajectory file')
    sys.exit()

l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system()

nid1 = int(sys.argv[3])
nid2 = int(sys.argv[4])

if len(sys.argv) == 6 :
    prune = int(sys.argv[5])
try:
    i1A = nid1
    i2A = nid2
    i1B = len(s._strands[1]._nucleotides) - nid1 - 1
    i2B = len(s._strands[1]._nucleotides) - nid2 - 1 

    #print >> sys.stderr, "Nucleotides", i1A, i2A, i1B, i2B
    print("#Nucleotides", i1A, i2A, i1B, i2B)

except:
    print("Supply nucleotides... Aborting")
    sys.exit(1)

L2 = 0.
l0 = 0.
Ll0 = 0.
Lmax = 1.025 * (i2A - i1A + 1)
niter = 0
correlations = [0.] * (i2A -i1A + 1)

correlations_counter = [0] * (i2A -i1A + 1)

read_confs = 1

while s:
    if(read_confs % prune != 0):
        read_confs += 1
        s = l.get_system()
        continue

    #base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
    for subcounter in range(0,i2A-i1A):
        firstA = s._strands[0]._nucleotides[i1A+subcounter]
        firstB = s._strands[1]._nucleotides[i1B-subcounter]     
        secondA = s._strands[0]._nucleotides[i1A+1+subcounter]
        secondB = s._strands[1]._nucleotides[i1B-1-subcounter]
         
        first_midpos = (firstA.get_pos_base() + firstB.get_pos_base()) / 2.0  
        second_midpos = (secondA.get_pos_base() + secondB.get_pos_base()) / 2.0  
      
        box = s._box  
        r01 = second_midpos - first_midpos 
        r01 -= box * np.rint (r01 / box)
        for j in range (0,i2A - i1A-subcounter):
            jA = j + i1A + subcounter
            #jB = len(s._strands[1]._nucleotides) - jA - 1
            #lastjA = s._strands[0]._nucleotides[jA]
            #lastjB = s._strands[1]._nucleotides[jB]
            #last_midpos_j = (lastjA.get_pos_base() + lastjB.get_pos_base()) / 2.0  
            #r0j = last_midpos_j - first_midpos
            #r0j -= box * np.rint(r0j / box)
            lj = get_lj(s,jA)
            ml0 = r01 / math.sqrt(np.dot(r01,r01))
            lj = lj / math.sqrt(np.dot(lj,lj))
            correlations[j] += np.dot(ml0,lj)
            #if j == 1:
            #    print 'Adding ', np.dot(ml0,lj)
            correlations_counter[j] += 1.
    #r0N = first.distance (last, PBC=False)
    #r01 = first.distance (second, PBC=False)
    
    l0 += np.sqrt (np.dot (r01, r01))
    
    s = l.get_system()
    niter += 1
    read_confs += 1

l0 /= float (niter)


print('#Configurations in total',niter)
print('# ',sys.argv[1])
print('# Correlation between lk and l_0:')
for j in range(len(correlations)-1):
    print(j,correlations[j]/correlations_counter[j])

