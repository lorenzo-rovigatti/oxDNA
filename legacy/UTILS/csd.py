#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
from utils import get_angle
import energies

MAX_NEIGHS = 2

# this incapsulates a double strand
class CSDParticle(object):
    def __init__(self, s1, s2):
        self.index = s1.index/2
        self.cluster = self.index
        self.nucleotides = s1._nucleotides + s2._nucleotides
        if len(self.nucleotides) != 2 and len(self.nucleotides) != 4:
            base.Logger.die("LorenzoReader.get_system() must be called with only_strand_ends")
        self.neighs = []

    def add_neigh(self, q):
        if q not in self.neighs: self.neighs.append(q)

def is_coaxial(p, q, box):
    dr = p.pos_stack - q.pos_stack
    dr -= box * np.rint (dr / box)
    dist = np.sqrt(np.dot(dr, dr))
    if dist < base.CXST_RCHIGH and dist > base.CXST_RCLOW:
        t4mt0 = get_angle(p._a3, q._a3) - base.CXST_THETA4_T0
        if t4mt0 < base.CXST_THETA4_TC and t4mt0 > -base.CXST_THETA4_TC:
            t1mt0 = np.pi - get_angle(p._a1, q._a1) - base.CXST_THETA1_T0
            t1mt02pi = 2 * np.pi - t1mt0
            if (t1mt0 < base.CXST_THETA1_TC and t1mt0 > -base.CXST_THETA1_TC) or (t1mt02pi < base.CXST_THETA1_TC and t1mt02pi > -base.CXST_THETA1_TC):
                drn = dr / dist
                t5mt0 = get_angle(q._a3, drn) - base.CXST_THETA5_T0
                t5mt0pi = np.pi - t5mt0
                if (t5mt0 < base.CXST_THETA5_TC and t5mt0 > -base.CXST_THETA5_TC) or (t5mt0pi < base.CXST_THETA5_TC and t5mt0pi > -base.CXST_THETA5_TC):
                    t6mt0 = get_angle(-p._a3, drn) - base.CXST_THETA6_T0
                    t6mt0pi = np.pi - t6mt0
                    if (t6mt0 < base.CXST_THETA6_TC and t6mt0 > -base.CXST_THETA6_TC) or (t6mt0pi < base.CXST_THETA6_TC and t6mt0pi > -base.CXST_THETA6_TC):
                        rback = p.pos_back - q.pos_back
                        rback -= box * np.rint (rback / box)
                        rback /= np.sqrt(np.dot(rback, rback))
                        cosphi3 = np.dot(drn, np.cross(rback, q._a1))
                        if cosphi3 > base.CXST_PHI3_XC: return True

    return False

def flip_neighs(particles, i, clusters):
    p = particles[i]
    for n in p.neighs:
        if n.cluster > p.cluster:
            clusters[n.cluster] -= 1
            n.cluster = p.cluster
            clusters[p.cluster] += 1

            flip_neighs(particles, n.index, clusters)

def get_csd(syst):
    syst.do_cells()
    N = len(syst._strands) / 2
    particles = [CSDParticle(syst._strands[i], syst._strands[i+1]) for i in range(0, 2*N, 2)]

    kc = [0, 0, 0]
    for p in particles:
        for n1 in p.nucleotides:
            # compute n1's cell index
            cs = np.array(np.floor((n1.cm_pos/syst._box - np.rint(n1.cm_pos / syst._box ) + 0.5) * (1.-base.FLT_EPSILON) * syst._box / syst._cellsides), dtype=int)
            for kc[0] in (-1, 0, 1):
                for kc[1] in (-1, 0, 1):
                    for kc[2] in (-1, 0, 1):
                        nc = [(cs[i] + kc[i] + syst._N_cells[i]) % syst._N_cells[i] for i in range(3)]
                        cella2 = nc[0] + nc[1] * syst._N_cells[0] + nc[2] * syst._N_cells[0] * syst._N_cells[1]
                        # loop over particles in cella2
                        n2 = syst._head[cella2]
                        while n2:
                            #if n2.strand > n1.strand and (is_coaxial(n1, n2, syst._box) or energies.cross_stacking(n1, n2, syst._box) < 0.):
                            if n2.strand > n1.strand and is_coaxial(n1, n2, syst._box):
                                q = particles[n2.strand/2]
                                p.add_neigh(q)
                                q.add_neigh(p)
                            
                            n2 = n2.next

    clusters = [1, ] * N
    csd = [0, ] * (N+1)
    
    for i in range(N):
        flip_neighs(particles, i, clusters)

    for i in range(N):
        csd[clusters[i]] += 1
    
    return csd

if len(sys.argv) < 3:
    base.Logger.log("Usage is %s configuration topology [N_skip]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(sys.argv) > 3:
    N_skip = int(sys.argv[3])
else: N_skip = 0

output = sys.argv[1] + ".csd"
E_output = output + ".E"
    
l = readers.LorenzoReader(sys.argv[1], sys.argv[2])
s = l.get_system(only_strand_ends=True, N_skip=N_skip)

counter = 0
N = len(s._strands) / 2
csd = [0, ] * (N+1)
E_bonds = [0, ] * (N+1)
E_N = [0, ] * (N+1)

# cycle over the trajectory
while s:
    ncsd = get_csd(s)
    nbonds = 0
    for i in range(1, N+1): 
        csd[i] += ncsd[i]
        nbonds += ncsd[i] * (i-1)

    E_bonds[nbonds] += s.E_pot
    E_N[nbonds] += 1

    base.Logger.log("Analyzed configuration", base.Logger.INFO)
    s = l.get_system(only_strand_ends=True)
    counter += 1

f = open (E_output, "w")
for i in range(0, len(E_bonds)):
    if E_N[i] > 0: f.write(str(i) + " " + str(E_bonds[i]/float(E_N[i])) + "\n")
f.close()

f = open (output, "w")
for i in range(1, len(csd)):
    if csd[i] > 0: f.write(str(i) + " " + str(csd[i]/float(counter)) + "\n")
f.close()

base.Logger.log("Output printed on '%s' and '%s'" % (output, E_output), base.Logger.INFO)

