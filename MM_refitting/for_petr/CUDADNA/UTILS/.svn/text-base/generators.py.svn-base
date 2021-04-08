import base
try:
    import numpy as np
except:
    import mynumpy as np
import math
import sys
import utils

# need to add this to avoid stupid warning about invalid division at line dir /= np.sqrt(np.dot(dir,dir))
try:
    np.seterr(invalid='ignore')
except:
    pass

BP = "bp"
DEGREES = "degrees"
DEG = "degrees"
RAD = "radiants"
RADIANTS = "radiants"
BASE_PAIRS = "bp"

class StrandGenerator (object):
    def generate(self, bp, sequence=None, start_pos=np.array([0, 0, 0]), dir=np.array([0, 0, 1]), perp=False, double=True, rot=0.):
        # we need numpy array for these
        start_pos = np.array(start_pos, dtype=float)
        dir = np.array(dir, dtype=float)
        if sequence == None:
            sequence = np.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence += np.random.randint(0, 4, n)
            base.Logger.log("sequence is too short, adding %d random bases" % n, Logger.WARNING)
            
        # create the sequence of the second strand as made of complementary bases
        sequence2 = [3-s for s in sequence]
        sequence2.reverse()
        
        # we need to find a vector orthogonal to dir
        dir_norm = np.sqrt(np.dot(dir,dir))
        if dir_norm < 1e-10:
            base.Logger.log("direction must be a valid vector, defaulting to (0, 0, 1)", base.Logger.WARNING)
            dir = np.array([0, 0, 1])
        else: dir /= dir_norm

        if perp is None or perp is False:
            v1 = np.random.random_sample(3)
            v1 -= dir * (np.dot(dir, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perp;
        
        # and we need to generate a rotational matrix
        R0 = utils.get_rotation_matrix(dir, rot)
        #R = get_rotation_matrix(dir, np.deg2rad(35.9))
        R = utils.get_rotation_matrix(dir, [1, BP])
        
        ns1 = base.Strand()

        a1 = v1
        a1 = np.dot (R0, a1)
        rb = np.array(start_pos)
        a3 = dir
        for i in range(bp):
            ns1.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence[i]))
            if i != bp-1:
                a1 = np.dot(R, a1)
                rb += a3 * base.BASE_BASE
            
        if double == True:
            a1 = -a1
            a3 = -dir
            R = R.transpose()
            ns2 = base.Strand()
            for i in range(bp):
                ns2.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence2[i]))
                a1 = np.dot(R, a1)
                rb += a3 * base.BASE_BASE
                
            return ns1, ns2
        else: return ns1

    def generate_or_sq(self, bp, sequence=None, start_pos=np.array([0,0,0]), dir=np.array([0, 0, 1]), perp=None, double=True, rot=0., angle=np.pi/180*33.75, length_change=0, region_begin=0, region_end=0):
        # we need numpy array for these
        start_pos = np.array(start_pos, dtype=float)
        dir = np.array(dir, dtype=float)
        if sequence == None:
            sequence = np.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence += np.random.randint(0, 4, n)
            base.Logger.log("sequence is too short, adding %d random bases" % n, base.Logger.WARNING)
        # angle should be an array, with a length 1 less than the # of base pairs
        if not isinstance(angle, np.ndarray):
            angle = np.ones(bp) * angle
        elif len(angle) != bp - 1:
            base.Logger.log("generate_or_sq: incorrect angle array length, should be 1 less than number of base pairs", base.Logger.CRITICAL)
        # create the sequence of the second strand as made of complementary bases
        sequence2 = [3-s for s in sequence]
        sequence2.reverse()
        
        # we need to find a vector orthogonal to dir
        dir_norm = np.sqrt(np.dot(dir,dir))
        if dir_norm < 1e-10:
            base.Logger.log("direction must be a valid vector, defaulting to (0, 0, 1)", base.Logger.WARNING)
            dir = np.array([0, 0, 1])
        else: dir /= dir_norm
        
        if perp is None:
            v1 = np.random.random_sample(3)
            v1 -= dir * (np.dot(dir, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perp;

        # and we need to generate a rotational matrix
        R0 = utils.get_rotation_matrix(dir, rot)
#        R = utils.get_rotation_matrix(dir, angle)
        #R = get_rotation_matrix(dir, [1, BP])
        
        ns1 = base.Strand()
        a1 = v1
        a1 = np.dot (R0, a1)
        rb = np.array(start_pos)
        a3 = dir
        for i in range(bp):
            ns1.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence[i]))
            if i != bp-1:
                R = utils.get_rotation_matrix(dir, angle[i])
                a1 = np.dot(R, a1)
                rb += a3 * base.BASE_BASE
                if length_change:
                    for j in range(len(length_change)):
                        if i >= region_begin[j] and i < region_end[j]:
                            if length_change[j]:
                                rb += a3 * base.BASE_BASE * (- (float(length_change[j])/(region_end[j] - region_begin[j])))
            
        if double == True:

            angle = np.flipud(angle)
            a1 = -a1
            a3 = -dir
            ns2 = base.Strand()
            
            for i in range(bp):
                
                ns2.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence2[i]))
                if i != bp - 1:
                    R = utils.get_rotation_matrix(dir,angle[i]).transpose()
                    a1 = np.dot(R, a1)
                    rb += a3 * base.BASE_BASE
                    if length_change:
                        for j in range(len(length_change)):
                            if bp - 2 - i >= region_begin[j] and bp - 2 - i < region_end[j]:
                                if length_change[j]:
                                    rb += a3 * base.BASE_BASE * (- (float(length_change[j])/(region_end[j] - region_begin[j])))
                                        
            return ns1, ns2
        else: return ns1

    def generate_double_offset(self, seqA, seqB, offset, start_pos=np.array([0,0,0]), dir=np.array([0,0,1]), perp=None, rot=0):
        if isinstance (seqA, str):
            seqa = [base.base_to_number[x] for x in seqA]
        else:
            seqa = seqA
        if isinstance (seqB, str):
            seqb = [base.base_to_number[x] for x in seqB]
        else:
            seqb = seqB
        
        bp = max (len(seqa), len(seqb) + offset)
        
        s1, s2 = self.generate(bp, None, start_pos, dir, False, True, 0.)

        s1 = s1.get_slice (0, len(seqa))

        if len(seqb) + offset > len(seqa):
            s2 = s2.get_slice (0, len(seqb)) # starts from opposite end
        else:
            s2 = s2.get_slice (bp - offset - len(seqb), len(seqb))
        
        #print >> sys.stderr, len(seqa), len(seqb), len(s1._nucleotides), len(s2._nucleotides)
        s1.set_sequence (seqa)
        s2.set_sequence (seqb)
        
        return s1, s2


class TetramerGenerator (object):
    seqs = ["CTACTATGGCGGGTGATAAAAACGGGAAGAGCATGCCCATCCACGATCG",
            "GGATGGGCATGCTCTTCCCGAACTCAACTGCCTGGTGATACGACGATCG",
            "CGTATCACCAGGCAGTTGAGAACATGCGAGGGTCCAATACCGACGATCG",
            "CGGTATTGGACCCTCGCATGAATTTATCACCCGCCATAGTAGACGATCG",]
    
    def generate(self, cm_pos=np.array([0, 0, 0])):
        seqs = [[base.base_to_number[x] for x in s] for s in TetramerGenerator.seqs]
        
        strands = [base.Strand() for i in range(4)]
        
        half_strand_len = base.BASE_BASE * 21.;
        sqrt1_3 = 1 / np.sqrt(3.)
        sqrt1_2 = 1 / np.sqrt(2.)
        cube_side = half_strand_len  * sqrt1_3;
        rb = cm_pos + np.array([cube_side, -cube_side, -cube_side])
        # this initial direction will give us a 0, -1, +1 direction in the center
        a1 = np.array([-0.438110, -0.815750, 0.377640]) 
        a1 /= np.sqrt(np.dot(a1, a1))
        a3 = np.array([-sqrt1_3, sqrt1_3, sqrt1_3])
        
        a3s = [[-sqrt1_3, sqrt1_3, -sqrt1_3],
               [sqrt1_3, sqrt1_3, sqrt1_3],
               [-sqrt1_3, -sqrt1_3, sqrt1_3],
               [sqrt1_3, -sqrt1_3, -sqrt1_3]]
        
        a1s = [[0, sqrt1_2, sqrt1_2],
               [0, -sqrt1_2, sqrt1_2],
               [0, -sqrt1_2, -sqrt1_2],
               [0, sqrt1_2, -sqrt1_2]]
        
        R = utils.get_rotation_matrix(a3, np.pi/180*(35.9))

        for s in range(4):
            for i in range(49):
                strands[s].add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS*a1, a1, a3, seqs[s][i]))
                if i == 21:
                    a1old = a1
                    a1 = np.array(a1s[s])
                    rback = rb - (base.CM_CENTER_DS - base.POS_BACK) * a1old
                    rback -= (a1 + a1old) * base.BASE_BASE
                    a3 = np.array(a3s[s])
                    rb = rback + (base.CM_CENTER_DS - base.POS_BACK) * a1
                    R = utils.get_rotation_matrix(a3, np.pi/180*(35.9))
                else:
                    a1 = np.dot(R, a1)
                    rb += a3 * base.BASE_BASE
                    
                # the next strand will begin from this base
                if i == 40:
                    na1 = -a1
                    na3 = -a3
                    nrb = np.copy(rb)
                    
            a1 = na1
            a3 = na3
            R = utils.get_rotation_matrix(a3, np.pi/180*(35.9))
            rb = nrb
            
        return strands
