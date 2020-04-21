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
    def generate(self, bp, sequence=None, start_pos=np.array([0, 0, 0]), dir=np.array([0, 0, 1]), perp=None, rot=0., double=True, circular=False, DELTA_LK=0, BP_PER_TURN=10.34, ds_start=None, ds_end=None, force_helicity=False):
        """
        Generate a strand of DNA.
            - linear, circular (circular)
            - ssDNA, dsDNA (double)
            - Combination of ss/dsDNA (ds_start, ds_end)
            Note: Relevent argument(s) in parentheses.

        Arguments:
        bp --- Integer number of bp/nt (required)
        sequence --- Array of integers or string. Should be same length as bp (default None)
            Default (None) generates a random sequence.
            Ex: [0,1,2,3,0]
            Ex: "AGCTA"
            See dictionary base.base_to_number for int/char conversion {0:'A'}
        start_pos --- Location to begin building the strand (default np.array([0, 0, 0]))
        dir --- a3 vector, orientation of the base (default np.array([0, 0, 1]))
        perp --- Sets a1 vector, the orientation of the backbone. (default False)
            Must be perpendicular to dir (as a1 must be perpendicular to a3)
            If perp is None or False, perp is set to a random orthogonal angle
        rot --- Rotation of first bp (default 0.)
        double --- Generate dsDNA (default True)
        circular --- Generate closed circular DNA (defalt False)
            Limitations...
            For ssDNA (double=False): bp >= 4
            For dsDNA (double=True) : bp >= 30
            Will throw warnings. Allowed, but use at your own risk.
        DELTA_LK --- Integer change in linking number from Lk0 (default 0)
            Only valid if circular==True
        BP_PER_TURN --- Base pairs per complete 2*pi helix turn. (default 10.34)
            Only valid if circular==True
        ds_start --- Index (from 0) to begin double stranded region (default None)
        ds_end --- Index (from 0) to end double stranded region (default None)
            Default is None, which is entirely dsDNA; sets ds_start = 0, ds_end=bp
            Ex: ds_start=0, ds_end=10 will create a double stranded region on bases
                range(0,10): [0,1,2,3,4,5,6,7,8,9]
            Note: To generate a nicked circular dsDNA, manually change state with
                  {Strand}.make_noncircular()
        force_helicity --- Force generation of helical strands. Use helicity by default
            for bp > 30. Warns from 18 to 29. Will crash oxDNA below 18. (default False)

        Note: Minimuim circular duplex is 18. Shorter circular strands disobey FENE.
        For shorter strands, circular ssDNA is generated in a circle instead of having
        imposed helicity.

        Examples:
        Generate ssDNA:
            generate(bp=4,sequence=[0,1,2,3],double=False,circular=False)
        Generate circular dsDNA with +2 Linking number:
            generate(bp=45,double=True,circular=True,DELTA_LK=2)
        Generate a circular ssDNA (45nt) with ssDNA (25nt) annealed to indices 0 to 24:
            generate(bp=45,double=True,circular=True,ds_start=0,ds_end=25)
        """
        # we need numpy array for these
        start_pos = np.array(start_pos, dtype=float)
        dir = np.array(dir, dtype=float)
        if isinstance(sequence, list):
            sequence = np.array(sequence)

        # Loads of input checking...
        if isinstance(sequence, str):
            try:
                sequence = [base.base_to_number[x] for x in sequence]
            except KeyError:
                base.Logger.die("Key Error: sequence is invalid" % n)
        if sequence is None:
            sequence = np.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence = np.append(sequence, np.random.randint(0, 4, n))
            base.Logger.log("sequence is too short, adding %d random bases" % n, base.Logger.WARNING)

        if circular == True and bp < 30:
            # 30 is about the cut off for circular dsDNA. Anything shorter will probably clash.
            # oxDNA can relax down to 18.
            # 4 is about the cut off for circular ssDNA. Use dsDNA cutoff for saftey.
            base.Logger.log("sequence is too short! Proceed at your own risk", base.Logger.WARNING)

        option_use_helicity = True
        if circular == True and bp < 30 and double == False:
            base.Logger.log("sequence is too short! Generating ssDNA without imposed helicity", base.Logger.WARNING)
            # Do not impose helcity to generate shorter circular ssDNA
            if not force_helicity:
                option_use_helicity = False

        if ds_start == None:
            ds_start = 0
        if ds_end == None:
            ds_end = bp
        if ds_start > ds_end:
            base.Logger.die("ds_end > ds_start")
        if  ds_end > bp:
            base.Logger.die("ds_end > bp")

        # we need to find a vector orthogonal to dir
        dir_norm = np.sqrt(np.dot(dir,dir))
        if dir_norm < 1e-10:
            base.Logger.log("direction must be a valid vector, defaulting to (0, 0, 1)", base.Logger.WARNING)
            dir = np.array([0, 0, 1])
        else:
            dir /= dir_norm

        if perp is None or perp is False:
            v1 = np.random.random_sample(3)
            v1 -= dir * (np.dot(dir, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perp;

        # Setup initial parameters
        ns1 = base.Strand()
        # and we need to generate a rotational matrix
        R0 = utils.get_rotation_matrix(dir, rot)
        #R = get_rotation_matrix(dir, np.deg2rad(35.9))
        R = utils.get_rotation_matrix(dir, [1, BP])
        a1 = v1
        a1 = np.dot (R0, a1)
        rb = np.array(start_pos)
        a3 = dir

        # Circular strands require a continuious deformation of the ideal helical pitch
        if circular == True:
            # Unit vector orthogonal to plane of torus
            # Note: Plane of torus defined by v1,dir
            torus_perp = np.cross(v1,dir)
            # bp-bp rotation factor to yield a smooth deformation along the torus
            smooth_factor = np.mod(bp,BP_PER_TURN)/float(bp) + 1
            # Angle between base pairs along torus
            angle = 2. * np.pi / float(bp)
            # Radius of torus
            radius = base.FENE_R0_OXDNA / math.sqrt(2. * (1. - math.cos(angle)));

        if circular == True and option_use_helicity:
            # Draw backbone in a helical spiral around a torus
            # Draw bases pointing to center of torus
            for i in range(bp):
                # Torus plane defined by dir and v1
                v_torus = v1 *base.BASE_BASE * math.cos(i * angle) + \
                        dir * base.BASE_BASE * math.sin(i * angle)
                rb += v_torus

                # a3 is tangent to the torus
                a3 = v_torus/np.linalg.norm(v_torus)
                R = utils.get_rotation_matrix(a3, [i * (round(bp/BP_PER_TURN) + DELTA_LK)/float(bp) * 360, DEGREES])

                # a1 is orthogonal to a3 and the torus normal
                a1 = np.cross (a3, torus_perp)

                # Apply the rotation matrix
                a1 = np.dot(R, a1)
                ns1.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence[i]))
            ns1.make_circular(check_join_len=True)
        elif circular == True and not option_use_helicity:
            for i in xrange(bp):
                rbx = math.cos (i * angle) * radius + 0.34 * math.cos(i * angle)
                rby = math.sin (i * angle) * radius + 0.34 * math.sin(i * angle)
                rbz = 0.
                rb = np.array([rbx, rby, rbz])
                a1x = math.cos (i * angle)
                a1y = math.sin (i * angle)
                a1z = 0.
                a1 = np.array([a1x, a1y, a1z])
                ns1.add_nucleotide(base.Nucleotide(rb, a1, np.array([0, 0, 1]), sequence[i]))
            ns1.make_circular(check_join_len=True)
        else:
            # Add nt in canonical double helix
            for i in range(bp):
                ns1.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence[i]))
                if i != bp-1:
                    a1 = np.dot(R, a1)
                    rb += a3 * base.BASE_BASE

        # Fill in complement strand
        if double == True:
            ns2 = base.Strand()
            for i in reversed(range(ds_start, ds_end)):
                # Note that the complement strand is built in reverse order
                nt = ns1._nucleotides[i]
                a1 = -nt._a1
                a3 = -nt._a3
                nt2_cm_pos = -(base.FENE_EPS + 2*base.POS_BACK) * a1 + nt.cm_pos
                ns2.add_nucleotide(base.Nucleotide(nt2_cm_pos, a1, a3, 3-sequence[i]))
            if ds_start == 0 and ds_end == bp and circular == True:
                ns2.make_circular(check_join_len=True)
            return ns1, ns2
        else:
            return ns1

    def generate_or_sq(self, bp, sequence=None, start_pos=np.array([0,0,0]), dir=np.array([0, 0, 1]), perp=None, double=True, rot=0., angle=np.pi/180*33.75, length_change=0, region_begin=0, region_end=0):
        if length_change and len(region_begin) != len(region_end):
            if (len(region_end) + 1) == len(region_begin):
                base.Logger.log("the lengths of begin (%d) and end (%d) arrays are mismatched; I will try to proceed by using the number of basepairs as the last element of the end array" % (len(region_begin), len(region_end)), base.Logger.WARNING)
                region_end.append(bp + 1)
            else:
                base.Logger.die("the lengths of begin (%d) and end (%d) arrays are unrecoverably mismatched" % (len(region_begin), len(region_end)))
        
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
    
    def generate_rw (self, sequence, start_pos=np.array([0.,0.,0.])):
        """
        Generate ssDNA as a random walk (high-energy configurations are possible):
            generate(bp=45,double=False,circular=False,random_walk=True)
        """
        # random walk generator
        base.Logger.log("Generating strand as a random walk. Remember to equilibrate the configuration with MC", base.Logger.WARNING)
        d = np.array ([0.7525, 0., 0.])
        pos = start_pos
        rw = []
        rw.append(pos)
        for i, s in enumerate(sequence[1:]):
            overlap = True
            while overlap:
                overlap = False
                R = utils.get_random_rotation_matrix()
                dd = np.dot (R, d)
                trypos = pos + np.dot (R, d);
                overlap = False
                for r in rw:
                    dd = trypos - r
                    if np.dot (dd, dd) < 0.40 * 0.40:
                        overlap = True
            pos = trypos
            rw.append (pos)
            #print >> sys.stderr, "# nucleotide", i + 1, "inserted"
        
        # we get the a1 vectors in a smart way
        a1s = []
        d = rw[1] - rw[0]
        a1s.append (d / np.sqrt (np.dot(d, d)))
        
        for i in range (1, len(rw) - 1):
            d = (rw[i + 1] + rw[i - 1]) * 0.5
            d = rw[i] - d
            a1s.append (d / np.sqrt(np.dot (d, d)))
        
        d = rw[len(rw)-1] - rw[len(rw)-2]
        a1s.append (d / np.sqrt (np.dot(d, d)))
        
        s = base.Strand()
        for i, r in enumerate(rw):
            a1, a2, a3 = utils.get_orthonormalized_base (a1s[i], utils.get_random_vector(), utils.get_random_vector()) 
            #print np.dot(a1, a1), np.dot(a2, a2), np.dot(a3, a3), np.dot(a1, a2), np.dot(a1, a3), np.dot(a2, a3)
            ## POS_BACK is negative!
            cm = r + a1s[i] * abs(base.POS_BACK)
            s.add_nucleotide (base.Nucleotide (cm, a1, a3, sequence[i]))
        
        return s
        
        # we have to initialize the a3s
       


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
