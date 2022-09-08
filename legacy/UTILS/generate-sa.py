#!/usr/bin/env python
"""
Standalone version of generate.py
Created for stable release end-users and machines without numpy support.
This file is only supported for stable release, and does not include more recent functionality.

(05 December 2012)
"""

import sys, os
try:
    import numpy as np
except:
    import mynumpy as np

try:
    box_side = float(sys.argv[1])
    infile = sys.argv[2]
except:
    print >> sys.stderr, "Usage: %s <%s> <%s>" % (sys.argv[0], "box size", "file with sequences")
    sys.exit(1)
box = np.array ([box_side, box_side, box_side])

try:
    inp = open (infile, 'r')
    inp.close()
except:
    print >> sys.stderr, "Could not open file '%s' for reading. Aborting" % infile
    sys.exit(2)

# return parts of a string
def partition(s, d):
    if d in s:
        sp = s.split(d, 1)
        return sp[0], d, sp[1]
    else:
        return s, "", ""

# every defined macro in model.h must be imported in this module
def import_model_constants():
    PI = np.pi
    model = os.path.join(os.path.dirname(__file__), "../src/model.h")
    f = open(model)
    for line in f.readlines():
        # line = line.strip().partition("//")[0].strip()
        line = (partition (line.strip (), "//")[0]).strip ()
        #macro = line.partition("#define ")[2].strip().split(" ", 1)
        macro = (partition (line, "#define ")[2]).strip().split(" ", 1)
        if len(macro) > 1:
            key, val = [x.strip() for x in macro]
            # the f needed by c to interpret numbers as floats must be removed
            # this could be a source of bugs
            val = val.replace("f", "")
            # this awful exec is needed in order to get the right results out of macro definitions
            exec "tmp = %s" % (val)
            globals()[key] = tmp
    f.close()

import_model_constants()

CM_CENTER_DS = POS_BASE + 0.2
BASE_BASE = 0.3897628551303122

RC2_BACK = EXCL_RC1**2
RC2_BASE = EXCL_RC2**2
RC2_BACK_BASE = EXCL_RC3**2

number_to_base = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T'}
base_to_number = {'A' : 0, 'a' : 0, 'G' : 1, 'g' : 1,
                  'C' : 2, 'c' : 2, 'T' : 3, 't' : 3}

positions = []
a1s = []
a3s = []
newpositions = []
newa1s = []
newa3s = []

def add_strands (mynewpositions, mynewa1s, mynewa3s):
    overlap = False

    #print len (mynewpositions), "@@@"

    for i in xrange(len(positions)):
        p = positions[i]
        pa1 = a1s[i]

        for j in xrange (len (mynewpositions)):
            q = mynewpositions[j]
            qa1 = mynewa1s[j]

            p_pos_back = p + pa1 * POS_BACK
            p_pos_base = p + pa1 * POS_BASE
            q_pos_back = q + qa1 * POS_BACK
            q_pos_base = q + qa1 * POS_BASE

            dr = p_pos_back - q_pos_back
            dr -= box * np.rint (dr / box)
            #print RC2_BACK, np.dot (dr, dr)
            if np.dot(dr, dr) < RC2_BACK:
                overlap = True

            dr = p_pos_base -  q_pos_base
            dr -= box * np.rint (dr / box)
            if np.dot(dr, dr) < RC2_BASE:
                overlap = True

            dr = p_pos_back - q_pos_base
            dr -= box * np.rint (dr / box)
            if np.dot(dr, dr) < RC2_BACK_BASE:
                overlap = True

            dr = p_pos_base - q_pos_back
            dr -= box * np.rint (dr / box)
            if np.dot(dr, dr) < RC2_BACK_BASE:
                overlap = True

            if overlap:
                return False

    if not overlap:
        for p in mynewpositions:
            positions.append(p)
        for p in mynewa1s:
            a1s.append (p)
        for p in mynewa3s:
            a3s.append (p)

    return True


def get_rotation_matrix(axis, anglest):
    # the argument anglest can be either an angle in radiants
    # (accepted types are float, int or np.float64 or np.float64)
    # or a tuple [angle, units] where angle a number and
    # units is a string. It tells the routine whether to use degrees,
    # radiants (the default) or base pairs turns
    if not isinstance (anglest, (np.float64, np.float32, float, int)):
        if len(anglest) > 1:
            if anglest[1] in ["degrees", "deg", "o"]:
                #angle = np.deg2rad (anglest[0])
                angle = (np.pi / 180.) * (anglest[0])
            elif anglest[1] in ["bp"]:
                angle = int(anglest[0]) * (np.pi / 180.) * (35.9)
            else:
                angle = float(anglest[0])
        else:
            angle = float(anglest[0])
    else:
        angle = float(anglest) # in degrees, I think

    axis = np.array(axis)
    axis /= np.sqrt(np.dot(axis, axis))

    ct = np.cos(angle)
    st = np.sin(angle)
    olc = 1. - ct
    x, y, z = axis

    return np.array([[olc*x*x+ct, olc*x*y-st*z, olc*x*z+st*y],
                    [olc*x*y+st*z, olc*y*y+ct, olc*y*z-st*x],
                    [olc*x*z-st*y, olc*y*z+st*x, olc*z*z+ct]])

class StrandGenerator (object):
    def generate(self, bp, sequence=None, start_pos=np.array([0, 0, 0]), dir=np.array([0, 0, 1]), perp=False, double=True, rot=0.):
        mynewpositions, mynewa1s, mynewa3s = [], [], []
        #assert (len (newpositions) == 0)
        # we need a numpy array for these
        start_pos = np.array(start_pos, dtype=float)
        dir = np.array(dir, dtype=float)
        if sequence == None:
            sequence = np.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence += np.random.randint(0, 4, n)
            print >> sys.stderr, "sequence is too short, adding %d random bases" % n

        # create the sequence of the second strand as made of complementary bases
        sequence2 = [3-s for s in sequence]
        sequence2.reverse()

        # we need to find a vector orthogonal to dir
        dir_norm = np.sqrt(np.dot(dir,dir))
        if dir_norm < 1e-10:
            print >> sys.stderr, "direction must be a valid vector, defaulting to (0, 0, 1)"
            dir = np.array([0, 0, 1])
        else: dir /= dir_norm

        if perp is None or perp is False:
            v1 = np.random.random_sample(3)
            v1 -= dir * (np.dot(dir, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perp;

        # and we need to generate a rotational matrix
        R0 = get_rotation_matrix(dir, rot)
        #R = get_rotation_matrix(dir, np.deg2rad(35.9))
        R = get_rotation_matrix(dir, [1, "bp"])

        #ns1 = base.Strand()

        a1 = v1
        a1 = np.dot (R0, a1)
        rb = np.array(start_pos)
        a3 = dir
        for i in range(bp):
            #ns1.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence[i]))
            rcdm = rb - CM_CENTER_DS * a1
            mynewpositions.append (rcdm)
            mynewa1s.append(a1)
            mynewa3s.append(a3)
            #print >> outfile, rcdm[0], rcdm[1], rcdm[2],
            #print >> outfile, a1[0], a1[1], a1[2],
            #print >> outfile, a3[0], a3[1], a3[2],
            #print >> outfile, 0., 0., 0., 0., 0., 0. # v and L
            if i != bp-1:
                a1 = np.dot(R, a1)
                rb += a3 * BASE_BASE

        if double == True:
            a1 = -a1
            a3 = -dir
            R = R.transpose()
            #ns2 = base.Strand()
            for i in range(bp):
                #ns2.add_nucleotide(base.Nucleotide(rb - base.CM_CENTER_DS * a1, a1, a3, sequence2[i]))
                rcdm = rb - CM_CENTER_DS * a1
                mynewpositions.append (rcdm)
                mynewa1s.append (a1)
                mynewa3s.append (a3)
                #print >> outfile, rcdm[0], rcdm[1], rcdm[2],
                #print >> outfile, a1[0], a1[1], a1[2],
                #print >> outfile, a3[0], a3[1], a3[2],
                #print >> outfile, 0., 0., 0., 0., 0., 0. # v and L
                a1 = np.dot(R, a1)
                rb += a3 * BASE_BASE

        #print "eccoce", len(mynewpositions)
        assert (len (mynewpositions) > 0)

        return [mynewpositions, mynewa1s, mynewa3s]


def read_strands(filename):
    """
    The main() function for this script
    Reads a text file with the following format:
    - Each line contains the sequence for a single strand (A,C,T,G)
    - Lines begining in DOUBLE produce double stranded DNA

    Ex: Two ssDNA (single stranded DNA)
    ATATATA
    GCGCGCG

    Ex: Two strands, one double stranded, the other single stranded.
    DOUBLE AGGGCT
    CCTGTA

    """
    # prendiamo le sequenze da un file; ogni riga uno strand; se e' un
    # double strand, ci deve essere DOUBLE davanti

    try:
        infile = open (filename)
    except:
        print >> sys.stderr, "Could not open file ", filename, "Aborting now"
        sys.exit(2)

    if len(sys.argv) > 1:
        side = float(sys.argv[1])
    else: side = 50

    # get the number of strands and nucleotides
    nstrands, nnucl = 0, 0

    lines = infile.readlines()
    for line in lines:
        line = line.upper().strip()
        if len(line) == 0:
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1]
            length = len(line)
            print >> sys.stderr, "## Found duplex of %i bases" % (length)
            nnucl += 2 * length
            nstrands += 2
        else:
            line = line.split()[0]
            print >> sys.stderr, "## Found single strand of %i bases" % (len(line))
            nnucl += len(line)
            nstrands += 1

    infile.seek(0)

    print >> sys.stderr, "## nstrands, nnucl = ", nstrands, nnucl

    # here we generate the topology file
    try:
        out = open ("generated.top", "w")
    except:
        print >> sys.stderr, "Could not open generated.top for writing. Aborting"
        sys.exit(4)

    print >> out, nnucl, nstrands
    myns, mynn = 1, 0
    lines = infile.readlines()
    for line in lines:
        line = line.upper().strip()
        if len(line) == 0:
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1].upper()
            print >> out, myns, line[0], -1, mynn + 1
            mynn += 1
            for i in xrange (1, len(line) - 1):
                print >> out, myns, line[i], mynn - 1, mynn + 1
                mynn += 1
            print >> out, myns, line[-1], mynn - 1, -1
            mynn += 1
            myns += 1

            # get the compl sequence in numbers
            seq = [3 - base_to_number[x] for x in line]
            # put it back in letters
            line = [number_to_base[x] for x in seq[::-1]]

            print >> out, myns, line[0], -1, mynn + 1
            mynn += 1
            for i in xrange (1, len(line) - 1):
                print >> out, myns, line[i], mynn - 1, mynn + 1
                mynn += 1
            print >> out, myns, line[-1], mynn - 1, -1
            mynn += 1
            myns += 1
        else:
            line = line.split()[0].upper()
            print >> out, myns, line[0], -1, mynn + 1
            mynn += 1
            for i in xrange (1, len(line) - 1):
                print >> out, myns, line[i], mynn - 1, mynn + 1
                mynn += 1
            print >> out, myns, line[-1], mynn -1, -1
            mynn += 1
            myns += 1
    out.close ()
    infile.seek (0)

    # generate the strands
    double = StrandGenerator()
    lines = infile.readlines()
    nlines = len(lines)
    i = 1
    for line in lines:
        #print "AAA", i, len(positions)
        line = line.upper().strip()
        # skip empty lines
        if len(line) == 0: continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1]
            seq = [base_to_number[x] for x in line]
            length = len(line)
            print >> sys.stderr, "## Adding duplex of %i bases" % (length)
            cdm = np.random.random_sample(3) * box
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))
            newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=True)
            while not add_strands(newpositions, newa1s, newa3s):
                cdm = np.random.random_sample(3) * box
                axis = np.random.random_sample(3)
                axis /= np.sqrt(np.dot(axis, axis))
                newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=True)
                print >> sys.stderr, "##  trying %i" % (i)
            #print >> sys.stderr, "##  done line %i / %i" % (i, nlines)
            print >> sys.stderr, "##  done line %i / %i, now at %i/%i" % (i, nlines, len(positions), nnucl)
        else:
            seq = [base_to_number[x] for x in line]
            cdm = np.random.random_sample(3) * box
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))
            print >> sys.stderr, "## Adding single strand of %i bases" % (len(line))
            newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=False)
            while not add_strands(newpositions, newa1s, newa3s):
                cdm = np.random.random_sample(3) * box
                axis = np.random.random_sample(3)
                axis /= np.sqrt(np.dot(axis, axis))
                newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=False)
                print >> sys.stderr, "##  trying  %i" % (i)
            print >> sys.stderr, "##  done line %i / %i, now at %i/%i" % (i, nlines, len(positions), nnucl)

        #print "AAA", i, len(positions)

        i += 1

    if not len(positions) == nnucl:
        print len(positions), nnucl
        raise AssertionError

    # here we generate the configuration file (coordinates)
    try:
        outfile = open ('generated.dat', 'w')
    except:
        print >> sys.stderr, "Could not open generated.dat for writing.  Aborting"
        sys.exit(5)

    print >> outfile, "t = 0"
    print >> outfile, "b = ", box_side, box_side, box_side
    print >> outfile, "E = 0. 0. 0."
    for i in xrange(nnucl):
        print >> outfile, positions[i][0], positions[i][1], positions[i][2],
        print >> outfile, a1s[i][0], a1s[i][1], a1s[i][2],
        print >> outfile, a3s[i][0], a3s[i][1], a3s[i][2],
        print >> outfile, 0., 0., 0., 0., 0., 0. # v and L

    outfile.close()
    print >> sys.stderr, "## ALL DONE. just generated 'generated.dat' and 'generated.top'"


read_strands (infile)

