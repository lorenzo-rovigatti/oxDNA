#!/usr/bin/env python3

import sys, os
try:
    import numpy as np
except:
    import mynumpy as np

try:
    box_side = float(sys.argv[1])
    infile = sys.argv[2]
except:
    print("Usage: %s <%s> <%s>" % (sys.argv[0], "box size", "file with sequences"), file=sys.stderr)
    sys.exit(1)
box = np.array ([box_side, box_side, box_side])

try:
    inp = open (infile, 'r')
    inp.close()
except:
    print("Could not open file '%s' for reading. Aborting" % infile, file=sys.stderr)
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
        line = (partition (line.strip (), "//")[0]).strip ()
        macro = (partition (line, "#define ")[2]).strip().split(" ", 1)
        if len(macro) > 1:
            key, val = [x.strip() for x in macro]
            # the f needed by c to interpret numbers as floats must be removed
            # this could be a source of bugs
            val = val.replace("f", "").replace("PI", str(np.pi))
            # this awful exec is needed in order to get the right results out of macro definitions
            exec("tmp = %s" % (val), globals())
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
    '''
    Populates positions, a1s, a3s, with values from their respective mynewpos, mynewa1s, 
    mynewa3s arrays
    Ensures each nt doesn't overlap with any other upon generation
    '''
    overlap = False

    for i in range(len(positions)):
        p = positions[i]
        pa1 = a1s[i]

        for j in range (len (mynewpositions)):
            q = mynewpositions[j]
            qa1 = mynewa1s[j]

            p_pos_back = p + pa1 * POS_BACK
            p_pos_base = p + pa1 * POS_BASE
            q_pos_back = q + qa1 * POS_BACK
            q_pos_base = q + qa1 * POS_BASE

            dr = p_pos_back - q_pos_back
            dr -= box * np.rint (dr / box)
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

        start_pos = np.array(start_pos, dtype=float)

        dir = np.array(dir, dtype=float)
        
        if sequence == None:
            sequence = np.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence += np.random.randint(0, 4, n)
            print("sequence is too short, adding %d random bases" % n, file=sys.stderr)

        # create the sequence of the second strand as made of complementary bases
        sequence2 = [3-s for s in sequence]
        sequence2.reverse()

        # we need to find a vector orthogonal to dir
        dir_norm = np.sqrt(np.dot(dir,dir))
        if dir_norm < 1e-10:
            print("direction must be a valid vector, defaulting to (0, 0, 1)", file=sys.stderr)
            dir = np.array([0, 0, 1])
        else: dir /= dir_norm

        if perp is None or perp is False:
            v1 = np.random.random_sample(3)
            v1 -= dir * (np.dot(dir, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perp
            
        # and we need to generate a rotational matrix
        R0 = get_rotation_matrix(dir, rot)
        #R = get_rotation_matrix(dir, np.deg2rad(35.9))
        R = get_rotation_matrix(dir, [1, "bp"])

        a1 = v1
        a1 = np.dot (R0, a1)
        rb = np.array(start_pos)
        a3 = dir
        for i in range(bp):
            rcdm = rb - CM_CENTER_DS * a1
            mynewpositions.append (rcdm)
            mynewa1s.append(a1)
            mynewa3s.append(a3)
            if i != bp-1:
                a1 = np.dot(R, a1)
                rb += a3 * BASE_BASE

        if double == True:
            a1 = -a1
            a3 = -dir
            R = R.transpose()
            for _ in range(bp):
                rcdm = rb - CM_CENTER_DS * a1
                mynewpositions.append (rcdm)
                mynewa1s.append (a1)
                mynewa3s.append (a3)
                a1 = np.dot(R, a1)
                rb += a3 * BASE_BASE

        assert (len (mynewpositions) > 0)

        return [mynewpositions, mynewa1s, mynewa3s]
    

    def generate_sticky(self, top_bp, bot_bp, bot_seq_start_nn, start_pos=np.array([0, 0, 0]), dir=np.array([0, 0, 1]), perp=False, rot=0.):
        mynewpositions, mynewa1s, mynewa3s = [], [], []

        start_pos = np.array(start_pos, dtype = float)

        dir = np.array(dir, dtype = float)

        # normal vector to the dir respect to top strand
        dir_norm = np.sqrt(np.dot(dir, dir))
        if dir_norm < 1e-10:
            print("STICKY: direction must be a valid vector, defaulting to (0,0,1)", file = sys.stderr)
            dir = np.array([0,0,1])
        else:
            dir /= dir_norm
        
        if perp is None or perp is False:
            v1 = np.random.random_sample(3)
            v1 -= dir * (np.dot(dir, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perp

        # rotational matrix
        R0 = get_rotation_matrix(dir, rot)
        R = get_rotation_matrix(dir, [1, "bp"])

        a1 = v1
        a1 = np.dot(R0, a1)
        rb = np.array(start_pos)
        a3 = dir

        # generates the position of all nt for the top stand seq
        for _ in range(top_bp):
            rcdm = rb - (CM_CENTER_DS * a1)
            mynewpositions.append(rcdm)
            mynewa1s.append(a1)
            mynewa3s.append(a3)
            a1 = np.dot(R, a1)
            rb += a3 * BASE_BASE
        
        # adjusting the position such that we start at the end of the bottom strand seq's
        # sticky end
        for _ in range(int(bot_seq_start_nn) - 1):
            rcdm = rb - (CM_CENTER_DS * a1)
            a1 = np.dot(R, a1)
            rb += a3 * BASE_BASE

        # reversing the direction of generation to meet anti-parallel condition
        a1 = -a1
        a3 = -dir
        R = R.transpose()

        # generates the position of all nt for the bottom strand seq
        for _ in range(bot_bp):
            rcdm = rb - CM_CENTER_DS * a1
            mynewpositions.append (rcdm)
            mynewa1s.append (a1)
            mynewa3s.append (a3)
            a1 = np.dot(R, a1)
            rb += a3 * (BASE_BASE)

        assert (len (mynewpositions) > 0)

        return [mynewpositions, mynewa1s, mynewa3s]


# helper function to avoid code repeition in the read_strands() function
def topo_file_helper(nt_count: int, strand_count: int, output_file, strand):
    print(strand_count, strand[0], -1, nt_count + 1, file=output_file)
    nt_count += 1
    for i in range (1, len(strand) - 1):
        print(strand_count, strand[i], nt_count - 1, nt_count + 1, file=output_file)
        nt_count += 1
    print(strand_count, strand[-1], nt_count - 1, -1, file=output_file)
    nt_count += 1
    strand_count += 1
    return nt_count, strand_count


# helper function to avoid code repeition in the read_strands() function
def initial_strand_gen_randomizer():
    cdm = np.random.random_sample(3) * box
    axis = np.random.random_sample(3)
    axis /= np.sqrt(np.dot(axis, axis))
    return cdm, axis


def read_strands(filename):
    """
    The main() function for this script
    Reads a text file with the following format:
    - Each line contains the sequence for a single strand (A,C,T,G)
    - Lines begining in DOUBLE produce double stranded DNA
    - Lines beginning in STICKY produce a dsDNA with sticky ends, sticky ends have the same length
        - We define sticky end as the end strand of a linear dsDNA without its complementary strand

    Ex: Two ssDNA (single stranded DNA)
    ATATATA
    GCGCGCG

    Ex: Two strands, one double stranded, the other single stranded.
    DOUBLE AGGGCT
    CCTGTA

    Ex: one sticky dsDNA
    STICKY 4 TTTTGAGACACG AAAA

    The number after STICKY denotes the length of the sticky end
    
    In the above example, the sticky ends are TTTT and AAAA (they must match in length), they
    do not have to be the same sequence
    
    The code will automatically generate the complementary sequences making up the second strand
    after the "sticky end" of the first strand
    
    The first strand in the above example is TTTTGAGACACG

    """

    try:
        infile = open (filename)
    except:
        print("Could not open file ", filename, "Aborting now", file=sys.stderr)
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
            print("## Found duplex of %i bases" % (length), file=sys.stderr)
            nnucl += 2 * length
            nstrands += 2
        elif line[:6] == 'STICKY':
            parsed_line = line.split()

            top_strand_seq = parsed_line[2]
            bot_strand_seq = parsed_line[3]
            top_strand_len = len(top_strand_seq)
            bot_strand_len = len(bot_strand_seq)

            len_complement = 0
            if int(parsed_line[1]) != 0:
                len_complement = top_strand_len - (int(parsed_line[1]))

            total_len = top_strand_len + bot_strand_len + len_complement

            print("## Found sticky duplex of %i total bases" % (total_len), file = sys.stderr)
            nnucl += total_len
            nstrands += 2
        else:
            line = line.split()[0]
            print("## Found single strand of %i bases" % (len(line)), file=sys.stderr)
            nnucl += len(line)
            nstrands += 1

    infile.seek(0)

    print("## nstrands, nnucl = ", nstrands, nnucl, file=sys.stderr)

    # here we generate the topology file
    try:
        out = open ("generated.top", "w")
    except:
        print("Could not open generated.top for writing. Aborting", file=sys.stderr)
        sys.exit(4)

    print(nnucl, nstrands, file=out)
    myns, mynn = 1, 0
    lines = infile.readlines()
    for line in lines:
        line = line.upper().strip()
        if len(line) == 0:
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1].upper()
            mynn, myns = topo_file_helper(mynn, myns, out, line)

            # get the compl sequence in numbers
            seq = [3 - base_to_number[x] for x in line]
            # put it back in letters
            line = [number_to_base[x] for x in seq[::-1]]

            mynn, myns = topo_file_helper(mynn, myns, out, line)

        elif line[:6] == 'STICKY':
            splited_line = line.split()
            
            bot_start_nn = int(splited_line[1])
            assert bot_start_nn > 0, "bot strand nt start # must > 0, maybe use double?"
            
            top_strand_seq = splited_line[2].upper()
            mynn, myns = topo_file_helper(mynn, myns, out, top_strand_seq)

            # bottom seq
            bot_strand_hang = list(splited_line[3].upper())
            assert len(bot_strand_hang) == bot_start_nn, "the length of the stickey ends must be equal"

            truncated_top_seq = splited_line[2][(bot_start_nn):]
            bot_strand_seq_in_numbers = [3 - base_to_number[x] for x in truncated_top_seq]
            bot_strand_seq = bot_strand_hang + [number_to_base[x] for x in bot_strand_seq_in_numbers[::-1]]

            mynn, myns = topo_file_helper(mynn, myns, out, bot_strand_seq)

        else:
            line = line.split()[0].upper()
            mynn, myns = topo_file_helper(mynn, myns, out, line)
    out.close ()
    infile.seek (0)

    # generate the strands
    double = StrandGenerator()
    lines = infile.readlines()
    nlines = len(lines)
    i = 1
    for line in lines:
        line = line.upper().strip()
        # skip empty lines
        if len(line) == 0: 
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1]
            seq = [base_to_number[x] for x in line]
            length = len(line)
            print("## Adding duplex of %i bases" % (length), file=sys.stderr)

            cdm, axis = initial_strand_gen_randomizer()
            newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=True)
            
            while not add_strands(newpositions, newa1s, newa3s):
                cdm, axis = initial_strand_gen_randomizer()
                newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=True)
                print("##  trying %i" % (i), file=sys.stderr)

            print("##  done line %i / %i, now at %i/%i" % (i, nlines, len(positions), nnucl), file=sys.stderr)
        elif line[:6] == 'STICKY':
            splited_line = line.split()

            bot_start_nn = splited_line[1]
            top_strand_seq = splited_line[2]
            bot_strand_seq = splited_line[3]

            top_bp_len = len(top_strand_seq)
            
            len_complement = 0
            if int(bot_start_nn) != 0:
                len_complement = top_strand_len - (int(bot_start_nn))

            bot_bp_len = len(bot_strand_seq) + len_complement
            total_len = top_strand_len + bot_bp_len

            print("## Adding sticky duplex of %i bases" % (total_len), file = sys.stderr)

            cdm, axis = initial_strand_gen_randomizer()
            newpositions, newa1s, newa3s = double.generate_sticky(top_bp_len, bot_bp_len, bot_start_nn, start_pos=cdm, dir=axis)
            
            while not add_strands(newpositions, newa1s, newa3s):
                cdm, axis = initial_strand_gen_randomizer()
                newpositions, newa1s, newa3s = double.generate_sticky(top_bp_len, bot_bp_len, bot_start_nn, start_pos=cdm, dir=axis)
                print("## tryin %i" % i, file = sys.stderr)

            print("## done line %i / %i, now at %i / %i" % (i, nlines, len(positions), nnucl), file = sys.stderr)
        else:
            seq = [base_to_number[x] for x in line]
            cdm, axis = initial_strand_gen_randomizer()
            print("## Adding single strand of %i bases" % (len(line)), file=sys.stderr)
            newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=False)
            
            while not add_strands(newpositions, newa1s, newa3s):
                cdm, axis = initial_strand_gen_randomizer()
                newpositions, newa1s, newa3s = double.generate(len(line), sequence=seq, dir=axis, start_pos=cdm, double=False)
                print("##  trying  %i" % (i), file=sys.stderr)

            print("##  done line %i / %i, now at %i/%i" % (i, nlines, len(positions), nnucl), file=sys.stderr)

        #print "AAA", i, len(positions)

        i += 1

    if not len(positions) == nnucl:
        print(len(positions), nnucl)
        raise AssertionError

    # here we generate the configuration file (coordinates)
    try:
        outfile = open ('generated.dat', 'w')
    except:
        print("Could not open generated.dat for writing.  Aborting", file=sys.stderr)
        sys.exit(5)

    # print("---")
    # print(positions)

    print("t = 0", file=outfile)
    print("b = ", box_side, box_side, box_side, file=outfile)
    print("E = 0. 0. 0.", file=outfile)
    for i in range(nnucl):
        print(positions[i][0], positions[i][1], positions[i][2], end=' ', file=outfile)
        print(a1s[i][0], a1s[i][1], a1s[i][2], end=' ', file=outfile)
        print(a3s[i][0], a3s[i][1], a3s[i][2], end=' ', file=outfile)
        print(0., 0., 0., 0., 0., 0., file=outfile) # v and L

    outfile.close()
    print("## ALL DONE. just generated 'generated.dat' and 'generated.top'", file=sys.stderr)


read_strands (infile)
