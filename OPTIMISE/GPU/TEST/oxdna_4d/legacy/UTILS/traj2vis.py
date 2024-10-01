#!/usr/bin/env python

import sys, os
try:
    import numpy as np
except:
    import mynumpy as np

GROOVE_ENV_VAR = 'OXDNA_GROOVE'

try:
    which = sys.argv[1]
    conffile = sys.argv[2]
    topfile = sys.argv[3]
except:
    print >> sys.stderr, "Usage: %s <%s> <%s> <%s>" % (sys.argv[0], "pdb|xyz", "trajectory", "topology")
    sys.exit(1)

if which not in ['pdb','xyz']:
    print >> sys.stderr, "choose either pdb or xyz in the second argument. Aborting"
    sys.exit (3)

pdb, xyz = False, False
if which == 'pdb':
    pdb = True
elif which == 'xyz':
    xyz = True
else:
    pass

box = np.array ([1.,1.,1.])

try:
    inp = open (conffile, 'r')
    inp.close()
except:
    print >> sys.stderr, "Could not open file '%s' for reading. Aborting" % conffile
    sys.exit(2)

try:
    inp = open (topfile, 'r')
    inp.close()
except:
    print >> sys.stderr, "Could not open file '%s' for reading. Aborting" % topfile
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

def get_rotation_matrix(axis, anglest):
    # the argument anglest can be either an angle in radiants
    # (accepted types are float, int or np.float64 or np.float64)
    # or a tuple [angle, units] where angle a number and 
    # units is a string. It tells the routine whether to use degrees,
    # radiants (the default) or base pairs turns
    if not isinstance (anglest, (np.float64, np.float32, float, int)):
        if len(anglest) > 1:
            if anglest[1] in ["degrees", "deg", "o"]:
                angle = np.deg2rad (anglest[0])
            elif anglest[1] in ["bp"]:
                angle = int(anglest[0]) * np.deg2rad(35.9)
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


# read topology
inp = open (topfile, 'r')
lines = inp.readlines()
inp.close()

nnucl, nstrands = [int(x) for x in lines[0].split()]

print >> sys.stderr, "## Topology: found %i strands, %i nucl" % (nstrands, nnucl)

strandid = [None for x in xrange(nnucl)]
basetype = [None for x in xrange(nnucl)]
nn3 = [None for x in xrange(nnucl)]
nn5 = [None for x in xrange(nnucl)]

i = 0
for line in lines[1:]:
    words = line.split()
    strandid[i] = int(words[0]) -1 
    basetype[i] = words[1]
    nn3[i] = int (words[2])
    nn5[i] = int (words[3])
    i += 1

# read conf
inp = open (conffile, 'r')

if xyz:
    outfile = conffile + '.xyz'
elif pdb:
    outfile = conffile + '.pdb'
out = open (outfile, 'w')

nconfs = 0
times = []
box = []

# for chimera output
strtypes = ["ALA","GLY","CYS","ARG","PHE","LYS","SER","PRO","VAL","ASN","ASP"]

while True:
    rcs = [None for x in xrange(nnucl)]
    a1s = [None for x in xrange(nnucl)]
    a2s = [None for x in xrange(nnucl)]
    a3s = [None for x in xrange(nnucl)]
    
    try:
        times.append (int(inp.readline().split()[2]))
        box = np.array([float(x) for x in inp.readline().split()[2:]])
        print "## conf %i, time = %i box =(%lf %lf %lf)" % (nconfs, times[nconfs], box[0], box[1], box[2])
    except:
        print >> sys.stderr, "## End of trajectory...", nconfs
        break
    
    inp.readline () # to remove the energy line 
    for i in xrange (nnucl):
        line = inp.readline()
        words = line.split()
        rcs[i] = np.array([float(x) for x in words[0:3]])
        a1s[i] = np.array([float(x) for x in words[3:6]])
        a3s[i] = np.array([float(x) for x in words[6:9]])
        a2s[i] = np.cross (a3s[i], a1s[i])

    # strand-wise periodic boundary conditions
    cdms = [np.array([0.,0.,0.]) for x in xrange(nstrands)]
    nin = [0 for x in xrange (nstrands)]
    for i in xrange (nnucl):
        cdms[strandid[i]] += rcs[i]
        nin[strandid[i]] += 1

    ocdm = np.array([0.,0.,0.])
    for i in xrange (nstrands):
        cdms[i] /= nin[i]
        ocdm += cdms[i] - box * np.rint (cdms[i] / box)
    ocdm /= float(nstrands)

    if xyz:
        print >> out, 2 * nnucl
        print >> out
        for i in xrange(nnucl):
            sid = strandid[i]
            #rnow = rcs[i] - cdms[strandid[i]]  + POS_BACK * a1s[i]
            if os.environ.get(GROOVE_ENV_VAR) == '1':
                rnow = rcs[i] - cdms[sid] + POS_MM_BACK1 * a1s[i] + POS_MM_BACK2 * a2s[i] + (cdms[sid] - ocdm) - box * np.rint ((cdms[sid] - ocdm) / box)
            else:
                rnow = rcs[i] - cdms[sid] + POS_BACK * a1s[i] + (cdms[sid] - ocdm) - box * np.rint ((cdms[sid] - ocdm) / box)
            print >> out, "C %lf %lf %lf" % (rnow[0], rnow[1], rnow[2])
            rnow = rcs[i] - cdms[sid] + POS_BASE * a1s[i] + (cdms[sid] - ocdm) - box * np.rint ((cdms[sid] - ocdm) / box)
            #rnow = rcs[i] - cdms[strandid[i]] + POS_BASE * a1s[i]
            print >> out, "O %lf %lf %lf" % (rnow[0], rnow[1], rnow[2])
    elif pdb:
        # header
        res = "HEADER    frame t= " + str(times[nconfs])+ " \nMODEL        0 \nREMARK ## 0,0\n" 
        for i in xrange (nnucl):
            sid = strandid[i]
            stringid = strtypes[((sid + 1) % len(strtypes))]

            # define position vector for each nucleotide element
            rnow = rcs[i] - cdms[sid] + (cdms[sid] - ocdm) - box * np.rint ((cdms[sid] - ocdm) / box)
            if os.environ.get(GROOVE_ENV_VAR) == '1':
                s1 = rnow + POS_MM_BACK1 * a1s[i] + POS_MM_BACK2 * a2s[i]
                s2 = rnow
                index_jump = 3
            else:
                s1 = rnow + POS_BACK * a1s[i]
                index_jump = 2
            s3 = rnow + (POS_BACK + 0.68) * a1s[i]

            # magic to get nice ellipse of the base particle
            I_b = np.matrix( ((1.1,0,0),(0,1.5,0),(0,0,2.2)) )
            
            R = np.matrix ((a1s[i], np.cross (a3s[i], a1s[i]), a3s[i]))
            I_l = R.T*I_b*R
            
            anis = np.multiply(-1,I_l)
            anis[0,0] = (I_l[1,1]+I_l[2,2]-I_l[0,0])/2.0
            anis[1,1] = (I_l[0,0]+I_l[2,2]-I_l[1,1])/2.0
            anis[2,2] = (I_l[0,0]+I_l[1,1]-I_l[2,2])/2.0

            # print the backbone site
            res += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (index_jump * i + 1,"A",stringid,'A',i,' ',s1[0],s1[1],s1[2],1,7.895)
            res += "ANISOU%5d %4s %3s %c%4d%c %7i%7i%7i%7i%7i%7i\n" % (index_jump * i + 1,"A",stringid,'A',i,' ',1000,1000,1000,0,0,0)
            # print the centre of mass site (if grooving is switched on)
            if os.environ.get(GROOVE_ENV_VAR) == '1':
              res += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (index_jump * i + 2,"B",stringid,'B',i,' ',s2[0],s2[1],s2[2],1,7.895)
              res += "ANISOU%5d %4s %3s %c%4d%c %7i%7i%7i%7i%7i%7i\n" % (index_jump * i + 2,"B",stringid,'B',i,' ',250,250,250,0,0,0)

            base = base_to_number[basetype[i]]
            if base == 0:
                atomtype = 'O'
            elif base == 1:
                atomtype = 'S'
            elif base == 2:
                atomtype = 'K'
            elif base == 3:
                atomtype = 'P'
            else:
                print >> sys.stderr, "Should not happen..."
                atomtype = 'H'

            # print the base site
            res += "ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (index_jump * i + 3, atomtype, stringid, 'C', i,' ',s3[0],s3[1],s3[2],1,6.316)
            res += "ANISOU%5d %4s %3s %c%4d%c %7i%7i%7i%7i%7i%7i\n" % (index_jump * i + 3, atomtype, stringid, 'C', i, ' ' , anis[0,0]*1000, anis[1,1]*1000, anis[2,2]*1000, anis[0,1]*1000, anis[0,2]*1000, anis[1,2]*1000)

        res += "REMARK  ######### \n\nTER \nENDML \n "
        out.write (res)

    nconfs += 1



inp.close()
out.close()

# generate commands.com to read with UCSF Chimera
if pdb:
    chimera_colors = ["sandy brown", "blue", "red", "green", "yellow", "plum", "sandy brown"]
    strtypes = ["ALA","GLY","CYS","TYR","ARG","PHE","LYS","SER","PRO","VAL","ASN","ASP"]
    mycolors = [chimera_colors[i % len(chimera_colors)] for i in range ( len (strtypes))]
    commands = []
    # cambiamo il colore dello sfondo
    commands.append ("set bg_color white")
    # distruggiamo i legami e li ricreiamo uno per uno per essere
    # sicuri che siano giusti
    commands.append ("~bond #0")
    # make the bonds within each nucleotide
    for i in xrange (nnucl):
        if os.environ.get(GROOVE_ENV_VAR) == '1':
            commands.append ("bond #0:%i.A:%i.B" % (i, i))
            commands.append ("bond #0:%i.B:%i.C" % (i, i))
        else:
            commands.append ("bond #0:%i" % (i))
    # make the bonds between nucleotide backbones
    for i in xrange (nnucl):
        if nn5[i] >= 0:
            commands.append ("bond #0:%i.A,%i.A" % (i, i + 1))
    
    # sistemiamo il colore degli strands (dopo quello delle basi)
    # DOPO aver ricreato i legami
    for i in range (len (strtypes)):
        commands.append ("color %s #0:%s" % (mycolors[i], strtypes[i]))
        commands.append ("col cyan #0:%s@O" % (strtypes[i]))
        commands.append ("col coral #0:%s@S" % (strtypes[i]))
        commands.append ("col yellow #0:%s@K" % (strtypes[i]))
        commands.append ("col cornflower blue #0:%s@P" % (strtypes[i]))
        commands.append ("bondcolor %s #0:%s" % (mycolors[i], strtypes[i]))
    
    # facciamo gli ellissoidi, no? visto che ci siamo
    commands.append ("aniso scale 0.75 smoothing 4")
    
    # sistemiamo la dimensione dei legami
    commands.append ("setattr m stickScale 0.6 #0")

    # e per il momento ci accontentiamo
    f = open ("chimera.com", "w")
    for c in commands:
        #print c
        print >> f, c
    f.close ()

