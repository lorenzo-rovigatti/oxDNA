import sys
try:
    import numpy as np
except:
    import mynumpy as np
import base
import generators as gen

# definiamo l'origami su una scacchiera...
width = 4
height = 8
#checkboard = np.matrix ([[ 0 for i in range(width)] for i in range(height)], np.int32)

checkboard = [[1 for i in range (width)] for i in range (height)]
checkboard[0] = [0, 1, 1, 0]

# major index on heigth
# each height step is a dsDNA
# each width step is a 16bp turn
print >> sys.stderr, "#(stderr) Will generate a %d nucleotide system" % (width * height * 32)

print checkboard

# for each row, we print out a dsDNA
box = np.max([2.5 * height, 17 * width])
box = 50

s = base.System([box, box, box])

g = gen.StrandGenerator()

for i in range (height):
    for j in range (width):
        if checkboard[i][j] > 0:
            print "#",
            # orientation along positive x if i is odd
            align = [float (pow (-1, i)), 0, 0]
            #align_perp = [0, 0, float (pow (-1, i))]
            align_perp = [0, 1, 0]
            rot = [j * 16, gen.BP]

            # start position of strand
            pos = [j * 16 * base.BASE_BASE, -i * (2. * base.RC2_BACK + base.CM_CENTER_DS + 0.01), 0]
            
            # fix starting position if align = [-1, 0, 0]
            if align[0] < 0:
                pos[0] += 16 * base.BASE_BASE;
                rot[0] = width * 16 - rot[0]
            #print pos

            s.add_strands(g.generate(16, dir=align, perp=align_perp, start_pos=pos, rot=rot), check_overlap=False)
        else:
            print "0",
    print

s.print_vmd_xyz_output (same_colors=False)
s.print_crepy_output ("caca.mgl", same_colors=False)
s.print_crepy_output ("caca2.mgl", same_colors=True)

def max (a, b):
    if a > b:
        return a
    else:
        return b

box = 4 * max(width, height)
s2 = base.System([box, box, box])
npoints = 0
for i in range (height):
    for j in range (width):
        if checkboard[i][j] == 1:
            npoints += 1

s2.add_strands(g.generate (npoints, dir=[0,1,0], double = False), check_overlap=False)

strand = s2._strands[0]
now = 0
for i in range (height):
    for j in range (width):
        if checkboard[i][j] == 1:
            n = strand._nucleotides[now]
            print "caca"
            print "prima", n.cm_pos
            n.cm_pos = [3 * j, - 3 * i, 0]
            n.cm_pos_box = n.cm_pos
            print "dopo", n.cm_pos
            now = now + 1



for n in strand._nucleotides:
    print n.cm_pos_box
    
s2.print_crepy_output ("ricaca.mgl", same_colors=False)

