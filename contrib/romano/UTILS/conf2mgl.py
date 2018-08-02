#!/usr/bin/env python

import plpatchy as plp
import sys
import os

if len(sys.argv) < 3:
    print 'Usage: %s input_file configuration_file [conf_id=0] [outname]' % (sys.argv[0])
    sys.exit(-1)
    
    
input_file = sys.argv[1]
conf_file = sys.argv[2]
topology_file = None


if os.getenv("PL_THICK") != None:
       thick_balls = True
else:
       thick_balls = False

if len(sys.argv) >= 4:
    conf_id = int(sys.argv[3])
else:
    conf_id = 0


if len(sys.argv) >= 5:
    outname = sys.argv[4]
else:
    outname = conf_file + '.mgl'

icosahedron = False
inf = open(input_file,'r')
for line in inf.readlines():
    line = line.strip()
    if len(line) > 1 and line[0] != '#':
        if 'topology' in line:
            topology_file = line.split('=')[1].strip()
        if 'shape' in line:
            icotype = line.split('=')[1].strip()
            if icotype == 'icosahedron':
                icosahedron = True

inf.close()

print >> sys.stderr, 'Input file: %s, topology_file: %s, conf_file: %s, output_file: %s' % (input_file,topology_file,conf_file,outname)

simulation = plp.PLPSimulation()
simulation.load_from_files(input_file,topology_file,conf_file,conf_id)

#print simulation.particles[0].cm_pos
#print simulation.particles[0].a1,  simulation.particles[0].a2,  simulation.particles[0].a3

if not thick_balls:
    simulation.set_radius(0.3)
simulation.bring_in_box(True)
simulation.export_to_mgl(outname,'w',icosahedron)


