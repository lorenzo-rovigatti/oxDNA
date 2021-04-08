#!/usr/bin/env python

import sys
import base
import readers
import numpy as np
import origami_utils as oru

trim = 3
up_to = 1e99

if (len(sys.argv) < 3):
    print 'Usage %s trajectory topology [output]' % sys.argv[0]
    sys.exit()

if len(sys.argv) > 3:
    fnradius = sys.argv[3]+"_radius.dat"
    fnrise = sys.argv[3]+"_rise.dat"
    fnradius_ser = sys.argv[3]+"_rad_series.dat"
else:
    fnradius = "radius.dat"
    fnrise = "rise.dat"
    fnrise_ser = "rad_series.dat"

reader = readers.LorenzoReader(sys.argv[1], sys.argv[2])
ss = reader.get_system()
origami = oru.Origami(ss, "virt2nuc")

frise_ser = open(fnrise_ser, "w")
radius_av = [[0 for jj in range(len(origami.vvib[ii]))] for ii in range(len(origami.vhelix_indices))]
rise_av = [[0 for jj in range(len(origami.vvib[ii])-1)] for ii in range(len(origami.vhelix_indices))]
conf_counter = 0
while ss:
    conf_counter += 1
    base.Logger.log("reading conf %d" % conf_counter, base.Logger.INFO)

    this_rad_av = 0
    for ii in range(len(origami.vhelix_indices)):
        for jj in range(len(origami.vvib[ii])):
            if not jj == len(origami.vvib[ii])-1: # rise is calculated between bps whereas radius is calculated for each bp
                rise_av[ii][jj] += origami.get_rise(ii, jj)
            this_radius = origami.get_radius(ii, jj)
            radius_av[ii][jj] += this_radius
            this_rad_av += this_radius

    frise_ser.write("%f\n" % this_rad_av)
    if conf_counter > up_to:
        break
    ss = reader.get_system()
    origami.update_system(ss)

frise_ser.close()
# find average
rise_av = [[yy/conf_counter for yy in xx] for xx in rise_av]
radius_av = [[yy/conf_counter for yy in xx] for xx in radius_av]

# print result
frad = open(fnradius, "w")
for ii in range(len(radius_av)):
    for jj in range(len(radius_av[ii])):
        frad.write("%d %d %f\n" %(ii,  jj, radius_av[ii][jj]))
frad.close()
base.Logger.log("wrote to file %s" % fnradius, base.Logger.INFO)

frise = open(fnrise, "w")
for ii in range(len(rise_av)):
    for jj in range(len(rise_av[ii])):
        frise.write("%d %d %f\n" %(ii,  jj, rise_av[ii][jj]))
frise.close()
base.Logger.log("wrote to file %s" % fnrise, base.Logger.INFO)

radavav = sum(radius_av[0][trim:-trim])/len(radius_av[0][trim:-trim])
riseavav = sum(rise_av[0][trim:-trim])/len(rise_av[0][trim:-trim])
print "sim units:", radavav, riseavav
print "nm:", radavav*0.8518, riseavav*0.8518
