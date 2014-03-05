#!/usr/bin/env python
# this script just reads the energy line in a trajectory and spits 
# something similar to an energy.dat file

import sys

if len (sys.argv) < 2:
    print >> sys.stderr, "Usage: ", sys.argv[0], "<trajectory_file>"
    sys.exit (-1)

try:
    inp = open (sys.argv[1], 'r')
except:
    print >> sys.stderr, "Could not open ", sys.argv[1], "Aborting"
    sys.exit (-2)

info = [[]]
line = inp.readline()
while line:
    if line.startswith ('t'):
        words = line.split ("=")
        t = float(words[1])
        line = inp.readline ()
        line = inp.readline ()
        words = line.split ("=")
        e = float(words[1].split()[0])
        info.append ([t, e])
    line = inp.readline ()

for t in info:
    try:
        print t[0], t[1]
    except:
        if t:
            print >> sys.stderr, "skipped", t
    # it t is empty, we dont care


