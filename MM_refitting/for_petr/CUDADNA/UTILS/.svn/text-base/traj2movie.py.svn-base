#!/usr/bin/env python

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import getopt
import subprocess as sp

def print_usage():
    print "USAGE:"
    print "\t%s configuration topology [--help] [--fps=3] [--width=500]" % sys.argv[0]
    print "\t[--heigh=500] [--crepy=crepy.py] [--crepyopts=''] [--skip=0]\n"
    sys.exit(1)

shortArgs = ''
longArgs = ['help', 'fps=', 'width=', 'height=', 'crepy=', 'crepyopts=', 'skip=']

skip = 0
fps = 3
width = 500
height = 500
crepy = "crepy.py"
crepyopts = ""

try:
    args, files = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)

    if len(files) != 2:
        print_usage ()

    for k in args:
        if k[0] == '--help': print_usage()
        elif k[0] == '--fps': fps = int(k[1])
        elif k[0] == '--width': width = int(k[1])
        elif k[0] == '--height': height = int(k[1])
        elif k[0] == '--crepy': crepy = k[1]
        elif k[0] == '--crepyopts': crepyopts = k[1]
        elif k[0] == '--skip': skip = int(k[1])
        else: pass

except Exception:
    print_usage()

l = readers.LorenzoReader(files[0], files[1])
s = l.get_system(N_skip=skip)
append = False
niter = 1
pngfiles = []
devnull = open("/dev/null", "w")
while s:
    base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
    s.print_ribbon_output("__ggg__.mgl", same_colors=True)
    
    a = sp.Popen ([crepy, "__ggg__.mgl", "--output=__ggg__.pov", "-o", crepyopts], stdout=devnull)
    a.wait()
    if a.returncode != 0:
        base.Logger.die ("crepy e' morto")
    
    pngfile = "__%s.png" % niter
    a = sp.Popen (["povray", "+A", "+r", "-D", "+H%s" % (str(height)), "+W%s" % (str(width)),"+O%s" % pngfile, "__ggg__.pov"], stdout=devnull, stderr=devnull)
    a.wait()
    if a.returncode != 0:
        base.Logger.die ("povray e' morto")

    pngfiles.append(pngfile)

    s = l.get_system(N_skip=skip)
    append = True
    niter += 1

devnull.close ()

base.Logger.log("Enconding movie...")
# costruiamo la lista
f = open ("__ginocchio__", "w")
for p in pngfiles:
    print >> f, p
f.close ()

print "mencoder mf://@__ginocchio__ -mf w=%s:h=%s:fps=%s:type=png -ovc xvid -xvidencopts bitrate=200 -o output.avi" % (str(width), str(height), str(fps))

#b = sp.Popen ("/usr/bin/mencoder mf://@__ginocchio__ -mf w=%s:h=%s:fps=%s:type=png -ovc xvid -xvidencopts bitrate=200 -o output.avi" % (str(width), str(height), str(fps)))
b = sp.Popen (["mencoder", "mf://@__ginocchio__", "-mf","w=%s:h=%s:fps=%s:type=png" % (str(width), str(height), str(fps)), "-ovc", "xvid", "-xvidencopts", "bitrate=200", "-o", "output.avi"])
#, "-ovc xvid -xvidencopts bitrate=200 -o output.avi"])
b.wait()


