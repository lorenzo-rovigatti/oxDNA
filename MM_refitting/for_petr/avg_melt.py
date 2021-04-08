#!/usr/bin/env python

import random

bases = ['A','C','G','T']
inverted = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A' }

def GenRanSeq(length):
        sek = ''
        compsek = ''
        for j in range(length):
                compsek += '*'
        for i in range(length):
                base = random.choice(bases)
                sek += base
                idx = length - 1  -i
                compsek = compsek[0:idx] + inverted[base] + compsek[idx+1:]
                #compsek[length - 1 - i] = inverted[base]
        return sek,compsek


import sys
import math

if len(sys.argv) != 3:
	print 'Usage: ./program LENGTH box'
	sys.exit(1)

box = float(sys.argv[2])
length = int(sys.argv[1]) - 1



molcon = (2/(8.5179*10.0**(-9)*box)**3)/(6.0221415*10**23)
saltcorr =  0.368*(length ) * math.log(0.5)

DH = (-8.2375*length + 1.2*2)
DS = (-22.0188*length + 1.2) + saltcorr

RealTemp = (DH )*1000 / (DS + 1.9859* math.log(molcon/4) )

c, s = GenRanSeq(length+1)

print c,s,RealTemp
