#!/usr/bin/env python

import sys, math, getopt
import math
import numpy as np


import copy

def get_Tm(temps,finfs):
	x = copy.deepcopy(finfs)
	y = copy.deepcopy(temps)
	x.reverse()
	y.reverse()
	xin = np.arange(0.1,1.,0.1)
	f = np.interp(xin,np.array(x), np.array(y))
	xnew = 0.5
	return f[4]	

def get_Width(temps,finfs,a,b):
	x = copy.deepcopy(finfs)
	y = copy.deepcopy(temps)
	x.reverse()
	y.reverse()
	xin = np.arange(0.1,1.,0.1)
	f = np.interp(xin,np.array(x), np.array(y))
	
	return f[1] - f[7]	
	

## option parsing
shortArgs = []
longArgs = []

args, files = getopt.getopt(sys.argv[1:], shortArgs, longArgs)

for file in files:
    try:
        inp = open(file, 'r')
        inp.close()
    except:
        print('Could not open file', file, 'Aborting.', file=sys.stderr)
        sys.exit()


MAX_N_T = 20

op_dim = 1
merged = []
partials = []
ratios = []


inp = open (files[0], 'r')
nvalues = len (inp.readlines ()) - 1
inp.seek(0)
inp.readline ()
ntemps = len (inp.readline().split()) - (op_dim + 1)
inp.seek(0)


histo = [[0 for i in range(ntemps)] for j in range(nvalues)]
guessed_t = False

for file in files:
    inp = open (file, 'r')
    partial = [[] for i in range (ntemps)]
    
    for line in inp.readlines():
        if not line.startswith('#'):
            words = line.split()
            #print len(words)
            for i in range (0, len (words) - (op_dim + 1)):
                partial[i].append (float (words[op_dim + 1 + i]))
        elif not guessed_t:
            temps = [float(x) * 3000. - 273.15 for x in line.split()[-(ntemps - 1):]] 
            guessed_t = True
    
    inp.close()
    partials.append (partial)

merged = [[0. for j in range (len (partials[0][0]))] for i in range (len (partials[0]))]

for i in range (len(partials)):
    for j in range (len(partials[i])): # j on order parameter value
        for k in range (len (partials[i][j])):
            merged[j][k] += partials[i][j][k]


# find the maximum for each column
max = [m[0] for m in merged]

for k in range(len(merged)):
    m = merged[k]
    for e in m:
        if e > max[k]:
            max[k] = e

for k in range(len(merged)):
    for l in range (len(merged[k])):
        merged[k][l] /= max[k]

'''
for i in range (nvalues):
    print "%4i" % i, 
    for j in range (ntemps):
        print "% 8.6e" % merged[j][i],
    print "@@"
'''

finfs = []
print("#### T Phi Phi-extra   ####")
i = 0
for m in merged[1:]:
    s = 0
    for e in m[1:]:
        s += e
    ratio = s/m[0]
    finf = 1. + 1./(2.*ratio) - math.sqrt(   (1. + 1./(2.*ratio))**2  - 1. )
    finfs.append(finf)

    print("%6.2f %9.7f %9.7f" %(temps[i], s/(m[0] + s), finf ))
    i += 1

melttemp = get_Tm(temps,finfs)
width = get_Width(temps,finfs,0.2,0.8)

print("## Tm =",melttemp, 'C =',melttemp+273.15,'K, width =',width, "K")

sys.exit()

