#!/usr/bin/env python

import math
import numpy as np


import copy

def get_Tm(temps,finfs):
	x = copy.deepcopy(finfs)
	y = copy.deepcopy(temps)
	x.reverse()
	y.reverse()
	#print temps
	#print finfs
	xin = np.arange(0.1,1.,0.1)
	f = np.interp(xin,np.array(x), np.array(y))
	xnew = 0.5
	#melttemp = f(xnew)
	#print f
	return f[4]	

def get_Width(temps,finfs,a,b):
	x = copy.deepcopy(finfs)
	y = copy.deepcopy(temps)
	x.reverse()
	y.reverse()
	#print temps
	#print finfs
	xin = np.arange(0.1,1.,0.1)
	f = np.interp(xin,np.array(x), np.array(y))
	
	#melttemp = f(xnew)
	#print xin[1], xin[7]
	#print f
	return f[1] - f[7]	


MAX_N_T = 20


def get_Tm_and_width(temps,histo) :
    
    
    # troviamo il massimo per ogni colonna

    maxx = [m[0] for m in histo]

    for k in range(len(histo)):
        m = histo[k]
        for e in m:
            if e > maxx[k]:
                maxx[k] = e

    for k in range(len(histo)):
        for l in range (len(histo[k])):
            histo[k][l] /= maxx[k]

    '''
    for i in range (nvalues):
        print "%4i" % i, 
        for j in range (ntemps):
            print "% 8.6e" % merged[j][i],
        print "@@"
    '''

    finfs = []
    print( "#### T Phi Phi-extra   ####" )
    i = 0
    for m in histo[0:]:
        s = 0
        for e in m[1:]:
            s += e
        ratio = s/m[0]
        finf = 1. + 1./(2.*ratio) - math.sqrt(   (1. + 1./(2.*ratio))**2  - 1. )
        finfs.append(finf)

        print( "%6.2f %9.7f %9.7f" %(temps[i], s/(m[0] + s), finf ) )
        i += 1

    melttemp = get_Tm(temps,finfs)
    width = get_Width(temps,finfs,0.2,0.8)

    print ("## Tm = ",melttemp, ' = ',melttemp+273.15,'K', ' width = ',width)
    
    return melttemp, width