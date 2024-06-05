#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 17:02:40 2024

@author: yqb22156
"""

import sys
import numpy as np


def estimate_new_wfile(hist,wfile_name) :
    

    max_entry = max(hist)
    w = []
    for i in hist :
        if i == 0:
            print("WARNING! There is a zero in the w file. Can't produce final w file")
            return False    
        w.append(max_entry/i)
        
    ofile = open(wfile_name,'w')
    
    for i in range(len(w)):
        print(str(i) + " " + str(w[i]),file=ofile)
        
    ofile.close()
    
    return True


if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 new_wfile.py Nreps")
    sys.exit()
    
Nreps = int(sys.argv[1])

hist0 = np.loadtxt("Rep0/last_hist.dat")

print(hist0)

hist = hist0[:,2]

print(hist)

for i in range(1,Nreps) :
    hist += np.loadtxt("Rep"+str(i)+"/last_hist.dat")[:,2]
    
print(hist)

wfile_name = "new_weights.txt"
estimate_new_wfile(hist,wfile_name)

    