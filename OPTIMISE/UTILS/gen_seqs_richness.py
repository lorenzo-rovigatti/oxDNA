#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 12:43:04 2024

@author: yqb22156

generates a random sequence with A_r, G_r, C_r and T_r richness


"""
import sys
import numpy as np
import random


if len(sys.argv) != 6 :
    print("Invalid syntax.")
    print(sys.argv[0]+ " Nbps A_rich G_rich C_rich T_rich")
    exit(1)
    
Nbps = int(sys.argv[1])
A_r = float(sys.argv[2])
G_r = float(sys.argv[3])
C_r = float(sys.argv[4])
T_r = float(sys.argv[5])

norm = A_r + G_r + C_r + T_r

A_r = A_r/norm
G_r = G_r/norm
C_r = C_r/norm
T_r = T_r/norm

seq = ""

counts = [0,0,0,0]

for i in range(Nbps) :
    rf = random.random()
    if rf < A_r :
        seq+='A'
        counts[0] += 1
    elif rf < A_r+G_r :
        seq+='G'
        counts[1] += 1
    elif rf < A_r+G_r+C_r:
        seq += 'C'
        counts[2] += 1
    else :
        seq += 'T'
        counts[3] += 1
    
ofile = open("seq.txt",'w')
print("double "+seq,file=ofile)
ofile.close()
#print("Base types counts:")
#print(counts)
    


