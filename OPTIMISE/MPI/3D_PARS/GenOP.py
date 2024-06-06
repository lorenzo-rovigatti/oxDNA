#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:58:30 2024

@author: yqb22156
"""

import sys

# READ SEQUENCE
if len(sys.argv) != 3 :
    print("Unknown argument format.")
    print("Usage: python3 GenOP.py Nbp config_file")
    print("Nbp = Number of base pairs")
    sys.exit()
    
    
    
    
cfile = open(sys.argv[2],'r')


in_j = 0
j_from_end = 0

for line in cfile.readlines() :
    vals = line.split()
    if len(vals) == 0:
        continue
    if vals[0][0] == '#':
        continue   
    if(vals[0] == 'IN_J'):
        in_j = int(vals[1])
    if(vals[0] == 'J_FROM_END'):
        j_from_end = int(vals[1])

    
    
if in_j > 5 :
    in_j = in_j - 3
if j_from_end > 5 :
    j_from_end = j_from_end- 3
    
Nbp = int(sys.argv[1])    
    
ofile =open("op.txt",'w')

print("{", file=ofile)
print("order_parameter = bond", file=ofile)
print("name = all_native_bonds", file=ofile)
for i in range(in_j,Nbp-j_from_end) :
    print("pair"+str(i-in_j+1)+ " = "+str(i)+", "+str(2*Nbp-i-1), file=ofile)
print("}", file=ofile)

ofile.close()




