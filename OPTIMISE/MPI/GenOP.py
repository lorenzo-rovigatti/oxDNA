#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:58:30 2024

@author: yqb22156
"""

import sys

# READ SEQUENCE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 GenOP.py Nbp")
    print("Nbp = Number of base pairs")
    sys.exit()
    
    
Nbp = int(sys.argv[1])    
    
ofile =open("op.txt",'w')

print("{", file=ofile)
print("order_parameter = bond", file=ofile)
print("name = all_native_bonds", file=ofile)
for i in range(Nbp) :
    print("pair"+str(i+1)+ " = "+str(i)+", "+str(2*Nbp-i-1), file=ofile)
print("}", file=ofile)

ofile.close()




