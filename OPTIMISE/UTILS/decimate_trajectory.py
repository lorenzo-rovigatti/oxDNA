#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 12:01:30 2024

@author: yqb22156
"""

import sys


# READ CONFIG FILE
if len(sys.argv) != 4 :
    print("Unknown argument format.")
    print("Usage: python3 decimate_trajectory.py trajectory_file nevery nbps")
    sys.exit()
    

ifile = open(sys.argv[1], 'r')
nevery = int(sys.argv[2])
nbps = int(sys.argv[3])

ofile = open("trajectory_decimated.dat", 'w')



counter = -1
counter_conf = -1
for line in ifile.readlines() :
    counter+=1
    
    counter_conf = int(counter / (nbps*2+3))
    if counter_conf%nevery == 0:
        print(line,file=ofile,end='')
        
        
ifile.close()
ofile.close()
        