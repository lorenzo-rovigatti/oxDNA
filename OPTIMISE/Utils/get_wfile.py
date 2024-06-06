#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:37:33 2024

@author: yqb22156
"""

import sys
import math


if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 get_wfile.py nbps")
    sys.exit()

A = 10593.6
nbps = int(sys.argv[1])
rate = 1.72735*8/nbps 


def w(x) :
    w = A*math.exp(-rate*(x-1))
    return (w)


ofile = open("wfile_n"+str(nbps)+".txt","w")


for i in range(nbps+1) :
    if i == 0:
        print("0 8.",file=ofile)
    else:
        print(str(i)+" "+ str(w(i)),file=ofile)