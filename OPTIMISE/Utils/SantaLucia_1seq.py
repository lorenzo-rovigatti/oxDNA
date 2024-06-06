#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:03:08 2024

@author: yqb22156
"""
import sys
import SantaLucia as SL


# READ SEQUENCE
if len(sys.argv) != 4 :
    print("Unknown argument format.")
    print("Usage: python3 SantaLucia.py sequence Ct Cs")
    print("Ct = total single strand concentration in M")
    print("Cs = salt concentration in M")
    print("For 1 duplex in a box of size l ox units, Ct = 2/l^3*2.6868 M")
    sys.exit()
    

seq = sys.argv[1]
Ct = float(sys.argv[2])
Cs = float(sys.argv[3])

Tm = SL.melting_temperature(seq,Ct,Cs)

print("Tm: " +str(Tm))