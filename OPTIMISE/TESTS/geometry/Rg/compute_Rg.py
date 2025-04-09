#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 12:46:13 2023

@author: yqb22156
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad

#from oxdna_to_internal_triadII import read_oxdna_trajectory_standard_order #Carlon style (see 2017 J. chem. phys)
from oxdna_to_internal_cgna import read_oxdna_trajectory_standard_order #Maddocks style (see theses)

in_j = 0 #ignore ends
fin_j = Njuns -0 #ignore ends
in_snap = 0 #ignore first in_snap snapshots (equilibration)

#read oxdna trajectory and compute internal coordinates
iname = 'trajectory.dat'
tname = 'generated.top'

ifile = open(iname,'r')
tfile = open(tname,'r')


traj = read_oxdna_trajectory_standard_order(ifile, tfile)

ifile.close()
tfile.close()

#Compute Rg and end to end

ave_r_v_t = []
ave_r2_v_t = []
ave_ete_v_t = []
#ave_ete2_v_t = []

for i in range(len(trajectory)) :
    ave_r = 0.
    ave_r2 = 0.
    size = len(trajectory[i])
    for j in range(size) :
    	r = np.linalg.norm(trajectory[i][j].base_pair1.pos)
        ave_r += r/size
        ave_r2 += r*r/size
    ave_r_v_t.append(ave_r)
    ave_r2_v_t.append(ave_r2)
    r = np.linalg.norm(trajectory[i][size-1].base_pair1.pos - trajectory[i][0].base_pair1.pos)
    ave_ete_v_t.apend(r)
    #ave_ete2_v_t.apend(r*r)
    
ofile = open("Rg_ete.txt","w")
for m in range(len(ave_r_v_t)) :
    Rg = math.sqrt(ave_r2_v_t[m]-ave_r_v_t[m]*ave_r_v_t[m])
    print(m,Rg,ave_ete_v_t[m],file=ofile)
    
ofile.close()
