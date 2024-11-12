#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:11:39 2023

@author: yqb22156
"""

import os 
path_opticode=os.path.dirname(os.path.realpath(__file__))


bases = ['A','C','G','T']

Nseq = 0 #number of sequences
seq = []
Njuns = [] #total number of junctions

inj = 0
jfe = 0
in_j = [] #ignore ends
fin_j = [] #ignore ends

# IDS:
# 0,1,2 = buckle, propeller, opening
# 3,4,5 = shear, stretch, stagger
# 6,7,8 = tilt, roll, twist
# 9,10,11 = shift, slide, rise

ids_inter_rot = [6,7,8]

ids = [] #coordinates to optimse
ids_gs = [] #coordinates to optimse. Ground state
ids_cov = [] #coordinates to optimse. Covariance
ids_lrs = [] #coordinates to optimise. Long range stiffness (q=0) #XXXNOTE NOT IMPLEMENTED YET

dimension = [] #total number of coordinates to optimise

Nreps = 1   #number of simulated replicas
in_snap = 0  #ignore first in_snap snapshots (equilibration)

opti_lp = False

ave = False #optimise average or SD
diag = False # diagonal target covariance


algo = "nelder-mead" #default
neva = 40
miter = 10

#params for L-BFGS-B
LBFGSB_eps = 0.01
LBFGSB_iprint = 1

Niter = 0

T = 0.1 #300K

weight_gs = 1.

modelh = path_opticode+"/model.h"
