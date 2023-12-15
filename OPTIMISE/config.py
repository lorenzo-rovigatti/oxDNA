#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:11:39 2023

@author: yqb22156
"""

bases = ['A','C','G','T']

seq = ""
Njuns = 0 #total number of junctions

in_j = 0 #ignore ends
fin_j = 0 #ignore ends

# IDS:
# 0,1,2 = buckle, propeller, opening
# 3,4,5 = shear, stretch, stagger
# 6,7,8 = tilt, roll, twist
# 9,10,11 = shift, slide, rise

ids_inter_rot = [6,7,8]

ids = [] #coordinates to optimse
ids_gs = [] #coordinates to optimse. Ground state
ids_cov = [] #coordinates to optimse. Covariance
ids_lrs = [] #coordinates to optimise. Long range stiffness (q=0)

dimension = 0 #total number of coordinates to optimise

Nreps = 1   #number of simulated replicas
in_snap = 0  #ignore first in_snap snapshots (equilibration)


par_codename = [] #names of the parameters to use for the optimisation
continuity_par_codenames = [] #auxiliary parameters needed for imposing continuity
continuity_par_values = [] #values of these parameters
used = [] #track used parameters when imposing continuity (avoid printing duplication)
par_dimension = 0


internal_coords = [] #stores arrays of internal coordinates for sampled configurations
energy_sampled = [] #stores the enegy of the sampled configurations

mu_sampled = [] #sampled mean
cov_sampled = [] #sampled covaraiance

target_mu = [] #target mean
target_cov = [] #target covariance
target_lb = 45.0 #target bending persistence length (in nm).
target_lt = 220.0 #target torsional persistence length (in nm). Note: long range lt is about the same 

target_C = 140 #target bending persistence length (in nm).
target_Ar = 40 #target torsional persistence length (in nm). Note: long range lt is about the same 


ave = True #optimise average or SD
diag = True # diagonal covariance

miter = 10

