#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:11:39 2023

@author: yqb22156
"""

import os
path_opticode=os.path.dirname(os.path.realpath(__file__)) #path of this file

print_energy_to_file = False #True if we want to print energy to file
print_coords_to_file = False #True if we want to print coordinates to file
read_energy_from_file = True #True if we want to read initial energy from file. Useful because par0 might not be the sampled parameters (e.g. in geometry-melting cycle)
read_coords_from_file = True #True if we want to read coords from file. Useful to avoid computing coordinates over and over in a geometry+melting cycle

bases = ['A','C','G','T']

Nseq_n5 = 0 #number of sequences
Nseq_n8 = 0 #number of sequences
Nseq_n15 = 0 #number of sequences
seq_n5 = [] #sequences
seq_n8 = []
seq_n15 = []

#These are useful for geometry only

Njuns_n5 = [] #total number of junctions
Njuns_n8 = [] #total number of junctions
Njuns_n15 = [] #total number of junctions

inj = 0
jfe = 0
in_j_n5 = [] #ignore ends
fin_j_n5 = [] #ignore ends
in_j_n8 = [] #ignore ends
fin_j_n8 = [] #ignore ends
in_j_n15 = [] #ignore ends
fin_j_n15 = [] #ignore ends

# IDS:
# 0,1,2 = buckle, propeller, opening
# 3,4,5 = shear, stretch, stagger
# 6,7,8 = tilt, roll, twist
# 9,10,11 = shift, slide, rise

ids = [] #coordinates to optimse
ids_gs = [] #coordinates to optimse. Ground state
ids_cov = [] #coordinates to optimse. Covariance
ids_lp = [6,7,8] #coordinates to optimise. Persistence length

dimension = [] #total number of coordinates to optimise

Nreps = 1   #number of simulated replicas
in_snap = 0  #ignore first in_snap snapshots (equilibration)

opti_lp = False
lp_comp_range = None #this is m in Enrico's paper
#lp_m = 3 #lp_comp_range = used junctions - lp_m -- OLD VERSION
lp_m = 11 #NEW VERSION. 11 should give a good estimation of lb and lt, both underestimated by about 5%.
lb_corr = 1.03
lt_corr = 1.15*1.05

ave = False #optimise average or SD
diag = False # diagonal target covariance


#these are not used
algo = "nelder-mead" #default
neva = 40
miter = 10

#params for L-BFGS-B
LBFGSB_eps = 0.01
LBFGSB_iprint = 1

Niter = 0

T = 0.1 #Temeperature is initialised at 300K. The correct temperature is accounted by rescaling stck_eps (and betha)

weight_gs = 1.

modelh = path_opticode+"/model.h" #model.h file

############# FOR MELTING TEMPERATURE ############

stck_fact_eps = 0.18

parallel_tempering = False #True if parallel_tempering was used

N_PT_reps_n5 = 8 #number of pt replicas n5
N_PT_reps_n8 = 8 #number of pt replicas n8
N_PT_reps_n15 = 8 #number of pt replicas n15

Nconfs_per_pt_rep = 0 #sampled confs per pt replica - in_snap
tot_time = 0 #for pt, this is the sim time of each replica
delta_time = 0 #time interval between snapshots
delta_print_en = 0 #time interval between printing energy
delta_split_en = 0 #time interval between printing split_energy

Diter_Trange = 0 #boh
nevery_uprewTs = 1 #update reweighting Ts every n iteration

Ct_n5 = [] #total single strand concentration (in M (mol/L)) -- for a single double helix: 2/box_size^3*2.6868 M, where the 2 accounts for the two strands.
Cs_n5 = [] #salt concentartion (in M)
Ct_n8 = [] #total single strand concentration (in M (mol/L)) -- for a single double helix: 2/box_size^3*2.6868 M, where the 2 accounts for the two strands.
Cs_n8 = [] #salt concentartion (in M)
Ct_n15 = [] #total single strand concentration (in M (mol/L)) -- for a single double helix: 2/box_size^3*2.6868 M, where the 2 accounts for the two strands.
Cs_n15 = [] #salt concentartion (in M)

boxes_n5 = [] #box size
boxes_n8 = []
boxes_n15 = []

current_mT_n5 = 0. #melting Temperature with current parameters
current_mT_n8 = 0.
current_mT_n15 = 0.
#current_hist_Ts = []

first_step_flag = False #Not used, I think.
#fin_wfiles = []

symm_stck = False #Not used anymore, I think.
