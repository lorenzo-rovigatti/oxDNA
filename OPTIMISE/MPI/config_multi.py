#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:11:39 2023

@author: yqb22156
"""

from mpi4py import MPI

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


par_codename = [] #names of the parameters to use for the optimisation
continuity_par_codenames = [] #auxiliary parameters needed for imposing continuity
continuity_par_values = [] #values of these parameters
used = [] #track used parameters when imposing continuity (avoid printing duplication)
par_dimension = 0


internal_coords = [] #stores arrays of internal coordinates for sampled configurations
energy_sampled = [] #stores the energy of the sampled configurations
energy_stk_sampled = [] #stores the stacking enegy of the sampled configurations

mu_sampled = [] #sampled means (list of arrays, one for each different sequence)
cov_sampled = [] #sampled covaraiances (list of matrices, one for each sequence)

target_mu = [] #target means
target_cov = [] #target covariances

opti_lp = False

target_lb = 45.0 #target average long range bending persistence length (in nm).
target_lt = 220.0 #target average long range torsional persistence length (in nm).

target_C = 140 #target average long range C modulus
target_Ar = 40 #target average long range Ar (A2) modulus


weight_gs = 1.

ave = False #optimise average or SD
diag = True # diagonal covariance


algo = "nelder-mead" #default
neva = 40
miter = 10

#params for L-BFGS-B
LBFGSB_eps = 0.01
LBFGSB_iprint = 1

Niter = 0


#MPI

comm = MPI.COMM_WORLD
comm_seq = None
leaders = []
comm_leaders = None
rank = 0
rank_seq = 0
rank_leaders = -1
seq_id = 0
rep_id = 0

update_rews = False

############# FOR MELTING TEMPERATURE ############


sampled_op_histo_b = [] #order parameter histogram biased
sampled_op_histo_un = [] #order parameter histogram unbiased
sampled_op_hist_all_Ts = [] #order parameter histogram for all temperatures


Diter_Trange = 0
nevery_uprewTs = 1

Ct = [] #total single strand concentration (in M (mol/L)) -- for a single double helix: 2/box_size^3*2.6868 M, where the 2 accounts for the two strands. 
Cs = [] #salt concentartion (in M)

target_mTs = []
current_mT = 0.


simTs = []   #simulated temperature
n_Ts = []  #rew. temperatures
DTs = []  #delta T between rew temperatures
rew_Ts = [] #reweighted temperatures
weights = []    #weights

hbs_sampled = []



