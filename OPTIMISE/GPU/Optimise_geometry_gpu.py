
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 12:04:54 2024

@author: yqb22156
"""

from oxdna_to_internal_wflip import read_oxdna_trajectory_standard_order
import parameters_list as parl
import torch
import functions as fun
import cost_function as cfun
import config as cg
import sys
import time

# READ CONFIG FILE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 optimise.py config_file")
    sys.exit()
    
start_time = time.time()

config_file = sys.argv[1]

PARS_LIST = parl.PARS_LIST
par_index = parl.par_index


###################################################################################################
############## INITIALISE PARAMETERS FROM model.h E SD PARAMETERS FILE ############################
###################################################################################################

#print(stck_fact_eps)
model_file = open("model.h",'r')
pars_from_modelh, vals_from_modelh = fun.read_vanilla_parameters(model_file)    
model_file.close()

SD_par_file = open("oxDNA_sequence_dependent_parameters_in.txt",'r')   
over_indices, over_vals, stck_fact_eps_read, stck_fact_eps = fun.read_pars_from_SD_file(SD_par_file)    
SD_par_file.close()
    
if stck_fact_eps_read :
    print("STCK_FACT_EPS read from SD parameters file")
else:
    print("WARNING: No STCK_FACT_EPS found in SD parameters file")

OXPS_zero, shifts = fun.init_oxpars(pars_from_modelh, vals_from_modelh, over_indices, over_vals, cg.T, stck_fact_eps)

###################################################################################################
############## READ OPTIM OPTIONS  ################################################################
###################################################################################################

#READ PARAMETERS
if fun.read_config(config_file) == False :
    sys.exit()

#########################################################################################################################
############## READ TRAJECTORY, COMPUTE OXDNA COORDINATES (i.e angles and distances) AND INTERNAL COORDINATES ###########
#########################################################################################################################


fene_r = []
stck_r = []
th4_bn = []
th5 = []
th6 = []
cosphi1 = []
cosphi2 = []
types_bn = []
hydr_r = []
th1 = []
th2 = []
th3 = []
th4_unbn = []
th7 = []
th8 = []
types_unbn = []

rclow, rchigh = fun.find_cuts_for_lists(OXPS_zero)

nseq = 2
nreps = 1

for l in range(nseq):
    for m in range(nreps) :
        tr_file = open("trajectory.dat",'r')
        topo_file = open("generated.top", 'r')
        
        
        #oxdna distances, types and angles
        fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, tr_file, topo_file)
        
        if m == 0:
            fene_r.append(fr)
            stck_r.append(sr)
            th4_bn.append(t4bn)
            th5.append(t5)
            th6.append(t6)
            cosphi1.append(cp1)
            cosphi2.append(cp2)
            types_bn.append(tbn)
            hydr_r.append(hr)
            th1.append(t1)
            th2.append(t2)
            th3.append(t3)
            th4_unbn.append(t4un)
            th7.append(t7)
            th8.append(t8)
            types_unbn.append(tun)
        else:
            fene_r[l].extend(fr)
            stck_r[l].extend(sr)
            th4_bn[l].extend(t4bn)
            th5[l].extend(t5)
            th6[l].extend(t6)
            cosphi1[l].extend(cp1)
            cosphi2[l].extend(cp2)
            types_bn[l].extend(tbn)
            hydr_r[l].extend(hr)
            th1[l].extend(t1)
            th2[l].extend(t2)
            th3[l].extend(t3)
            th4_unbn[l].extend(t4un)
            th7[l].extend(t7)
            th8[l].extend(t8)
            types_unbn[l].extend(tun)
        
        tr_file.close()
        topo_file.close()
        
        tr_file = open("trajectory.dat",'r')
        topo_file = open("generated.top", 'r')
        
        #internal coordinates
        traj = read_oxdna_trajectory_standard_order(tr_file, topo_file)
        
        if l == 0 :
            fun.store_internal_coord(traj,l,cg.ids,fun.in_j[l],fun.fin_j[l],cg.in_snap,True)
        else  :
            fun.store_internal_coord(traj,l,cg.ids,cg.in_j[l],cg.fin_j[l],cg.in_snap,False)
        
        tr_file.close()
        topo_file.close()


###################################################################################################
############## SETUP TENSORS FOR COST FUNCTION ####################################################
###################################################################################################

#create all tensors on the gpu. Change this to easily swap between gpu and cpu
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

cfun.init_tensors(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                  th4_unbn, th7, th8, types_unbn, shifts, OXPS_zero)

 
#build tensors imposing continuity
cfun.build_continuity_tensors()
#build masks (for selecting optim parameters and marginalising mu and cov) and symm tensors (for imposing symmetries)
cfun.build_masks_and_symm_tensors()
#create reduced targets and compute cov^{-1}
cfun.reduce_targets()
#compute initial energy and modulation factors
cfun.compute_initial_energy()



timefile = open("runtime.txt",'w')
print("Run time [s]: " + str(time.time() - start_time),file=timefile)
timefile.close()