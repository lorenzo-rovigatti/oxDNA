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
import numpy as np
from scipy import optimize


# READ CONFIG FILE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 optimise.py config_file")
    sys.exit()

start_time = time.time()

config_file = sys.argv[1]

PARS_LIST = parl.PARS_LIST
par_index = parl.par_index


#X = torch.tensor([0,1,2])
#print(X.device)


###################################################################################################
############## INITIALISE PARAMETERS FROM model.h E SD PARAMETERS FILE ############################
###################################################################################################

print(torch.cuda.memory_summary(device=None, abbreviated=False))
torch.cuda.empty_cache()

#print(stck_fact_eps)

model_file = open(cg.modelh,'r')
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

"""
test_file = open("opti_p_test.txt",'w')
for l in range(len(cfun.OPT_PAR_LIST)) :
    print(cfun.OPT_PAR_LIST[l],file=test_file)

test_file.close()

test_file = open("opti_p_test1.txt",'w')
print(OXPS_zero[45],file = test_file)
print(OXPS_zero[1],file = test_file)

print(OXPS_zero[34],file = test_file)
test_file.close()
"""

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
types_unbn_33 = []
types_unbn_55 = []

rclow, rchigh = fun.find_cuts_for_lists(OXPS_zero)
print("cuts: "+str(rclow)+" "+str(rchigh))

cg.Nreps = 1
#fun.e3s = [[]*cg.Nseq]
for l in range(cg.Nseq): fun.e3s.append([])
#print(fun.e3s)

for l in range(cg.Nseq):
    for m in range(cg.Nreps) :
        tr_file = open("Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat",'r')
        topo_file = open("Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')


        #oxdna distances, types and angles
        fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun_33, tun_55 = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, tr_file, topo_file)

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
            types_unbn_33.append(tun_33)
            types_unbn_55.append(tun_55)
        else:
            for z in range(len(fr)): fene_r[l].append(fr[z])
            for z in range(len(sr)): stck_r[l].sppend(sr[z])
            for z in range(len(t4bn)): th4_bn[l].append(t4bn[z])
            for z in range(len(t5)): th5[l].append(t5[z])
            for z in range(len(t6)): th6[l].append(t6[z])
            for z in range(len(cp1)): cosphi1[l].append(cp1[z])
            for z in range(len(cp2)): cosphi2[l].append(cp2[z])
            for z in range(len(tbn)): types_bn[l].append(tbn[z])
            for z in range(len(hr)): hydr_r[l].append(hr[z])
            for z in range(len(t1)): th1[l].append(t1[z])
            for z in range(len(t2)): th2[l].append(t2[z])
            for z in range(len(t3)): th3[l].append(t3[z])
            for z in range(len(t4un)): th4_unbn[l].append(t4un[z])
            for z in range(len(t7)): th7[l].append(t7[z])
            for z in range(len(t8)): th8[l].append(t8[z])
            for z in range(len(tun_33)): types_unbn_33[l].append(tun_33[z])
            for z in range(len(tun_55)): types_unbn_55[l].append(tun_55[z])

        tr_file.close()
        topo_file.close()

        tr_file = open("Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat",'r')
        topo_file = open("Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

        #internal coordinates
        traj = read_oxdna_trajectory_standard_order(tr_file, topo_file)

        if l == 0 :
            fun.store_internal_coord(traj,l,cg.ids,cg.in_j[l],cg.fin_j[l],cg.in_snap,True)
            fun.store_normals(traj,l,cg.in_j[l],cg.fin_j[l],cg.in_snap,True)
        else  :
            fun.store_internal_coord(traj,l,cg.ids,cg.in_j[l],cg.fin_j[l],cg.in_snap,False)
            fun.store_normals(traj,l,cg.in_j[l],cg.fin_j[l],cg.in_snap,False)

        tr_file.close()
        topo_file.close()
"""
print("Flattening propeller")
print(cfun.internal_coords[0][0])
print(cfun.internal_coords[0][10])
print(cfun.target_mu[0])
if 1 in cg.ids_gs :
    fun.make_it_flat_GS([1])
print(cfun.internal_coords[0][0])
print(cfun.internal_coords[0][10])
print(cfun.target_mu[0])
print(cfun.target_mu[1])
"""

#Plot sampled coordinates
mus, covs = fun.ave_and_cov_sampled()
for l in range(cg.Nseq):
    fun.plot_gs_sampled(mus[l],l)
    fun.plot_std_sampled(covs[l],l)

cfun.save_mu = fun.Sequence_ave_GS(mus)
print("Average GS")
print(cfun.save_mu)

costb,cosot = fun.plengths_angles()

#make unbnd tensor square. Extra unbnd pairs have zero interaction energy.
max_ints = 0
for l in range(cg.Nseq) :
    for j in range(len(types_unbn_33[l])):
        if len(types_unbn_33[l][j]) > max_ints:
           max_ints = len(types_unbn_33[l][j])
print("max unbn pairs: "+str(max_ints))
for l in range(cg.Nseq) :
    for j in range(len(types_unbn_33[l])):
        for z in range(len(types_unbn_33[l][j]), max_ints):
            types_unbn_33[l][j].append(0)
            types_unbn_55[l][j].append(0)
            hydr_r[l][j].append(0.)

            th1[l][j].append(0.)
            th2[l][j].append(0.)
            th3[l][j].append(0.)

            th4_unbn[l][j].append(0.)
            th7[l][j].append(0.)
            th8[l][j].append(0.)

#Test: plot unbn_types

bases = ['A', 'C', 'G', 'T']

def to_letters(types):
    string = ""
    for l in range(len(types)):
        TY = types[l]
        ty0 = TY%4
        ty1 = (TY//4)%4
        ty2 = (TY//4//4)%4
        ty3 = (TY//4//4//4)%4
        ty = bases[ty0]+bases[ty1]+bases[ty2]+bases[ty3]
        string+=ty+" "

    return string

test_file = open("unbn_types_33.txt",'w')

for l in range(cg.Nseq) :
    print("SEQ"+str(l),file=test_file)
    print("0",file=test_file)
    print(to_letters(types_unbn_33[l][0]),file=test_file)
    print("10",file=test_file)
    print(to_letters(types_unbn_33[l][10]),file=test_file)
    print("20",file=test_file)
    print(to_letters(types_unbn_33[l][20]),file=test_file)

test_file.close()


test_file = open("unbn_types_55.txt",'w')

for l in range(cg.Nseq) :
    print("SEQ"+str(l),file=test_file)
    print("0",file=test_file)
    print(to_letters(types_unbn_55[l][0]),file=test_file)
    print("10",file=test_file)
    print(to_letters(types_unbn_55[l][10]),file=test_file)
    print("20",file=test_file)
    print(to_letters(types_unbn_55[l][20]),file=test_file)

test_file.close()
###################################################################################################
############## SETUP TENSORS FOR COST FUNCTION ####################################################
###################################################################################################

#create all tensors on the gpu. Change this to easily swap between gpu and cpu
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

if torch.cuda.is_available() :
    print("Running optimisation on gpu.")
else:
    print("gpu not available. Running optimisation on cpu.")


cfun.init_tensors(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                  th4_unbn, th7, th8, types_unbn_33, types_unbn_55, costb, cosot, shifts, OXPS_zero)

#build tensors imposing continuity
cfun.build_continuity_tensors()

#build masks (for selecting optim parameters and marginalising mu and cov) and symm tensors (for imposing symmetries)
cfun.build_masks_and_symm_tensors()

print("Continuity tensors created")

#print(cfun.LP_RED_IDS)
#print(cfun.LP_JK_IDS)
#print(cfun.LP_SUM_IDS_ROWS)
#print(cfun.LP_SUM_IDS_COLS)

#create reduced targets and compute cov^{-1}
cfun.reduce_targets()

print("Targets reduced")


#print(cfun.IDS_AVE_COV_SUM)
#print(cfun.IDS_AVE_COV_EXPAND)
#print(cfun.AVE_COV_RED_TARGET_COV)


#compute initial energy and modulation factors
cfun.compute_initial_energy()

print("Computed initial energy")


###################################################################################################
############## OPTIMISE ###########################################################################
###################################################################################################

#read initial optim parameters values
ids = []
for i in range(len(cfun.OPT_PAR_LIST)) :
    id = cfun.OPT_PAR_LIST[i][0]
    ty = cfun.OPT_PAR_LIST[i][1]
    ids.append(id*256+ty)


IDS_OP = torch.tensor(ids,device=device)
OPTI_PAR = torch.gather(torch.reshape(cfun.CURR_PARS,(-1,)),0,IDS_OP) #tensor with initial values of the opti parameters

print("OPTIMISING")

#We use scipy to handle the optimisation algorithm
#While we compute the cost function on the gpu with pytorch
#All the tensors necessary to the computation of the cost function are on the gpu

cfun.PAR0 = torch.clone(cfun.CURR_PARS) #clone the initial values of the parameters to the gpu
TMP = torch.tensor(OPTI_PAR,device='cpu') #copy opti parameters to the cpu.

#X0 is a numpy array with the optimised parameters; it is initialised and handled by scypi.optim on the cpu
#every time the cost function is computed, it is cloned to the device (gpu)
X0 = TMP.numpy()

#We impose bondaries on the optimisation parameters, to avoid pushing the reweighting too much
#The boundaries depend on the specific parameters (e.g. FENE_R0 can vary by +-3% max)
low_bond = torch.tensor(TMP, device='cpu').numpy()
up_bond = torch.tensor(TMP, device='cpu').numpy()

for n in range(len(low_bond)) :
    if cfun.OPT_PAR_LIST[n][0] == 1:
        lb = low_bond[n]*0.98
        ub = up_bond[n]*1.02
        if lb > 0.66: low_bond[n] = lb
        else: low_bond[n] = 0.66
        if ub > 0.82: up_bond[n] = ub
        else: up_bond[n] = 0.82
    elif cfun.OPT_PAR_LIST[n][0] == 2:
        low_bond[n] = low_bond[n]*0.98
        up_bond[n] = up_bond[n]*1.02
    elif cfun.OPT_PAR_LIST[n][0] == 45:
        lb = low_bond[n]*0.98
        ub = up_bond[n]*1.02
        if lb > 0.33: low_bond[n] = lb
        else: low_bond[n] = 0.33
        if ub > 0.42: up_bond[n] = ub
        else: up_bond[n] = 0.42
    elif cfun.OPT_PAR_LIST[n][0] == 101 or cfun.OPT_PAR_LIST[n][0] == 140:
        low_bond[n] = 2.96709
        up_bond[n] = 3.31609
    elif cfun.OPT_PAR_LIST[n][0] == 55:

        lb = low_bond[n]*0.5
        ub = up_bond[n]*2
        if lb > 0.6 : low_bond[n] = lb
        else : low_bond[n] = 0.6
        if ub < 3.0 : up_bond[n] = ub
        else : up_bond[n] = 4.0

    elif  cfun.OPT_PAR_LIST[n][0] == 60:

        lb = low_bond[n]*0.5
        ub = up_bond[n]*2
        if lb > 0.5 : low_bond[n] = lb
        else : low_bond[n] = 0.5
        if ub < 3.0 : up_bond[n] = ub
        else : up_bond[n] = 3.0

    else:
        low_bond[n] = low_bond[n]*0.5
        up_bond[n] = up_bond[n]*2.

bnd = optimize.Bounds(low_bond,up_bond)

#Compute average delta. Average value of FENE_DELTA is kept fixed, so that lt is kept at ~220 nm
cfun.AVE_DELTA = torch.mean(cfun.PAR0[2])
print("Ave delta:", cfun.AVE_DELTA)

print("S0: "+str(cfun.COST(X0)))
print("lb0: ")
print(cfun.LB)
print("lt0: ")
print(cfun.LT)

#callback is called at the end of each optimisation step. It prints the ratio x/x0
# and the current value of the cost function. x = current solution, x0 = initial parameters
def Callback(sol):
    print("x/x0: ")
    tmp = []
    for i in range(len(sol)):
        tmp.append(sol[i]/X0[i])
    print(tmp)
    print("S: "+str(cfun.COST(sol)))
    print("lb: ")
    print(cfun.LB)
    print("lt: ")
    print(cfun.LT)
    print("Ave delta fixed at: ")
    print(cfun.CURR_AVE_DELTA)


time_opti_start = time.time()

####### THIS LINE RUNS THE OPTIMISAION #######################
sol = optimize.minimize(cfun.COST,X0, method='L-BFGS-B', callback=Callback, bounds=bnd, options={'maxfun':10000,'iprint': 1})

S = cfun.COST(sol.x)
print("Final value of the cost function: "+str(S))


#printing runtime
timefile = open("runtime.txt",'w')
print("Total run time [s]: " + str(time.time() - start_time),file=timefile)
print("Optimisation run time [s]: " + str(time.time() - time_opti_start),file=timefile)
timefile.close()

#printing final SD file
in_SD_par_file = open("oxDNA_sequence_dependent_parameters_in.txt",'r')
fun.print_final_pfile(sol.x,in_SD_par_file)
#fun.print_final_pfile(OPTI_PAR,in_SD_par_file)
in_SD_par_file.close()

print("DONE")

