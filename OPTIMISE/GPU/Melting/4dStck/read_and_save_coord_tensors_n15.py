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
from scipy import optimize
import resource

import torch.multiprocessing as mp

#from mpi4py import MPI


#1 thread per sequence
#mpi_size = cg.comm.Get_size()
#mpi_rank = cg.comm.Get_rank()

#mpi_seq_id = mpi_rank

#mpi_comm_world = MPI.COMM_WORLD


def print_memory_usage():
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"Memory usage: {memory / 1024:.2f} MB")


# READ CONFIG FILE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 optimise.py config_file")
    sys.exit()

start_time = time.time()

config_file = sys.argv[1]

PARS_LIST = parl.PARS_LIST
par_index = parl.par_index

###XXXTODO ANDREA: check this out
cg.first_step_flag = True

#X = torch.tensor([0,1,2])
#print(X.device)

#torch.set_printoptions(precision=6)
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


OXPS_zero, shifts = fun.init_oxpars(pars_from_modelh, vals_from_modelh, over_indices, over_vals, 0.1, stck_fact_eps)

###################################################################################################
############## READ OPTIM OPTIONS  ################################################################
###################################################################################################

#READ PARAMETERS
if fun.read_config(config_file) == False :
    sys.exit()


cg.print_energy_to_file = False #True if we want to print energy to file
cg.print_coords_to_file = True #True if we want to print coordinates to file
cg.read_energy_from_file = False #True if we want to read initial energy from file. Useful because par0 might not be the sampled parameters (e.g. in geometry-melting cycle)
cg.read_coords_from_file = False

cfun.convert_Ts_to_ox_units()

print("Converting temepratures to oxdna units")

print("n5 sim Ts:")
print(cfun.sim_Ts_n5)

print("n5 rew Ts:")
print(cfun.rew_Ts_n5)

#compute debye huckel lambda before reading trajectories: we need to know this to correctly figure out the cutoff
cfun.compute_debye_huckel_lambda()
"""
print("lambdas")
print(cfun.dh_l_n5)
print(cfun.dh_l_n8)
print(cfun.dh_l_n15)


test_file = open("opti_p_test.txt",'w')
for l in range(len(cfun.OPT_PAR_LIST)) :
    print(cfun.OPT_PAR_LIST[l],file=test_file)

test_file.close()
"""
#test_file = open("opti_p_test1.txt",'w')
#print(OXPS_zero[45],file = test_file)
#print(OXPS_zero[1],file = test_file)

#print(OXPS_zero[34],file = test_file)
#test_file.close()


#create all tensors on the gpu. Change this to easily swap between gpu and cpu
device = torch.device("cpu")

print("Memory usage after read optim options:")
print_memory_usage()

#########################################################################################################################
############## READ TRAJECTORY, COMPUTE OXDNA COORDINATES (i.e angles and distances) AND INTERNAL COORDINATES ###########
#########################################################################################################################

#################
### nbps = 15 ####
#################


bases = ['A', 'C', 'G', 'T', 'E']

def type_to_base4(TY) :
    ty0 = TY%5
    ty1 = (TY//5)%5
    ty2 = (TY//5//5)%5
    ty3 = (TY//5//5//5)%5

    return str(ty0)+str(ty1)+str(ty2)+str(ty3)

def read_n15_seq(id) :

    l = int(id)

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
    debye_huckel_r = []
    debye_huckel_types = []
    debye_huckel_charge_cut = []


    rclow , rchigh, rcut_dh_n5, rcut_dh_n8, rcut_dh_n15 = fun.find_cuts_for_lists(OXPS_zero)
    print("cuts: "+str(rclow)+" "+str(rchigh))
    print("DH cuts n5: ", rcut_dh_n5)
    print("DH cuts n8: ", rcut_dh_n8)
    print("DH cuts n15: ", rcut_dh_n15)

    print("Reading seq " + str(l) + " n15")

    N_pts = 1
    if cg.parallel_tempering : N_pts = cg.N_PT_reps_n15

    for rp in range(N_pts) :
        for m in range(cg.Nreps) :

            tr_file_name = "n15/Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat"
            topo_file = open("n15/Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

            if cg.parallel_tempering :
                tr_file_name = "n15/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_trajectory.dat"

            tr_file = open(tr_file_name, 'r')

            #oxdna distances, types and angles
            fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun33, tun55, dh_r, dh_ty, dh_chcut = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, rcut_dh_n15[l], tr_file, topo_file, cg.boxes_n15[l])

            if m == 0 and rp == 0:
                fene_r = fr
                stck_r = sr
                th4_bn = t4bn
                th5 = t5
                th6 = t6
                cosphi1 = cp1
                cosphi2 = cp2
                types_bn = tbn
                hydr_r = hr
                th1 = t1
                th2 = t2
                th3 = t3
                th4_unbn = t4un
                th7 = t7
                th8 = t8
                types_unbn_33 = tun33
                types_unbn_55 = tun55
                debye_huckel_r = dh_r
                debye_huckel_types = dh_ty
                debye_huckel_charge_cut = dh_chcut
            else:
                fene_r.extend(fr)
                stck_r.extend(sr)
                th4_bn.extend(t4bn)
                th5.extend(t5)
                th6.extend(t6)
                cosphi1.extend(cp1)
                cosphi2.extend(cp2)
                types_bn.extend(tbn)
                hydr_r.extend(hr)
                th1.extend(t1)
                th2.extend(t2)
                th3.extend(t3)
                th4_unbn.extend(t4un)
                th7.extend(t7)
                th8.extend(t8)
                types_unbn_33.extend(tun33)
                types_unbn_55.extend(tun55)
                debye_huckel_r.extend(dh_r)
                debye_huckel_types.extend(dh_ty)
                debye_huckel_charge_cut.extend(dh_chcut)


            tr_file.close()
            topo_file.close()

    return fene_r, stck_r, th4_bn, th5, th6 , cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3, th4_unbn, th7, th8, types_unbn_33, types_unbn_55, debye_huckel_r, debye_huckel_types, debye_huckel_charge_cut


fene_r_all = []
stck_r_all = []
th4_bn_all = []
th5_all = []
th6_all = []
cosphi1_all = []
cosphi2_all = []
types_bn_all = []
hydr_r_all = []
th1_all = []
th2_all = []
th3_all = []
th4_unbn_all = []
th7_all = []
th8_all = []
types_unbn_33_all = []
types_unbn_55_all = []
debye_huckel_r_all = []
debye_huckel_types_all = []
debye_huckel_charge_cut_all = []


if __name__ == '__main__':

    seq_ids = torch.arange(cg.Nseq_n15)

    with mp.Pool(cg.Nseq_n15) as pool:

        results = pool.map(read_n15_seq, seq_ids)

        #print(len(results))

    for i in range(len(results)):
        fene_r_all.append(results[i][0])
        stck_r_all.append(results[i][1])
        th4_bn_all.append(results[i][2])
        th5_all.append(results[i][3])
        th6_all.append(results[i][4])
        cosphi1_all.append(results[i][5])
        cosphi2_all.append(results[i][6])
        types_bn_all.append(results[i][7])
        hydr_r_all.append(results[i][8])
        th1_all.append(results[i][9])
        th2_all.append(results[i][10])
        th3_all.append(results[i][11])
        th4_unbn_all.append(results[i][12])
        th7_all.append(results[i][13])
        th8_all.append(results[i][14])
        types_unbn_33_all.append(results[i][15])
        types_unbn_55_all.append(results[i][16])
        debye_huckel_r_all.append(results[i][17])
        debye_huckel_types_all.append(results[i][18])
        debye_huckel_charge_cut_all.append(results[i][19])


#make unbnd tensor square. Extra unbnd pairs have zero interaction energy.
max_ints = 0
print("Len unbn")
for l in range(cg.Nseq_n15) :
    for j in range(len(types_unbn_33_all[l])):
        if len(types_unbn_33_all[l][j]) > max_ints:
            max_ints = len(types_unbn_33_all[l][j])
print("max unbn pairs: "+str(max_ints))
for l in range(cg.Nseq_n15) :
    for j in range(len(types_unbn_33_all[l])):
        for z in range(len(types_unbn_33_all[l][j]), max_ints):
            types_unbn_33_all[l][j].append(0)
            types_unbn_55_all[l][j].append(0)
            hydr_r_all[l][j].append(0.)

            th1_all[l][j].append(0.)
            th2_all[l][j].append(0.)
            th3_all[l][j].append(0.)

            th4_unbn_all[l][j].append(0.)
            th7_all[l][j].append(0.)
            th8_all[l][j].append(0.)

max_ints = 0

print("Len debye huckle")
for l in range(cg.Nseq_n15) :
    for j in range(len(debye_huckel_types_all[l])):
        if len(debye_huckel_types_all[l][j]) > max_ints:
           max_ints = len(debye_huckel_types_all[l][j])
print("max debye huckle pairs: "+str(max_ints))
for l in range(cg.Nseq_n15) :
    for j in range(len(debye_huckel_types_all[l])):
        for z in range(len(debye_huckel_types_all[l][j]), max_ints):
            debye_huckel_types_all[l][j].append(0)
            debye_huckel_r_all[l][j].append(100.)
            debye_huckel_charge_cut_all[l][j].append(1)


for l in range(len(hydr_r_all)):
    for j in range(len(hydr_r_all[l])):
        if len(hydr_r_all[l][j]) != len(hydr_r_all[0][0]): print(l,j, "Len hydr_r is weird; ", len(hydr_r_all[l][j]))
print("Check lengths:")
if len(fene_r_all) > 0 : print("fene_r: "+str(len(fene_r_all))+", "+str(len(fene_r_all[0]))+", "+ str(len(fene_r_all[0][0])))
if len(hydr_r_all) > 0 : print("hydr_r: "+str(len(hydr_r_all))+", "+str(len(hydr_r_all[0]))+", "+ str(len(hydr_r_all[0][0])))
if len(debye_huckel_r_all) > 0 : print("debye_huckel_r: "+str(len(debye_huckel_r_all))+", "+str(len(debye_huckel_r_all[0]))+", "+ str(len(debye_huckel_r_all[0][0])))


print("Memory usage after reading n15 data:")
print_memory_usage()


ofile = open("test_fener.dat", 'w')
print(fene_r_all,file=ofile)
ofile.close()


cfun.init_tensors_n15(device,fene_r_all, stck_r_all, th4_bn_all, th5_all, th6_all, cosphi1_all, cosphi2_all, types_bn_all, hydr_r_all, th1_all, th2_all, th3_all,\
                  th4_unbn_all, th7_all, th8_all, types_unbn_33_all, types_unbn_55_all, debye_huckel_r_all, debye_huckel_types_all, debye_huckel_charge_cut_all)

print("Memory usage after initialising n15 tensors:")
print_memory_usage()

cfun.print_dists_and_angles_n15()

#for l in types_bn[0][0] :
#    print(type_to_base4(l))

del fene_r_all
del stck_r_all
del th4_bn_all
del th5_all
del th6_all
del cosphi1_all
del cosphi2_all
del types_bn_all
del hydr_r_all
del th1_all
del th2_all
del th3_all
del th4_unbn_all
del th7_all
del th8_all
del types_unbn_33_all
del types_unbn_55_all

print("Memory usage after deleting n15 lists:")
print_memory_usage()


print("DONE")
