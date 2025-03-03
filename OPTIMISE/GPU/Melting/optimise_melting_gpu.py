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


cfun.convert_Ts_to_ox_units()

print("Converting temepratures to oxdna units")

print("n5 sim Ts:")
print(cfun.sim_Ts_n5)

print("n5 rew Ts:")
print(cfun.rew_Ts_n5)

test_file = open("opti_p_test.txt",'w')
for l in range(len(cfun.OPT_PAR_LIST)) :
    print(cfun.OPT_PAR_LIST[l],file=test_file)

test_file.close()

test_file = open("opti_p_test1.txt",'w')
print(OXPS_zero[45],file = test_file)
print(OXPS_zero[1],file = test_file)

print(OXPS_zero[34],file = test_file)
test_file.close()


#create all tensors on the gpu. Change this to easily swap between gpu and cpu
device = torch.device("cpu")

print("Memory usage after read optim options:")
print_memory_usage()

#########################################################################################################################
############## READ TRAJECTORY, COMPUTE OXDNA COORDINATES (i.e angles and distances) AND INTERNAL COORDINATES ###########
#########################################################################################################################

#################
### nbps = 5 ####
#################

for l in range(cg.Nseq_n5):

    en_off_1 = []
    hbs_s_1 = []

    for m in range(cg.Nreps) :
        split_en_file = open("n5/Seq"+str(l)+"/Rep"+str(m)+"/split_energy.dat", 'r')
        en_file = open("n5/Seq"+str(l)+"/Rep"+str(m)+"/energy.dat", 'r')

        nline = 0
        for line in split_en_file.readlines():
            if nline % 5 == 0:
                if nline / 5 <= cg.in_snap :
                    nline += 1
                    continue
                vals = line.strip().split()
                if int(vals[0]) % 100000 != 0 :
                    print("Something weird when reading en_offset")
                off = float(vals[2]) + float(vals[4]) + float(vals[7]) + float(vals[8])  #excl + coaxial + Debye
                en_off_1.append(off)
                #print(vals[0],off)
            nline += 1

        nline = 0
        for line in en_file.readlines():
            if nline % 10 == 0:
                if nline / 10 <= cg.in_snap :
                    nline += 1
                    continue
                vals = line.strip().split()
                if int(vals[0]) % 100000 != 0 :
                    print(vals[0])
                    print("Something weird when reading hbs_sampled")
                hbs_s_1.append(int(vals[5]))  #hbs
                #print(vals[0],vals[5])
            nline += 1

        split_en_file.close()
        en_file.close()

    cfun.en_offset_n5.append(en_off_1)
    cfun.hbs_sampled_n5.append(hbs_s_1)


if cg.read_coords_from_file :
    cfun.init_tensors_from_file_n5(device, shifts, OXPS_zero)

else :
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

    for l in range(cg.Nseq_n5):

        for m in range(cg.Nreps) :
            tr_file = open("n5/Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat",'r')
            topo_file = open("n5/Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

            #oxdna distances, types and angles
            fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun33, tun55 = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, tr_file, topo_file, cg.boxes_n5[l])

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
                types_unbn_33.append(tun33)
                types_unbn_55.append(tun55)
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
                types_unbn_33[l].extend(tun33)
                types_unbn_55[l].extend(tun55)

            tr_file.close()
            topo_file.close()

    #make unbnd tensor square. Extra unbnd pairs have zero interaction energy.

    max_ints = 0
    print("Len unbn")
    for l in range(cg.Nseq_n5) :
        for j in range(len(types_unbn_33[l])):
            if len(types_unbn_33[l][j]) > max_ints:
               max_ints = len(types_unbn_33[l][j])
    print("max unbn pairs: "+str(max_ints))
    for l in range(cg.Nseq_n5) :
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


    print("Check lengths:")
    print("fene_r: "+str(len(fene_r))+", "+str(len(fene_r[0]))+", "+ str(len(fene_r[0][0])))
    print("hydr_r: "+str(len(hydr_r))+", "+str(len(hydr_r[0]))+", "+ str(len(hydr_r[0][0])))


    print("Memory usage after reading n5 data:")
    print_memory_usage()


    cfun.init_tensors_n5(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                      th4_unbn, th7, th8, types_unbn_33, types_unbn_55, shifts, OXPS_zero)

    print("Memory usage after initialising n5 tensors:")
    print_memory_usage()

    del fene_r
    del stck_r
    del th4_bn
    del th5
    del th6
    del cosphi1
    del cosphi2
    del types_bn
    del hydr_r
    del th1
    del th2
    del th3
    del th4_unbn
    del th7
    del th8
    del types_unbn_33
    del types_unbn_55

    print("Memory usage after deleting n5 lists:")
    print_memory_usage()


#################
### nbps = 8 ####
#################

for l in range(cg.Nseq_n8):

    en_off_1 = []
    hbs_s_1 = []

    for m in range(cg.Nreps) :
        split_en_file = open("n8/Seq"+str(l)+"/Rep"+str(m)+"/split_energy.dat", 'r')
        en_file = open("n8/Seq"+str(l)+"/Rep"+str(m)+"/energy.dat", 'r')

        nline = 0
        for line in split_en_file.readlines():
            if nline % 5 == 0:
                if nline / 5 <= cg.in_snap :
                    nline += 1
                    continue
                vals = line.strip().split()
                if int(vals[0]) % 100000 != 0 :
                    print("Something weird when reading en_offset")
                off = float(vals[2]) + float(vals[4]) + float(vals[7]) + float(vals[8])  #excl + coaxial + Debye
                en_off_1.append(off)
                #print(vals[0],off)
            nline += 1

        nline = 0
        for line in en_file.readlines():
            if nline % 10 == 0:
                if nline / 10 <= cg.in_snap :
                    nline += 1
                    continue
                vals = line.strip().split()
                if int(vals[0]) % 100000 != 0 :
                    print(vals[0])
                    print("Something weird when reading hbs_sampled")
                hbs_s_1.append(int(vals[5]))  #hbs
                #print(vals[0],vals[5])
            nline += 1

        split_en_file.close()
        en_file.close()

    cfun.en_offset_n8.append(en_off_1)
    cfun.hbs_sampled_n8.append(hbs_s_1)


if cg.read_coords_from_file :
    cfun.init_tensors_from_file_n8(device, shifts)

else :

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

    for l in range(cg.Nseq_n8):

        for m in range(cg.Nreps) :

            tr_file = open("n8/Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat",'r')
            topo_file = open("n8/Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

            #oxdna distances, types and angles
            fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun33, tun55 = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, tr_file, topo_file, cg.boxes_n8[l])

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
                types_unbn_33.append(tun33)
                types_unbn_55.append(tun55)
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
                types_unbn_33[l].extend(tun33)
                types_unbn_55[l].extend(tun55)

            tr_file.close()
            topo_file.close()

    #make unbnd tensor square. Extra unbnd pairs have zero interaction energy.

    max_ints = 0
    for l in range(cg.Nseq_n8) :
        for j in range(len(types_unbn_33[l])):
            if len(types_unbn_33[l][j]) > max_ints:
               max_ints = len(types_unbn_33[l][j])
    print("max unbn pairs: "+str(max_ints))
    for l in range(cg.Nseq_n8) :
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


    print("Check lengths:")
    print("fene_r: "+str(len(fene_r))+", "+str(len(fene_r[0]))+", "+ str(len(fene_r[0][0])))
    print("hydr_r: "+str(len(hydr_r))+", "+str(len(hydr_r[0]))+", "+ str(len(hydr_r[0][0])))

    cfun.init_tensors_n8(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                      th4_unbn, th7, th8, types_unbn_33, types_unbn_55, shifts)


    del fene_r
    del stck_r
    del th4_bn
    del th5
    del th6
    del cosphi1
    del cosphi2
    del types_bn
    del hydr_r
    del th1
    del th2
    del th3
    del th4_unbn
    del th7
    del th8
    del types_unbn_33
    del types_unbn_55


"""

##################
### nbps = 15 ####
##################


for l in range(cg.Nseq_n15):

    en_off_1 = []
    hbs_s_1 = []

    for m in range(cg.Nreps) :
        tr_file = open("n15/Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat",'r')
        topo_file = open("n15/Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

        nline = 0
        for line in split_en_file.readlines():
            if nline % 5 == 0:
                if nline / 5 <= cg.in_snap :
                    nline += 1
                    continue
                vals = line.strip().split()
                if int(vals[0]) % 100000 != 0 :
                    print("Something weird when reading en_offset")
                off = float(vals[2]) + float(vals[4]) + float(vals[7]) + float(vals[8])  #excl + coaxial + Debye
                en_off_1.append(off)
                #print(vals[0],off)
            nline += 1

        nline = 0
        for line in en_file.readlines():
            if nline % 10 == 0:
                if nline / 10 <= cg.in_snap :
                    nline += 1
                    continue
                vals = line.strip().split()
                if int(vals[0]) % 100000 != 0 :
                    print(vals[0])
                    print("Something weird when reading hbs_sampled")
                hbs_s_1.append(int(vals[5]))  #hbs
                #print(vals[0],vals[5])
            nline += 1

        split_en_file.close()
        en_file.close()

    cfun.en_offset_n15.append(en_off_1)
    cfun.hbs_sampled_n15.append(hbs_s_1)



if cg.read_coords_from_file :
    cfun.init_tensors_from_file_n15(device, shifts)

else :

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

    for l in range(cg.Nseq_n15):

        for m in range(cg.Nreps) :
            tr_file = open("n15/Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat",'r')
            topo_file = open("n15/Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

            #oxdna distances, types and angles
            fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun33, tun55 = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, tr_file, topo_file, cg.boxes_n15[l])

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
                types_unbn_33.append(tun33)
                types_unbn_55.append(tun55)
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
                types_unbn_33[l].extend(tun33)
                types_unbn_55[l].extend(tun55)

            tr_file.close()
            topo_file.close()

    #make unbnd tensor square. Extra unbnd pairs have zero interaction energy.

    max_ints = 0
    for l in range(cg.Nseq_n15) :
        for j in range(len(types_unbn_33[l])):
            if len(types_unbn_33[l][j]) > max_ints:
               max_ints = len(types_unbn_33[l][j])
    print("max unbn pairs: "+str(max_ints))
    for l in range(cg.Nseq_n15) :
        for j in range(len(types_unbn_33[l])):
            for z in range(len(types_unbn_33[l][j]), max_ints):
                types_unbn_33[l][j].append(0)
                types_unbn_33[l][j].append(0)
                hydr_r[l][j].append(0.)

                th1[l][j].append(0.)
                th2[l][j].append(0.)
                th3[l][j].append(0.)

                th4_unbn[l][j].append(0.)
                th7[l][j].append(0.)
                th8[l][j].append(0.)


    print("Check lengths:")
    print("fene_r: "+str(len(fene_r))+", "+str(len(fene_r[0]))+", "+ str(len(fene_r[0][0])))
    print("hydr_r: "+str(len(hydr_r))+", "+str(len(hydr_r[0]))+", "+ str(len(hydr_r[0][0])))

    cfun.init_tensors_n15(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                      th4_unbn, th7, th8, types_unbn_33, types_unbn_55, shifts)

    del fene_r
    del stck_r
    del th4_bn
    del th5
    del th6
    del cosphi1
    del cosphi2
    del types_bn
    del hydr_r
    del th1
    del th2
    del th3
    del th4_unbn
    del th7
    del th8
    del types_unbn_33
    del types_unbn_55



"""


print("Memory usage after reading all data:")
print_memory_usage()


###################################################################################################
############## SETUP TENSORS FOR COST FUNCTION ####################################################
###################################################################################################


if cg.print_coords_to_file:
    ofile = open("dists_and_angles_n5.txt", 'w')
    cfun.print_dists_and_angles_n5(ofile)
    ofile.close()
    cfun.print_dists_and_angles_n8()

ofile = open("melting_tensors_n5.txt", 'w')
cfun.print_melting_things_n5(ofile)
ofile.close()

#build tensors imposing continuity
cfun.build_continuity_tensors()

print("Continuity tensors cool.")

#build masks (for selecting optim parameters and marginalising mu and cov) and symm tensors (for imposing symmetries)
cfun.build_symm_tensors()

print("Masks cool")

#create reduced targets and compute cov^{-1}
#cfun.reduce_targets()

#print("Targets reduced")

print("Memory usage after building masks all data:")
print_memory_usage()

#compute initial energy and modulation factors

cfun.compute_energy_n5()
cfun.compute_energy_n8()
#cfun.compute_initial_energy_n15()

print("Memory usage after computing in energy:")
print_memory_usage()

print("Computed initial energy")

ofile = open("energy_in_all_n5.txt", 'w')
ofile_ave = open("energy_in_ave_n5.txt", 'w')

cfun.print_energy_n5(ofile,ofile_ave)

ofile.close()
ofile_ave.close()

ofile = open("energy_in_all_n8.txt", 'w')
ofile_ave = open("energy_in_ave_n8.txt", 'w')

cfun.print_energy_n8(ofile,ofile_ave)

ofile.close()
ofile_ave.close()

print("Memory usage before optim-1:")
print_memory_usage()

#read initial optim parameters values

ids = []
for i in range(len(cfun.OPT_PAR_LIST)) :
    id = cfun.OPT_PAR_LIST[i][0]
    ty = cfun.OPT_PAR_LIST[i][1]
    ids.append(id*256+ty)


IDS_OP = torch.tensor(ids,device=device)
OPTI_PAR = torch.gather(torch.reshape(cfun.CURR_PARS,(-1,)),0,IDS_OP)


time_cfun = time.time()

list_params = []

print("OPTIMISING")

cfun.PAR0 = torch.clone(cfun.CURR_PARS)

TMP = torch.tensor(OPTI_PAR,device='cpu')

X0 = TMP.numpy()

low_bond = torch.tensor(TMP, device='cpu').numpy()
up_bond = torch.tensor(TMP, device='cpu').numpy()

for n in range(len(low_bond)) :
    if cfun.OPT_PAR_LIST[n][0] == 4:
        lb = low_bond[n]*0.5
        ub = up_bond[n]*1.5
        if lb > 0.65: low_bond[n] = lb
        else: low_bond[n] = 0.65
        if ub > 1.4: up_bond[n] = ub
        else: up_bond[n] = 1.4
    elif cfun.OPT_PAR_LIST[n][0] == 44:
        lb = low_bond[n]*0.5
        ub = up_bond[n]*1.5
        if lb > 1.0: low_bond[n] = lb
        else: low_bond[n] = 1.0
        if ub > 2.0: up_bond[n] = ub
        else: up_bond[n] = 2.0
    elif cfun.OPT_PAR_LIST[n][0] == 77 or cfun.OPT_PAR_LIST[n][0] == 116:
        low_bond[n] = 0
        up_bond[n] = 76
    else:
        low_bond[n] = low_bond[n]*0.5
        up_bond[n] = up_bond[n]*2.

bnd = optimize.Bounds(low_bond,up_bond)


print("Memory usage before optim-2:")
print_memory_usage()

print("S0: "+str(cfun.COST(X0)))


print("Target mTs n5:")
print(cfun.target_Tms_n5)

print("Initial mTs n5:")
print(cfun.current_mT_n5)

def Callback(sol):

    cg.Niter += 1
    cg.Diter_Trange += 1

    cg.update_rews = True

    print("x: ")
    tmp = []
    for i in range(len(sol)):
        tmp.append(sol[i]/X0[i])
    print(tmp)
    print("S: "+str(cfun.COST(sol)))

    print("Target mTs n5:")
    print(cfun.target_Tms_n5)

    print("Current mTs n5")
    print(cfun.current_mT_n5)

    print("Memory usage at end of iteration:")
    print_memory_usage()

time_opti_start = time.time()

print("Memory usage before optim-3:")
print_memory_usage()

####### THIS LINE RUNS THE OPTIMISAION #######################
sol = optimize.minimize(cfun.COST,X0, method='L-BFGS-B', callback=Callback, bounds=bnd, options={'maxiter':7,'iprint': 1})

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

