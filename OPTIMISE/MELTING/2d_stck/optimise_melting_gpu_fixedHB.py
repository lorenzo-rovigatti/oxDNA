
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

torch.set_printoptions(precision=6)
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

#compute debye huckel lambda before reading trajectories: we need to know this to correctly figure out the cutoff
cfun.compute_debye_huckel_lambda()

print("lambdas")
print(cfun.dh_l_n5)
print(cfun.dh_l_n8)
print(cfun.dh_l_n15)

test_file = open("opti_p_test.txt",'w')
for l in range(len(cfun.OPT_PAR_LIST)) :
    print(cfun.OPT_PAR_LIST[l],file=test_file)

test_file.close()

#test_file = open("opti_p_test1.txt",'w')
#print(OXPS_zero[45],file = test_file)
#print(OXPS_zero[1],file = test_file)

#print(OXPS_zero[34],file = test_file)
#test_file.close()


#create all tensors on the gpu. Change this to easily swap between gpu and cpu
#device = torch.device("cpu")
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

print("Memory usage after read optim options:")
print_memory_usage()

#########################################################################################################################
############## READ TRAJECTORY, COMPUTE OXDNA COORDINATES (i.e angles and distances) AND INTERNAL COORDINATES ###########
#########################################################################################################################

nevery_en_n5 = int(cg.delta_time_n5/cg.delta_print_en_n5)
nevery_split_n5 = int(cg.delta_time_n5/cg.delta_split_en_n5)
nevery_en_n8 = int(cg.delta_time_n8/cg.delta_print_en_n8)
nevery_split_n8 = int(cg.delta_time_n8/cg.delta_split_en_n8)
nevery_en_n15 = int(cg.delta_time_n15/cg.delta_print_en_n15)
nevery_split_n15 = int(cg.delta_time_n15/cg.delta_split_en_n15)

if nevery_en_n5 == 0 or nevery_split_n5 == 0:
    if nevery_en_n5 == 0:
        print("Energy was printed less frequently than snapshots were sampled (n5).")
        print("Cannot read order parameter value for all snapshots (n5)")
        print("Aborting.")
    if nevery_split_n5 == 0:
        print("Split_energy was printed less frequently than snapshots were sampled (n5).")
        print("Cannot read constant energy terms (e.g. coaxial) for all snapshots (n5)")
        print("Aborting.")

    exit(1)

if nevery_en_n8 == 0 or nevery_split_n8 == 0:
    if nevery_en_n8 == 0:
        print("Energy was printed less frequently than snapshots were sampled (n8).")
        print("Cannot read order parameter value for all snapshots (n8)")
        print("Aborting.")
    if nevery_split_n8 == 0:
        print("Split_energy was printed less frequently than snapshots were sampled (n8).")
        print("Cannot read constant energy terms (e.g. coaxial) for all snapshots (n8)")
        print("Aborting.")

    exit(1)

if nevery_en_n15 == 0 or nevery_split_n15 == 0:
    if nevery_en_n15 == 0:
        print("Energy was printed less frequently than snapshots were sampled (n15).")
        print("Cannot read order parameter value for all snapshots (n15)")
        print("Aborting.")
    if nevery_split_n15 == 0:
        print("Split_energy was printed less frequently than snapshots were sampled (n15).")
        print("Cannot read constant energy terms (e.g. coaxial) for all snapshots (n15)")
        print("Aborting.")

    exit(1)



#################
### nbps = 5 ####
#################


bases = ['A', 'C', 'G', 'T', 'E']

def type_to_base4(TY) :
    ty0 = TY%5
    ty1 = (TY//5)%5
    ty2 = (TY//5//5)%5
    ty3 = (TY//5//5//5)%5

    return str(ty0)+str(ty1)+str(ty2)+str(ty3)


for l in range(cg.Nseq_n5):

    en_off_1 = []
    hbs_s_1 = []

    N_pts = 1
    if cg.parallel_tempering : N_pts = cg.N_PT_reps_n5

    for rp in range(N_pts) :    #first we loop over pt replicas, then repetitions: we are stacking configuraions of same PT temperature together
        for m in range(cg.Nreps) :

            split_en_file_name = "n5/Seq"+str(l)+"/Rep"+str(m)+"/split_energy.dat"
            en_file_name = "n5/Seq"+str(l)+"/Rep"+str(m)+"/energy.dat"

            if cg.parallel_tempering :

                split_en_file_name = "n5/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_split_energy.dat"
                en_file_name = "n5/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_energy.dat"

            split_en_file = open(split_en_file_name, 'r')
            en_file = open(en_file_name, 'r')

            nline = 0
            for line in split_en_file.readlines():
                if nline % nevery_split_n5 == 0:
                    if nline / nevery_split_n5 <= cg.in_snap :
                        nline += 1
                        continue
                    vals = line.strip().split()
                    if int(vals[0]) % cg.delta_time_n5 != 0 :
                        print("Something weird when reading en_offset. Snap time is off")
                    off = float(vals[2]) + float(vals[4]) + float(vals[7])
                    if cg.debye_huckel == False: off += float(vals[8])  #excl (2,4) + coaxial (7) + Debye (8)
                    en_off_1.append(off*10) #off is per nucleotinde!
                    #print(vals[0],off)
                nline += 1

            nline = 0
            for line in en_file.readlines():
                if nline % nevery_en_n5 == 0:
                    if nline / nevery_en_n5 <= cg.in_snap :
                        nline += 1
                        continue
                    vals = line.strip().split()
                    if int(vals[0]) % cg.delta_time_n5 != 0 :
                        print("Something weird when reading hbs_sampled. Snap time is off")
                    hbs_s_1.append(int(vals[5]))  #hbs
                    #print(vals[0],vals[5])
                nline += 1

            split_en_file.close()
            en_file.close()

    cfun.en_offset_n5.append(en_off_1)
    cfun.hbs_sampled_n5.append(hbs_s_1)

    ofile_hsb=open("hsb_n5.dat", 'w')
    print(cfun.hbs_sampled_n5, file=ofile_hsb)
    ofile_hsb.close()

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
    debye_huckel_r = []
    debye_huckel_types = []
    debye_huckel_charge_cut = []

    rclow , rchigh, rcut_dh_n5, rcut_dh_n8, rcut_dh_n15 = fun.find_cuts_for_lists(OXPS_zero)
    print("cuts: "+str(rclow)+" "+str(rchigh))
    print("DH cuts n5: ", rcut_dh_n5)
    print("DH cuts n8: ", rcut_dh_n8)
    print("DH cuts n15: ", rcut_dh_n15)

    for l in range(cg.Nseq_n5):

        print("Reading seq " + str(l) + " n5")

        N_pts = 1
        if cg.parallel_tempering : N_pts = cg.N_PT_reps_n5

        for rp in range(N_pts) :
            for m in range(cg.Nreps) :

                tr_file_name = "n5/Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat"
                topo_file = open("n5/Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

                if cg.parallel_tempering :
                    tr_file_name = "n5/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_trajectory.dat"

                tr_file = open(tr_file_name, 'r')

                _print = False
                if l == 0 and m == 0 and rp == 0:
                    _print = False

                #oxdna distances, types and angles
                fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun33, tun55, dh_r, dh_ty, dh_chcut = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, rcut_dh_n5[l], tr_file, topo_file, cg.boxes_n5[l], _print)

                if m == 0 and rp == 0:
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
                    debye_huckel_r.append(dh_r)
                    debye_huckel_types.append(dh_ty)
                    debye_huckel_charge_cut.append(dh_chcut)
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
                    debye_huckel_r[l].extend(dh_r)
                    debye_huckel_types[l].extend(dh_ty)
                    debye_huckel_charge_cut[l].extend(dh_chcut)


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

    max_ints = 0

    print("Len debye huckle")
    for l in range(cg.Nseq_n5) :
        for j in range(len(debye_huckel_types[l])):
            if len(debye_huckel_types[l][j]) > max_ints:
               max_ints = len(debye_huckel_types[l][j])
    print("max debye huckle pairs: "+str(max_ints))
    for l in range(cg.Nseq_n5) :
        for j in range(len(debye_huckel_types[l])):
            for z in range(len(debye_huckel_types[l][j]), max_ints):
                debye_huckel_types[l][j].append(0)
                debye_huckel_r[l][j].append(100.)
                debye_huckel_charge_cut[l][j].append(1)

    print("Check lengths:")
    if len(fene_r) > 0 : print("fene_r: "+str(len(fene_r))+", "+str(len(fene_r[0]))+", "+ str(len(fene_r[0][0])))
    if len(hydr_r) > 0 : print("hydr_r: "+str(len(hydr_r))+", "+str(len(hydr_r[0]))+", "+ str(len(hydr_r[0][0])))
    if len(dh_r) > 0 : print("debye_huckel_r: "+str(len(debye_huckel_r))+", "+str(len(debye_huckel_r[0]))+", "+ str(len(debye_huckel_r[0][0])))


    print("Memory usage after reading n5 data:")
    print_memory_usage()


    ofile = open("test_fener.dat", 'w')
    print(fene_r,file=ofile)
    ofile.close()


    cfun.init_tensors_n5(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                      th4_unbn, th7, th8, types_unbn_33, types_unbn_55, debye_huckel_r, debye_huckel_types, debye_huckel_charge_cut, shifts, OXPS_zero)

    print("Memory usage after initialising n5 tensors:")
    print_memory_usage()

    #for l in types_bn[0][0] :
    #    print(type_to_base4(l))

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

    N_pts = 1
    if cg.parallel_tempering : N_pts = cg.N_PT_reps_n8

    for rp in range(N_pts) :
        for m in range(cg.Nreps) :

            split_en_file_name = "n8/Seq"+str(l)+"/Rep"+str(m)+"/split_energy.dat"
            en_file_name = "n8/Seq"+str(l)+"/Rep"+str(m)+"/energy.dat"

            if cg.parallel_tempering :

                split_en_file_name = "n8/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_split_energy.dat"
                en_file_name = "n8/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_energy.dat"

            split_en_file = open(split_en_file_name, 'r')
            en_file = open(en_file_name, 'r')


            nline = 0
            for line in split_en_file.readlines():
                if nline % nevery_split_n8 == 0:
                    if nline / nevery_split_n8 <= cg.in_snap :
                        nline += 1
                        continue
                    vals = line.strip().split()
                    if int(vals[0]) % cg.delta_time_n8 != 0 :
                        print("Something weird when reading en_offset. Snap time is off")
                    off = float(vals[2]) + float(vals[4]) + float(vals[7])
                    if cg.debye_huckel == False: off += float(vals[8])  #excl + coaxial + Debye
                    en_off_1.append(off*16) #off is per nucleotinde!
                    #print(vals[0],off)
                nline += 1

            nline = 0
            for line in en_file.readlines():
                if nline % nevery_en_n8 == 0:
                    if nline / nevery_en_n8 <= cg.in_snap :
                        nline += 1
                        continue
                    vals = line.strip().split()
                    if int(vals[0]) % cg.delta_time_n8 != 0 :
                        print("Something weird when reading hbs_sampled. Snap time is off")
                    hbs_s_1.append(int(vals[5]))  #hbs
                    #print(vals[0],vals[5])
                nline += 1

            split_en_file.close()
            en_file.close()

    cfun.en_offset_n8.append(en_off_1)
    cfun.hbs_sampled_n8.append(hbs_s_1)



if cg.read_coords_from_file :
    cfun.init_tensors_from_file_n8(device)

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
    debye_huckel_r = []
    debye_huckel_types = []
    debye_huckel_charge_cut = []

    for l in range(cg.Nseq_n8):

        print("Reading seq " + str(l) + " n8")

        N_pts = 1
        if cg.parallel_tempering : N_pts = cg.N_PT_reps_n8

        for rp in range(N_pts) :
            for m in range(cg.Nreps) :

                tr_file_name = "n8/Seq"+str(l)+"/Rep"+str(m)+"/trajectory.dat"
                topo_file = open("n8/Seq"+str(l)+"/Rep"+str(m)+"/generated.top", 'r')

                if cg.parallel_tempering :
                    tr_file_name = "n8/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_trajectory.dat"

                tr_file = open(tr_file_name, 'r')

                fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun33, tun55, dh_r, dh_ty, dh_chcut = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, rcut_dh_n8[l], tr_file, topo_file, cg.boxes_n8[l])

                if m == 0 and rp == 0:
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
                    debye_huckel_r.append(dh_r)
                    debye_huckel_types.append(dh_ty)
                    debye_huckel_charge_cut.append(dh_chcut)
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
                    debye_huckel_r[l].extend(dh_r)
                    debye_huckel_types[l].extend(dh_ty)
                    debye_huckel_charge_cut[l].extend(dh_chcut)

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


    max_ints = 0
    print("Len debye huckel")
    for l in range(cg.Nseq_n8) :
        for j in range(len(debye_huckel_types[l])):
            if len(debye_huckel_types[l][j]) > max_ints:
               max_ints = len(debye_huckel_types[l][j])
    print("max debye huckel pairs: "+str(max_ints))
    for l in range(cg.Nseq_n8) :
        for j in range(len(debye_huckel_types[l])):
            for z in range(len(debye_huckel_types[l][j]), max_ints):
                debye_huckel_types[l][j].append(0)
                debye_huckel_r[l][j].append(100.)
                debye_huckel_charge_cut[l][j].append(1)

    print("Check lengths:")
    if len(fene_r) > 0 : print("fene_r: "+str(len(fene_r))+", "+str(len(fene_r[0]))+", "+ str(len(fene_r[0][0])))
    if len(hydr_r) > 0 : print("hydr_r: "+str(len(hydr_r))+", "+str(len(hydr_r[0]))+", "+ str(len(hydr_r[0][0])))

    cfun.init_tensors_n8(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                      th4_unbn, th7, th8, types_unbn_33, types_unbn_55, debye_huckel_r, debye_huckel_types, debye_huckel_charge_cut)


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


##################
### nbps = 15 ####
##################

for l in range(cg.Nseq_n15):

    en_off_1 = []
    hbs_s_1 = []

    N_pts = 1
    if cg.parallel_tempering : N_pts = cg.N_PT_reps_n15

    for rp in range(N_pts) :
        for m in range(cg.Nreps) :

            split_en_file_name = "n15/Seq"+str(l)+"/Rep"+str(m)+"/split_energy.dat"
            en_file_name = "n15/Seq"+str(l)+"/Rep"+str(m)+"/energy.dat"

            if cg.parallel_tempering :

                split_en_file_name = "n15/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_split_energy.dat"
                en_file_name = "n15/Seq"+str(l)+"/Rep"+str(m)+"/mpi_"+str(rp)+"_energy.dat"

            split_en_file = open(split_en_file_name, 'r')
            en_file = open(en_file_name, 'r')

            nline = 0
            for line in split_en_file.readlines():
                if nline % nevery_split_n15 == 0:
                    if nline / nevery_split_n15 <= cg.in_snap :
                        nline += 1
                        continue
                    vals = line.strip().split()
                    if int(vals[0]) % cg.delta_time_n15 != 0 :
                        print("Something weird when reading en_offset. Snap time is off")
                    off = float(vals[2]) + float(vals[4]) + float(vals[7])
                    if cg.debye_huckel == False: off += float(vals[8])  #excl + coaxial + Debye
                    en_off_1.append(off*30) #off is per nucleotinde!
                    #print(vals[0],off)
                nline += 1

            nline = 0
            for line in en_file.readlines():
                if nline % nevery_en_n15 == 0:
                    if nline / nevery_en_n15 <= cg.in_snap :
                        nline += 1
                        continue
                    vals = line.strip().split()
                    if int(vals[0]) % cg.delta_time_n15 != 0 :
                        print("Something weird when reading hbs_sampled. Snap time is off")
                    hbs_s_1.append(int(vals[5]))  #hbs
                    #print(vals[0],vals[5])
                nline += 1

            split_en_file.close()
            en_file.close()

    cfun.en_offset_n15.append(en_off_1)
    cfun.hbs_sampled_n15.append(hbs_s_1)


if cg.read_coords_from_file :
    cfun.init_tensors_from_file_n15(device)

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
    debye_huckel_r = []
    debye_huckel_types = []
    debye_huckel_charge_cut = []

    for l in range(cg.Nseq_n15):

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

                fr, sr, t4bn, t5, t6, cp1, cp2, tbn, hr, t1, t2, t3, t4un, t7, t8, tun33, tun55, dh_r, dh_ty, dh_chcut = fun.read_oxdna_trajectory_dist_and_angles(rclow, rchigh, rcut_dh_n15[l], tr_file, topo_file, cg.boxes_n15[l])

                if m == 0 and rp == 0:
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
                    debye_huckel_r.append(dh_r)
                    debye_huckel_types.append(dh_ty)
                    debye_huckel_charge_cut.append(dh_chcut)
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
                    debye_huckel_r[l].extend(dh_r)
                    debye_huckel_types[l].extend(dh_ty)
                    debye_huckel_charge_cut[l].extend(dh_chcut)

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
                types_unbn_55[l][j].append(0)
                hydr_r[l][j].append(0.)

                th1[l][j].append(0.)
                th2[l][j].append(0.)
                th3[l][j].append(0.)

                th4_unbn[l][j].append(0.)
                th7[l][j].append(0.)
                th8[l][j].append(0.)

    max_ints = 0
    print("Len debye huckel")
    for l in range(cg.Nseq_n15) :
        for j in range(len(debye_huckel_types[l])):
            if len(debye_huckel_types[l][j]) > max_ints:
               max_ints = len(debye_huckel_types[l][j])
    print("max debye huckel pairs: "+str(max_ints))
    for l in range(cg.Nseq_n15) :
        for j in range(len(debye_huckel_types[l])):
            for z in range(len(debye_huckel_types[l][j]), max_ints):
                debye_huckel_types[l][j].append(0)
                debye_huckel_r[l][j].append(100.)
                debye_huckel_charge_cut[l][j].append(1)


    for l in range(cg.Nseq_n15) :
        for j in range(len(types_unbn_33[l])):
            if j < 40: print(len(types_unbn_33[l][j]))

    print("Check lengths:")
    if len(fene_r) > 0 : print("fene_r: "+str(len(fene_r))+", "+str(len(fene_r[0]))+", "+ str(len(fene_r[0][0])))
    if len(hydr_r) > 0 : print("hydr_r: "+str(len(hydr_r))+", "+str(len(hydr_r[0]))+", "+ str(len(hydr_r[0][0])))

    cfun.init_tensors_n15(device,fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                      th4_unbn, th7, th8, types_unbn_33, types_unbn_55, debye_huckel_r, debye_huckel_types, debye_huckel_charge_cut)

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


print("Memory usage after reading all data:")
print_memory_usage()


###################################################################################################
############## SETUP TENSORS FOR COST FUNCTION ####################################################
###################################################################################################

#extend lambda tensor for deby huckel
if cg.debye_huckel:

   cfun.extend_lambda_for_debye_huckel()
   print("Extended lamda. n5 size:")
   print(cfun.DH_LAMBDA_EXT_n5.shape)

#extend sim_Ts tensors if using parallel tempering sim_Ts[seq][pt_replica]->sim_Ts[seq][conf]
if cg.parallel_tempering: cfun.extend_Ts_and_weights_for_PT()

if cg.print_coords_to_file:
    ofile = open("dists_and_angles_n5.txt", 'w')
    cfun.print_dists_and_angles_n5(ofile)
    ofile.close()
    cfun.print_dists_and_angles_n8()
    cfun.print_dists_and_angles_n15()


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

if cg.Nseq_n5 > 0 : cfun.compute_energy_n5()
if cg.Nseq_n8 > 0 : cfun.compute_energy_n8()
if cg.Nseq_n15 > 0 : cfun.compute_initial_energy_n15()

print("Memory usage after computing in energy:")
print_memory_usage()

print("Computed initial energy")

ofile = open("energy_in_all_n5.txt", 'w')
ofile_ave = open("energy_in_ave_n5.txt", 'w')

print("# of sequences n5", cg.Nseq_n5)

if cg.Nseq_n5 > 0 : cfun.print_energy_n5(ofile,ofile_ave)

ofile.close()
ofile_ave.close()

#print("Exit +debug.")
#exit(1)

if cg.Nseq_n8 > 0 : cfun.print_energy_n8()
if cg.Nseq_n15 > 0 : cfun.print_energy_n15()

print("Memory usage before optim-1:")
print_memory_usage()

#read initial optim parameters values

ids = []
for i in range(len(cfun.OPT_PAR_LIST)) :
    id = cfun.OPT_PAR_LIST[i][0]
    ty = cfun.OPT_PAR_LIST[i][1]
    ids.append(id*625+ty)

#now that we computed the energy, let's take the T factor in eps stack away, so that in pars are correct.

cfun.CURR_PARS[par_index[44]] = cfun.CURR_PARS[par_index[44]]/ (1.0 - cg.stck_fact_eps + (0.1 * 9.0 * cg.stck_fact_eps))


IDS_OP = torch.tensor(ids,device=device)
OPTI_PAR = torch.gather(torch.reshape(cfun.CURR_PARS,(-1,)),0,IDS_OP)


torch.set_printoptions(profile="full")

print("STACKING STRENGTH")
print(cfun.CURR_PARS[44])

torch.set_printoptions(profile="default")

time_cfun = time.time()

list_params = []

print("OPTIMISING")

cfun.PAR0 = torch.clone(cfun.CURR_PARS)


for i in range(5):
   for j in range(5):
       #AT
       ty = i+j*125+0*5+3*25
       cfun.PAR0[4][ty] =
       ty = i+j*125+3*5+0*25
       cfun.PAR0[4][ty] =

       #CG
       ty = i+j*125+1*5+2*25
       cfun.PAR0[4][ty] =
       ty = i+j*125+2*5+1*25
       cfun.PAR0[4][ty] =

TMP = torch.tensor(OPTI_PAR,device='cpu')

X0 = TMP.numpy()

low_bond = torch.tensor(TMP, device='cpu').numpy()
up_bond = torch.tensor(TMP, device='cpu').numpy()

for n in range(len(low_bond)) :
    if cfun.OPT_PAR_LIST[n][0] == 4:
        lb = low_bond[n]*0.5
        ub = up_bond[n]*1.5
        if lb > 0.60: low_bond[n] = lb
        else: low_bond[n] = 0.60
        if ub < 1.4: up_bond[n] = ub
        else: up_bond[n] = 1.4
    elif cfun.OPT_PAR_LIST[n][0] == 44:
        lb = low_bond[n]*0.25
        ub = up_bond[n]*2.0
        low_bond[n] = lb
        up_bond[n] = ub
        if lb > 1.2: low_bond[n] = lb
        else: low_bond[n] = 1.3
        if ub < 1.95: up_bond[n] = ub
        else: up_bond[n] = 1.95
    elif cfun.OPT_PAR_LIST[n][0] == 77 or cfun.OPT_PAR_LIST[n][0] == 116:
        low_bond[n] = 0
        up_bond[n] = 76.1
    else:
        low_bond[n] = low_bond[n]*0.5
        up_bond[n] = up_bond[n]*2.

bnd = optimize.Bounds(low_bond,up_bond)


print("Memory usage before optim-2:")
print_memory_usage()

print("S0: "+str(cfun.COST(X0)))

torch.set_printoptions(profile="full")

if cg.Nseq_n5 > 0:
    print("Initial mTs n5:")
    print(cfun.current_mT_n5)

    print("Target mTs n5:")
    print(cfun.TARGET_mTs_n5)

    print("Delta n5:")
    print(cfun.current_mT_n5-cfun.TARGET_mTs_n5)

if cg.Nseq_n8 > 0:
    print("Initial mTs n8:")
    print(cfun.current_mT_n8)

    print("Target mTs n8:")
    print(cfun.TARGET_mTs_n8)

    print("Delta n8:")
    print(cfun.current_mT_n8-cfun.TARGET_mTs_n8)

if cg.Nseq_n15 > 0:
    print("Initial mTs n15:")
    print(cfun.current_mT_n15)

    print("Target mTs n15:")
    print(cfun.TARGET_mTs_n15)

    print("Delta n15:")
    print(cfun.current_mT_n15-cfun.TARGET_mTs_n15)


torch.set_printoptions(profile="default")


def Callback(sol):

    cg.Niter += 1
    cg.Diter_Trange += 1

    cg.update_rews = True

    print("x_norm: ")
    tmp = []
    for i in range(len(sol)):
        tmp.append(sol[i]/X0[i])
    print(tmp)

    print("x: ")
    tmp = []
    for i in range(len(sol)):
        tmp.append(sol[i])
    print(tmp)


    print("S: "+str(cfun.COST(sol)))

    torch.set_printoptions(profile="full")

    if cg.Nseq_n5 > 0:
        print("mTs n5 at iteration "+str(cg.Niter)+":")
        print(cg.good_n5)
        print(cfun.current_mT_n5)

        print("Target mTs n5:")
        print(cfun.TARGET_mTs_n5)

        print("Delta n5:")
        print(cfun.current_mT_n5-cfun.TARGET_mTs_n5)

    if cg.Nseq_n8 > 0:
        print("mTs n8 at iteration "+str(cg.Niter)+":")
        print(cg.good_n8)
        print(cfun.current_mT_n8)

        print("Target mTs n8:")
        print(cfun.TARGET_mTs_n8)

        print("Delta n8:")
        print(cfun.current_mT_n8-cfun.TARGET_mTs_n8)


    if cg.Nseq_n15 > 0:
        print("mTs n15 at iteration "+str(cg.Niter)+":")
        print(cg.good_n15)
        print(cfun.current_mT_n15)

        print("Target mTs n15:")
        print(cfun.TARGET_mTs_n15)

        print("Delta n15:")
        print(cfun.current_mT_n15-cfun.TARGET_mTs_n15)

    torch.set_printoptions(profile="default")

    print("Memory usage at end of iteration:")
    print_memory_usage()

time_opti_start = time.time()

print("Memory usage before optim-3:")
print_memory_usage()

####### THIS LINE RUNS THE OPTIMISAION #######################
sol = optimize.minimize(cfun.COST,X0, method='L-BFGS-B', callback=Callback, bounds=bnd, options={'maxiter':50,'iprint': 1})


#change slightly parameters at boundaries to avoid problems with next iteration
par_fin = sol.x

for n in range(len(par_fin)) :
    if (par_fin[n]-low_bond[n])/(low_bond[n]+0.00000001) <= 0.0001: par_fin[n]= par_fin[n]*1.0001
    if (up_bond[n]-par_fin[n])/up_bond[n] <= 0.0001: par_fin[n] = par_fin[n]*0.9999


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
in_SD_par_file.close()
in_SD_par_file = open("oxDNA_sequence_dependent_parameters_in.txt",'r')
fun.print_final_pfile_AllFromOpt(sol.x,in_SD_par_file)
in_SD_par_file.close()

print("DONE")
