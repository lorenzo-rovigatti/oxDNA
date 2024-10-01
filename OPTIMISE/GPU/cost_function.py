#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:00:47 2024

@author: yqb22156
"""

import parameters_list as parl
import torch
import functions as fun
import config as cg

#target_lb = 45.0 #target average long range bending persistence length (in nm).
#target_lt = 220.0 #target average long range torsional persistence length (in nm).

#target_C = 140 #target average long range C modulus
#target_Ar = 40 #target average long range Ar (A2) modulus

PARS_LIST = parl.PARS_LIST
par_index = parl.par_index

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

FENE_R = []
STCK_R = []
TH4_BN = []
TH5 = []
TH6 = []
COSPHI1 = []
COSPHI2 = []

TYPES_BN = []


HYDR_R = []
TH1 = []
TH2 = []
TH3 = []
TH4_UNBN = []
TH7 = []
TH8 = []

TYPES_UNBN = []

SHIFT_STCK = []
SHIFT_HYDR = []

PARS_IN = []
CURR_PARS = []
ZEROS_BN = []
ONES_BN = []
ZEROS_UNBN = []


EN_FENE_IN = []
EN_STCK_IN = []
STCK_MOD_IN = []
EN_HYDR_IN = []
HYDR_MOD_IN = []
EN_CRST_33_IN = []
CRST_33_MOD_IN = []
EN_CRST_55_IN = []
CRST_55_MOD_IN = []

UPDATE_MASK = []
SYMM_LIST = [] #list with parameters to symmetrise - dim 1 is par index, dim 2 type
SYMM_LIST_SYMM = [] #list of symmetric parameters - SYMM_LIST_SYMM[i][j] is the symmetric counterpart of SYMM_LIST[i][j]


internal_coords = []
INTERNAL_COORDS = []

target_mu = []
target_cov = []
TARGET_MU = []
TARGET_COV = []

COV_RED_COV_MASK = []
MU_RED_COV_MASK = []
MU_RED_MU_MASK = []

COV_RED_TARGET_COV = [] 
COV_RED_TARGET_COV_m1 = []
MU_RED_TARGET_MU = [] 
MU_RED_TARGET_COV = [] 
MU_RED_TARGET_COV_m1 = []

#initailise tensors with coordinates and 
def init_tensors(dev, fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                  th4_unbn, th7, th8, types_unbn, shifts, OXPS_zero) :
    
    global device

    global FENE_R
    global STCK_R
    global TH4_BN
    global TH5
    global TH6
    global COSPHI1
    global COSPHI2
    global TYPES_BN


    global HYDR_R
    global TH1
    global TH2
    global TH3
    global TH4_UNBN
    global  TH7
    global TH8
    global TYPES_UNBN

    global SHIFT_STCK
    global SHIFT_HYDR

    global PARS_IN
    global CURR_PARS
    global ZEROS_BN
    global ONES_BN
    global ZEROS_UNBN


    global EN_FENE_IN
    global EN_STCK_IN
    global STCK_MOD_IN
    global EN_HYDR_IN
    global HYDR_MOD_IN
    global EN_CRST_33_IN
    global CRST_33_MOD_IN
    global EN_CRST_55_IN
    global CRST_55_MOD_IN
    
    global UPDATE_MASK
    
    global internal_coords
    global INTERNAL_COORDS
    
    global target_mu
    global TARGET_MU
    
    global target_cov
    global TARGET_COV


    device = dev
    
    #initialise tensors
    FENE_R = torch.tensor(fene_r, device=device)
    STCK_R = torch.tensor(stck_r, device=device)
    TH4_BN = torch.tensor(th4_bn, device=device)
    TH5 = torch.tensor(th5, device=device)
    TH6 = torch.tensor(th6, device=device)
    COSPHI1 = torch.tensor(cosphi1, device=device)
    COSPHI2 = torch.tensor(cosphi2, device=device)

    TYPES_BN = torch.tensor(types_bn,device=device)


    HYDR_R = torch.tensor(hydr_r, device=device)
    TH1 = torch.tensor(th1, device=device)
    TH2 = torch.tensor(th2, device=device)
    TH3 = torch.tensor(th3, device=device)
    TH4_UNBN = torch.tensor(th4_unbn, device=device)
    TH7 = torch.tensor(th7, device=device)
    TH8 = torch.tensor(th8, device=device)

    TYPES_UNBN = torch.tensor(types_unbn,device=device)

    SHIFT_STCK = torch.tensor(shifts[1], device=device)
    SHIFT_HYDR = torch.tensor(shifts[0], device=device)

    PARS_IN = torch.tensor(OXPS_zero,device=device)
    CURR_PARS = torch.tensor(OXPS_zero,device=device)
    ZEROS_BN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    ONES_BN = torch.ones(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    ZEROS_UNBN = torch.zeros(len(types_unbn),len(types_unbn[0]),len(types_bn[0][0]),device=device)

    EN_FENE_IN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    EN_STCK_IN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    STCK_MOD_IN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    EN_HYDR_IN = torch.zeros(len(types_unbn),len(types_unbn[0]),len(types_bn[0][0]),device=device)
    HYDR_MOD_IN = torch.zeros(len(types_unbn),len(types_unbn[0]),len(types_bn[0][0]),device=device)
    EN_CRST_33_IN = torch.zeros(len(types_unbn),len(types_unbn[0]),len(types_bn[0][0]),device=device)
    CRST_33_MOD_IN = torch.zeros(len(types_unbn),len(types_unbn[0]),len(types_bn[0][0]),device=device)
    EN_CRST_55_IN = torch.zeros(len(types_unbn),len(types_unbn[0]),len(types_bn[0][0]),device=device)
    CRST_55_MOD_IN = torch.zeros(len(types_unbn),len(types_unbn[0]),len(types_bn[0][0]),device=device)
    
    UPDATE_MASK = torch.zeros(len(PARS_LIST),256,device=device) # 1 if parameter is optimised, 0 otherwise
    
    INTERNAL_COORDS = torch.tensor(internal_coords,device=device)
    
    TARGET_MU = torch.tensor(target_mu,device=device)
    INTERNAL_COORDS = torch.tensor(internal_coords,device=device)


def F1(R,EPS,R0,A,BLOW,BHIGH,RLOW,RHIGH,RCLOW,RCHIGH,SHIFT,ZERO) :
    f1_1 = EPS*torch.square( 1.-torch.exp(-(R-R0)*A) )-SHIFT
    f1_2 = EPS*BLOW*torch.square( (R-RCLOW) )  #blow
    f1_3 = EPS*BHIGH*torch.square( (R-RCHIGH) )   #bhigh
    f1_4 = ZERO
    
    #apply if conditions
    f1 = torch.where(R<RCHIGH,f1_3,f1_4)
    f1 = torch.where(R<=RHIGH,f1_1,f1)
    f1 = torch.where(R<RLOW,f1_2,f1)
    f1 = torch.where(R<RCLOW,f1_4,f1)
    
    return f1

def F2(R, K, R0, RC, BLOW, BHIGH, RLOW, RHIGH, RCLOW, RCHIGH, ZERO):
    f2_1 = K*0.5*(torch.square( R-R0 )-torch.square( RC-R0 ))
    f2_2 = K*BLOW*torch.square( (R-RCLOW) )  #blow
    f2_3 = K*BHIGH*torch.square( (R-RCHIGH) )   #bhigh
    f2_4 = ZERO
   
    f2 = torch.where(R<RCHIGH,f2_3,f2_4)
    f2 = torch.where(R<=RHIGH,f2_1,f2)
    f2 = torch.where(R<RLOW,f2_2,f2)
    f2 = torch.where(R<RCLOW,f2_4,f2)
    
    return f2


def F4(TH, TH0, A, B, TS, TC, ZERO) :
    f4_1 = 1.-A*torch.square( (TH-TH0) )
    f4_2 = B*torch.square( (TH0-TC-TH) )    #low
    f4_3 = B*torch.square( (TH0+TC-TH) )   #high
    f4_4 = ZERO
    
    f4 = torch.where(TH<TH0+TC,f4_3,f4_4)
    f4 = torch.where(TH<=TH0+TS,f4_1,f4)
    f4 = torch.where(TH<TH0-TS,f4_2,f4)
    f4 = torch.where(TH<TH0-TC,f4_4,f4)

    return f4

def F5(COSPHI, A, B, XC, XS, ZERO, ONE):
    #cosphi1, cosphi2 angular modulations
    f5_1 = 1.-A*torch.square( COSPHI )
    f5_2 = B*torch.square( (XC-COSPHI ) )    #low
    f5_3 = ONE   #high
    f5_4 = ZERO
    
    f5 = torch.where(COSPHI>=0,f5_3,f5_1)
    f5 = torch.where(COSPHI<XS,f5_2,f5)
    f5 = torch.where(COSPHI<XC,f5_4,f5)
    
    return f5


#given a configuration, computes the oxdna potential energy.
#only FENE, hydrogen, stacking and cross-stacking
#takes tensors with distances and angles for every interacting pair of every configurations of multiple trajectories and computes energy with parameterset = PARS


#TENSOR FOR TREATING PARAMETRS UPDATE


f1_R0_ID = []
f1_A_ID = []
f1_RC_ID = []
f1_BL_ID = []
f1_BH_ID = []
f1_RL_ID = []
f1_RH_ID = []
f1_RCL_ID = []
f1_RCH_ID = []

f4_A_ID = []
f4_B_ID = []
f4_TS_ID = []
f4_TC_ID = []

def build_continuity_tensors() :
    
    global f1_R0_ID
    global f1_A_ID
    global f1_RC_ID
    global f1_BL_ID
    global f1_BH_ID
    global f1_RL_ID
    global f1_RH_ID
    global f1_RCL_ID
    global f1_RCH_ID
    
    global f4_A_ID
    global f4_B_ID
    global f4_TS_ID
    global f4_TC_ID
    
    #f1
    f1_r0_id = [5,46]
    f1_a_id = [6,47]
    f1_rc_id = [7,48]
    f1_bl_id = [8,49]
    f1_bh_id = [9,50]
    f1_rl_id = [10,51]
    f1_rh_id = [11,52]
    f1_rcl_id = [12,53]
    f1_rch_id = [13,54]
    
    #f4
    f4_a_id = [15,20,25,30,35,40,56,61,66,88,93,98,103,108,113,127,132,137,142,147,152]
    f4_b_id = [16,21,26,31,36,41,57,62,67,89,94,99,104,109,114,128,133,138,143,148,153]
    f4_ts_id = [17,22,27,32,37,42,58,63,68,90,95,100,105,110,115,129,134,139,144,149,154]
    f4_tc_id = [18,23,28,33,38,43,59,64,69,91,96,101,106,111,116,130,135,140,145,150,155]

    f1_R0_ID = torch.tensor(f1_r0_id,device=device)
    f1_A_ID = torch.tensor(f1_a_id,device=device)
    f1_RC_ID = torch.tensor(f1_rc_id,device=device)
    f1_BL_ID = torch.tensor(f1_bl_id,device=device)
    f1_BH_ID = torch.tensor(f1_bh_id,device=device)
    f1_RL_ID = torch.tensor(f1_rl_id,device=device)
    f1_RH_ID = torch.tensor(f1_rh_id,device=device)
    f1_RCL_ID = torch.tensor(f1_rcl_id,device=device)
    f1_RCH_ID = torch.tensor(f1_rch_id,device=device)
    
    f4_A_ID = torch.tensor(f4_a_id,device=device)
    f4_B_ID = torch.tensor(f4_b_id,device=device)
    f4_TS_ID = torch.tensor(f4_ts_id,device=device)
    f4_TC_ID = torch.tensor(f4_tc_id,device=device)
    
    return 



OPT_PAR_LIST = [] # dim 1 is par id, dim2 type

def add_opti_par(name) :
    
    global OPT_PAR_LIST 
    
    vals_cn = name.split("_")
    if len(vals_cn) == 0:
        return
    if vals_cn[0] == "STCK" and len(vals_cn) == 3:
            ty=fun.base_to_id(vals_cn[1])*4+fun.base_to_id(vals_cn[2])*4*4
            OPT_PAR_LIST.append( [PARS_LIST.index("STCK_EPS"),ty])
    elif vals_cn[0] == "HYDR" and len(vals_cn) == 3:
                ty=fun.base_to_id(vals_cn[1])*4+fun.base_to_id(vals_cn[2])*4*4   #from base 4 to base 10
                OPT_PAR_LIST.append( [PARS_LIST.index("HYDR_EPS"),ty] )
    elif vals_cn[0] == "HYDR" or vals_cn[0] == "CRST" :
        par_name = vals_cn[0]
        for i in range(1,len(vals_cn)-2):
            par_name += "_"+vals_cn[i]
        ty=fun.base_to_id(vals_cn[len(vals_cn)-2])*4+fun.base_to_id(vals_cn[len(vals_cn)-1])*4*4  #from base 4 to base 10
        OPT_PAR_LIST.append( [PARS_LIST.index(par_name),ty] )
        
    elif vals_cn[0] == "STCK" or vals_cn[0] == "FENE" :
        par_name = vals_cn[0]
        for i in range(1,len(vals_cn)-4):
            par_name += "_"+vals_cn[i]
        ty=fun.base_to_id(vals_cn[len(vals_cn)-4])+fun.base_to_id(vals_cn[len(vals_cn)-3])*4+fun.base_to_id(vals_cn[len(vals_cn)-2])*4*4+fun.base_to_id(vals_cn[len(vals_cn)-1])*4*4*4    #from base 4 to base 10
        OPT_PAR_LIST.append( PARS_LIST.index(par_name),ty )
        
    return

#UPDATE_MASK is the the tensor which tells which parameters are optimised
#SYMM_LIST and SYMM_LIST_SYMM are used to impose symmetries. i.e. SYMM_LIST_SYMM[i,j] = SYMM_LIST [i,j]
def build_masks_and_symm_tensors() :
    
    global OPT_PAR_LIST
    global UPDATE_MASK
    global device
    global SYMM_LIST
    global SYMM_LIST_SYMM
    
    #reduction masks
    
    global COV_RED_COV_MASK
    global MU_RED_COV_MASK
    global MU_RED_MU_MASK
    
    l = 0
    
    red_n = (cg.fin_j[l]-cg.in_j[l]+1)*(len(cg.ids)-len(cg.ids_cov))
    
    
    COV_RED_COV_MASK = torch.tensor(cg.dimension[l]-red_n,cg.dimension[l]-red_n,2,device=device)
    
    nrow=-1

    for i in range(cg.dimension[l]):
        if cg.ids[i%len(cg.ids)]  in cg.ids_cov:
    
            nrow+=1
            ncol=-1
            for j in range(cg.dimension[l]):
                if cg.ids[j%len(cg.ids)]  in cg.ids_cov:
                    ncol+=1
                    COV_RED_COV_MASK[nrow,ncol] = [i,j]
    
    
    red_n = (cg.fin_j[l]-cg.in_j[l]+1)*(len(cg.ids)-len(cg.ids_mu))
    
    
    MU_RED_COV_MASK = torch.tensor(cg.dimension[l]-red_n,cg.dimension[l]-red_n,2,device=device)
    
    nrow=-1

    for i in range(cg.dimension[l]):
        if cg.ids[i%len(cg.ids)]  in cg.ids_mu:
    
            nrow+=1
            ncol=-1
            for j in range(cg.dimension[l]):
                if cg.ids[j%len(cg.ids)]  in cg.ids_mu:
                    ncol+=1
                    MU_RED_COV_MASK[nrow,ncol] = [i,j]
                    
                    
    MU_RED_MU_MASK = torch.tensor(cg.dimension[l]-red_n,cg.dimension[l]-red_n,device=device)
    
    nrow=-1

    for i in range(cg.dimension[l]):
        if cg.ids[i%len(cg.ids)]  in cg.ids_mu:
    
            nrow+=1

            MU_RED_MU_MASK[nrow] = i                    
    
    
    sl = []
    sls = []
    
    for i in range(len(OPT_PAR_LIST)) :
        ID = OPT_PAR_LIST[0]
        TY = OPT_PAR_LIST[1]
        UPDATE_MASK[ID,TY] = 1
        if ID in parl.compl_symm_ids_2d :
            ty3 = TY%4
            ty2 = (TY//4)%4
            ty1 = (TY//4//4)%4
            ty0 = (TY//4//4//4)%4
            for l in range(4) :
                for m in range(4) :
                    TY_S = m+(4-ty2)*4+(4-ty1)*4*4+l*4*4*4
                    sl.append(ID,TY)
                    sls.append(ID,TY_S)
        if ID in parl.compl_symm_ids_2d_no_TT_AA :
            ty3 = TY%4
            ty2 = (TY//4)%4
            ty1 = (TY//4//4)%4
            ty0 = (TY//4//4//4)%4
            if (ty1 == 0 and ty3 == 0) or (ty1 == 3 and ty3 == 3):
                continue
            for l in range(4) :
                for m in range(4) :
                    TY_S = m+(4-ty2)*4+(4-ty1)*4*4+l*4*4*4
                    sl.append(ID,TY)
                    sls.append(ID,TY_S)
                    
        if ID in parl.compl_symm_ids :
            sl.append(ID,TY)
            ty3 = TY%4
            ty2 = (TY//4)%4
            ty1 = (TY//4//4)%4
            ty0 = (TY//4//4//4)%4
            TY_S = (4-ty3)+(4-ty2)*4+(4-ty1)*4*4+(4-ty0)*4*4*4
            sls.append(ID,TY_S)
            if ID in parl.is_th2 :
                sl.append(ID,TY)
                sls.append(parl.is_th3[parl.is_th2.index[ID]],TY)
                sl.append(ID,TY)
                sls.append(parl.is_th3[parl.is_th2.index[ID]],TY_S)
            if ID in parl.is_th5 :
                sl.append(ID,TY)
                sls.append(parl.is_th6[parl.is_th5.index[ID]],TY)
                sl.append(ID,TY)
                sls.append(parl.is_th6[parl.is_th5.index[ID]],TY_S)
            if ID in parl.is_th7 :
                sl.append(ID,TY)
                sls.append(parl.is_th8[parl.is_th7.index[ID]],TY)
                sl.append(ID,TY)
                sls.append(parl.is_th8[parl.is_th7.index[ID]],TY_S)
                
        if ID in parl.perm_symm_ids_2d :
            ty3 = TY%4
            ty2 = (TY//4)%4
            ty1 = (TY//4//4)%4
            ty0 = (TY//4//4//4)%4
            for l in range(4) :
                for m in range(4) :
                    TY_S = m+(ty2)*4+(ty1)*4*4+l*4*4*4
                    sl.append(ID,TY)
                    sls.append(ID,TY_S)
                    if ID in parl.is_th2 :
                        sl.append(ID,TY)
                        sls.append(parl.is_th3[parl.is_th2.index[ID]],TY)
                        sl.append(ID,TY)
                        sls.append(parl.is_th3[parl.is_th2.index[ID]],TY_S)
                    if ID in parl.is_th5 :
                        sl.append(ID,TY)
                        sls.append(parl.is_th6[parl.is_th5.index[ID]],TY)
                        sl.append(ID,TY)
                        sls.append(parl.is_th6[parl.is_th5.index[ID]],TY_S)
                    if ID in parl.is_th7 :
                        sl.append(ID,TY)
                        sls.append(parl.is_th8[parl.is_th7.index[ID]],TY)
                        sl.append(ID,TY)
                        sls.append(parl.is_th8[parl.is_th7.index[ID]],TY_S)
        if ID in parl.perm_symm_ids :
            sl.append(ID,TY)
            ty3 = TY%4
            ty2 = (TY//4)%4
            ty1 = (TY//4//4)%4
            ty0 = (TY//4//4//4)%4
            TY_S = ty3+ty2*4+ty1*4*4+ty0*4*4*4
            sls.append(ID,TY_S)
            if ID in parl.is_th2 :
                sl.append(ID,TY)
                sls.append(parl.is_th3[parl.is_th2.index[ID]],TY)
                sl.append(ID,TY)
                sls.append(parl.is_th3[parl.is_th2.index[ID]],TY_S)
            if ID in parl.is_th5 :
                sl.append(ID,TY)
                sls.append(parl.is_th6[parl.is_th5.index[ID]],TY)
                sl.append(ID,TY)
                sls.append(parl.is_th6[parl.is_th5.index[ID]],TY_S)
            if ID in parl.is_th7 :
                sl.append(ID,TY)
                sls.append(parl.is_th8[parl.is_th7.index[ID]],TY)
                sl.append(ID,TY)
                sls.append(parl.is_th8[parl.is_th7.index[ID]],TY_S)
     
            
     
    SYMM_LIST = torch.tensor(sl,device=device)
    SYMM_LIST_SYMM = torch.tensor(sls,device=device)
    return


#PREPARE REDUCED TARGETS AND COV^{-1}

def reduce_targets() :
    
    global TARGET_MU
    global TARGET_COV
    global COV_RED_TARGET_COV
    global COV_RED_TARGET_COV_m1
    global MU_RED_TARGET_MU
    global MU_RED_TARGET_COV
    global MU_RED_TARGET_COV_m1
    
    COV_RED_TARGET_COV = TARGET_COV[COV_RED_COV_MASK]
    MU_RED_TARGET_COV = TARGET_COV[MU_RED_COV_MASK]
    MU_RED_TARGET_MU = TARGET_MU[MU_RED_MU_MASK]
    
    sx,sy,sz = COV_RED_TARGET_COV.size()
    COV_RED_TARGET_COV_m1 = torch.tensor(sx,sy,sz,device=device)
    
    sx,sy,sz = MU_RED_TARGET_COV.size()
    MU_RED_TARGET_COV_m1 = torch.tensor(sx,sy,sz,device=device)
    
    for l in range(cg.Nseq) :
        COV_RED_TARGET_COV_m1[l] = torch.linalg.inv( COV_RED_TARGET_COV )
        MU_RED_TARGET_COV_m1[l] = torch.linalg.inv( MU_RED_TARGET_COV )
    
    return

#COMPUTE INITIAL ENERGY AND MODULATIONS

def compute_initial_energy() :
    
    global EN_FENE_IN
    global EN_STCK_IN
    global STCK_MOD_IN
    global EN_HYDR_IN
    global HYDR_MOD_IN
    global EN_CRST_33_IN
    global CRST_33_MOD_IN
    global EN_CRST_55_IN
    global CRST_55_MOD_IN

    #FENE
    EN_FENE_IN = -PARS_IN[par_index[0]][TYPES_BN]/2.*torch.log( 1.-torch.square( FENE_R-PARS_IN[par_index[1]][TYPES_BN] )/PARS_IN[par_index[3]][TYPES_BN])
    
    #STACKING
    #radial part
    #piecewise energy: compute all possibilities, and filter with where
    
    f1 = F1(STCK_R, PARS_IN[par_index[44]][TYPES_BN], PARS_IN[par_index[45]][TYPES_BN], PARS_IN[par_index[46]][TYPES_BN], PARS_IN[par_index[47]][TYPES_BN],PARS_IN[par_index[48]][TYPES_BN],\
            PARS_IN[par_index[49]][TYPES_BN], PARS_IN[par_index[50]][TYPES_BN], PARS_IN[par_index[51]][TYPES_BN], PARS_IN[par_index[52]][TYPES_BN], PARS_IN[par_index[53]][TYPES_BN], SHIFT_STCK[TYPES_BN], ZEROS_BN)
    
    #th4, th5, th6 Angular modulations
    
    
    f4_th4 = F4(TH4_BN,PARS_IN[par_index[54]][TYPES_BN],PARS_IN[par_index[55]][TYPES_BN],PARS_IN[par_index[56]][TYPES_BN],\
                PARS_IN[par_index[57]][TYPES_BN],PARS_IN[par_index[58]][TYPES_BN],ZEROS_BN)
    
    f4_th5 = F4(TH5,PARS_IN[par_index[59]][TYPES_BN],PARS_IN[par_index[60]][TYPES_BN],PARS_IN[par_index[61]][TYPES_BN],\
                PARS_IN[par_index[62]][TYPES_BN],PARS_IN[par_index[63]][TYPES_BN],ZEROS_BN)
    
    f4_th6 = F4(TH6,PARS_IN[par_index[64]][TYPES_BN],PARS_IN[par_index[65]][TYPES_BN],PARS_IN[par_index[66]][TYPES_BN],\
                PARS_IN[par_index[67]][TYPES_BN],PARS_IN[par_index[68]][TYPES_BN],ZEROS_BN)
    
    #cosphi1, cosphi2 angular modulations
    
    f5_phi1 = F5(COSPHI1,PARS_IN[par_index[69]][TYPES_BN],PARS_IN[par_index[70]][TYPES_BN],PARS_IN[par_index[71]][TYPES_BN],PARS_IN[par_index[72]][TYPES_BN],ZEROS_BN,ONES_BN)
    
    f5_phi2 = F5(COSPHI2,PARS_IN[par_index[73]][TYPES_BN],PARS_IN[par_index[74]][TYPES_BN],PARS_IN[par_index[75]][TYPES_BN],PARS_IN[par_index[76]][TYPES_BN],ZEROS_BN,ONES_BN)
    
    EN_STCK_IN = f1*f4_th4*f4_th5*f4_th6*f5_phi1*f5_phi2
    STCK_MOD_IN = f1*f4_th4*f4_th5*f4_th6 #this is the product of the modulations we are optimising!
     
     
    #HYDROGEN
    
    f1 = F1(HYDR_R, PARS_IN[par_index[4]][TYPES_UNBN], PARS_IN[par_index[5]][TYPES_UNBN], PARS_IN[par_index[6]][TYPES_UNBN], PARS_IN[par_index[7]][TYPES_UNBN],PARS_IN[par_index[8]][TYPES_UNBN],\
            PARS_IN[par_index[9]][TYPES_UNBN], PARS_IN[par_index[10]][TYPES_UNBN], PARS_IN[par_index[11]][TYPES_UNBN], PARS_IN[par_index[12]][TYPES_UNBN], PARS_IN[par_index[13]][TYPES_UNBN], SHIFT_HYDR[TYPES_UNBN], ZEROS_UNBN)
    
    f4_th1 = F4(TH1,PARS_IN[par_index[14]][TYPES_UNBN],PARS_IN[par_index[15]][TYPES_UNBN],PARS_IN[par_index[16]][TYPES_UNBN],\
                PARS_IN[par_index[17]][TYPES_UNBN],PARS_IN[par_index[18]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th2 = F4(TH2,PARS_IN[par_index[19]][TYPES_UNBN],PARS_IN[par_index[20]][TYPES_UNBN],PARS_IN[par_index[21]][TYPES_UNBN],\
                PARS_IN[par_index[22]][TYPES_UNBN],PARS_IN[par_index[23]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th3 = F4(TH3,PARS_IN[par_index[24]][TYPES_UNBN],PARS_IN[par_index[25]][TYPES_UNBN],PARS_IN[par_index[26]][TYPES_UNBN],\
                PARS_IN[par_index[27]][TYPES_UNBN],PARS_IN[par_index[28]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th4 = F4(TH4_UNBN,PARS_IN[par_index[29]][TYPES_UNBN],PARS_IN[par_index[30]][TYPES_UNBN],PARS_IN[par_index[31]][TYPES_UNBN],\
                PARS_IN[par_index[32]][TYPES_UNBN],PARS_IN[par_index[33]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th7 = F4(TH7,PARS_IN[par_index[34]][TYPES_UNBN],PARS_IN[par_index[35]][TYPES_UNBN],PARS_IN[par_index[36]][TYPES_UNBN],\
                PARS_IN[par_index[37]][TYPES_UNBN],PARS_IN[par_index[38]][TYPES_UNBN],ZEROS_UNBN)
    
    f4_th8 = F4(TH8,PARS_IN[par_index[39]][TYPES_UNBN],PARS_IN[par_index[40]][TYPES_UNBN],PARS_IN[par_index[41]][TYPES_UNBN],\
                PARS_IN[par_index[42]][TYPES_UNBN],PARS_IN[par_index[43]][TYPES_UNBN],ZEROS_UNBN)
    
    EN_HYDR_IN = f1*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    HYDR_MOD_IN  = f4_th4
    
    #CROSS STACKING
    #33
    
    f2 = F2(PARS_IN[par_index[77]][TYPES_UNBN], PARS_IN[par_index[78]][TYPES_UNBN], PARS_IN[par_index[79]][TYPES_UNBN], PARS_IN[par_index[80]][TYPES_UNBN], PARS_IN[par_index[81]][TYPES_UNBN],\
            PARS_IN[par_index[82]][TYPES_UNBN], PARS_IN[par_index[83]][TYPES_UNBN], PARS_IN[par_index[84]][TYPES_UNBN], PARS_IN[par_index[85]][TYPES_UNBN], ZEROS_UNBN)
    
    f4_th1 = F4(TH1,PARS_IN[par_index[86]][TYPES_UNBN],PARS_IN[par_index[87]][TYPES_UNBN],PARS_IN[par_index[88]][TYPES_UNBN],\
                PARS_IN[par_index[89]][TYPES_UNBN],PARS_IN[par_index[90]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th2 = F4(TH2,PARS_IN[par_index[91]][TYPES_UNBN],PARS_IN[par_index[92]][TYPES_UNBN],PARS_IN[par_index[93]][TYPES_UNBN],\
                PARS_IN[par_index[94]][TYPES_UNBN],PARS_IN[par_index[95]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th3 = F4(TH3,PARS_IN[par_index[96]][TYPES_UNBN],PARS_IN[par_index[97]][TYPES_UNBN],PARS_IN[par_index[98]][TYPES_UNBN],\
                PARS_IN[par_index[99]][TYPES_UNBN],PARS_IN[par_index[100]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th4 = F4(TH4_UNBN,PARS_IN[par_index[101]][TYPES_UNBN],PARS_IN[par_index[102]][TYPES_UNBN],PARS_IN[par_index[103]][TYPES_UNBN],\
                PARS_IN[par_index[104]][TYPES_UNBN],PARS_IN[par_index[105]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th7 = F4(TH7,PARS_IN[par_index[106]][TYPES_UNBN],PARS_IN[par_index[107]][TYPES_UNBN],PARS_IN[par_index[108]][TYPES_UNBN],\
                PARS_IN[par_index[109]][TYPES_UNBN],PARS_IN[par_index[110]][TYPES_UNBN],ZEROS_UNBN)
    
    f4_th8 = F4(TH8,PARS_IN[par_index[111]][TYPES_UNBN],PARS_IN[par_index[112]][TYPES_UNBN],PARS_IN[par_index[113]][TYPES_UNBN],\
                PARS_IN[par_index[114]][TYPES_UNBN],PARS_IN[par_index[115]][TYPES_UNBN],ZEROS_UNBN)
    
    
    EN_CRST_33_IN = f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    CRST_33_MOD_IN  = f4_th4
    
    
    f2 = F2(PARS_IN[par_index[116]][TYPES_UNBN], PARS_IN[par_index[117]][TYPES_UNBN], PARS_IN[par_index[118]][TYPES_UNBN], PARS_IN[par_index[119]][TYPES_UNBN], PARS_IN[par_index[120]][TYPES_UNBN],\
            PARS_IN[par_index[121]][TYPES_UNBN], PARS_IN[par_index[122]][TYPES_UNBN], PARS_IN[par_index[123]][TYPES_UNBN], PARS_IN[par_index[124]][TYPES_UNBN], ZEROS_UNBN)
    
    f4_th1 = F4(TH1,PARS_IN[par_index[125]][TYPES_UNBN],PARS_IN[par_index[126]][TYPES_UNBN],PARS_IN[par_index[127]][TYPES_UNBN],\
                PARS_IN[par_index[128]][TYPES_UNBN],PARS_IN[par_index[129]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th2 = F4(TH2,PARS_IN[par_index[130]][TYPES_UNBN],PARS_IN[par_index[131]][TYPES_UNBN],PARS_IN[par_index[132]][TYPES_UNBN],\
                PARS_IN[par_index[133]][TYPES_UNBN],PARS_IN[par_index[134]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th3 = F4(TH3,PARS_IN[par_index[135]][TYPES_UNBN],PARS_IN[par_index[136]][TYPES_UNBN],PARS_IN[par_index[137]][TYPES_UNBN],\
                PARS_IN[par_index[138]][TYPES_UNBN],PARS_IN[par_index[139]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th4 = F4(TH4_UNBN,PARS_IN[par_index[140]][TYPES_UNBN],PARS_IN[par_index[141]][TYPES_UNBN],PARS_IN[par_index[142]][TYPES_UNBN],\
                PARS_IN[par_index[143]][TYPES_UNBN],PARS_IN[par_index[144]][TYPES_UNBN],ZEROS_UNBN)
        
    f4_th7 = F4(TH7,PARS_IN[par_index[145]][TYPES_UNBN],PARS_IN[par_index[146]][TYPES_UNBN],PARS_IN[par_index[147]][TYPES_UNBN],\
                PARS_IN[par_index[148]][TYPES_UNBN],PARS_IN[par_index[149]][TYPES_UNBN],ZEROS_UNBN)
    
    f4_th8 = F4(TH8,PARS_IN[par_index[150]][TYPES_UNBN],PARS_IN[par_index[151]][TYPES_UNBN],PARS_IN[par_index[152]][TYPES_UNBN],\
                PARS_IN[par_index[153]][TYPES_UNBN],PARS_IN[par_index[154]][TYPES_UNBN],ZEROS_UNBN)
    
    EN_CRST_55_IN = f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    CRST_55_MOD_IN = f4_th4

    return
    
def compute_rew_factor(PARS) :
       
    #FENE
    EN_FENE_REW = -PARS[par_index[0]][TYPES_BN]/2.*torch.log( 1.-torch.square( FENE_R-PARS[par_index[1]][TYPES_BN] )/PARS[par_index[3]][TYPES_BN])
    
    
    #STACKING
    
    f1 = F1(STCK_R, PARS[par_index[44]][TYPES_BN], PARS[par_index[45]][TYPES_BN], PARS[par_index[46]][TYPES_BN], PARS[par_index[47]][TYPES_BN],PARS[par_index[48]][TYPES_BN],\
            PARS[par_index[49]][TYPES_BN], PARS[par_index[50]][TYPES_BN], PARS[par_index[51]][TYPES_BN], PARS[par_index[52]][TYPES_BN], PARS[par_index[53]][TYPES_BN], SHIFT_STCK[TYPES_BN], ZEROS_BN)
    
    #th4, th5, th6 Angular modulations
    
    
    f4_th4 = F4(TH4_BN,PARS[par_index[54]][TYPES_BN],PARS[par_index[55]][TYPES_BN],PARS[par_index[56]][TYPES_BN],\
                PARS[par_index[57]][TYPES_BN],PARS[par_index[58]][TYPES_BN],ZEROS_BN)
    
    f4_th5 = F4(TH5,PARS[par_index[59]][TYPES_BN],PARS[par_index[60]][TYPES_BN],PARS[par_index[61]][TYPES_BN],\
                PARS[par_index[62]][TYPES_BN],PARS[par_index[63]][TYPES_BN],ZEROS_BN)

    f4_th6 = F4(TH6,PARS[par_index[64]][TYPES_BN],PARS[par_index[65]][TYPES_BN],PARS[par_index[66]][TYPES_BN],\
                PARS[par_index[67]][TYPES_BN],PARS[par_index[68]][TYPES_BN],ZEROS_BN)
    
    #cosphi1, cosphi2 angular modulations

    STCK_MOD_REW = f1*f4_th4*f4_th5*f4_th6
     
    #HYDROGEN
    f4_th4 = F4(TH4_UNBN,PARS[par_index[29]][TYPES_UNBN],PARS[par_index[30]][TYPES_UNBN],PARS[par_index[31]][TYPES_UNBN],\
                PARS[par_index[32]][TYPES_UNBN],PARS[par_index[33]][TYPES_UNBN],ZEROS_UNBN)
        
    HYDR_MOD_REW = f4_th4

    #CROSS STACKING
    #33
    
    f4_th4 = F4(TH4_UNBN,PARS[par_index[101]][TYPES_UNBN],PARS[par_index[102]][TYPES_UNBN],PARS[par_index[103]][TYPES_UNBN],\
                PARS[par_index[104]][TYPES_UNBN],PARS[par_index[105]][TYPES_UNBN],ZEROS_UNBN)
        
    CRST_33_MOD_REW = f4_th4

        
    f4_th4 = F4(TH4_UNBN,PARS[par_index[140]][TYPES_UNBN],PARS[par_index[141]][TYPES_UNBN],PARS[par_index[142]][TYPES_UNBN],\
                PARS[par_index[143]][TYPES_UNBN],PARS[par_index[144]][TYPES_UNBN],ZEROS_UNBN)
        
    CRST_55_MOD_REW = f4_th4
        
    return EN_FENE_REW,  STCK_MOD_REW, HYDR_MOD_REW, CRST_33_MOD_REW, CRST_55_MOD_REW



def COST(PARS_OPTI) :
    
    global CURR_PARS


    #UPDATE PARAMETERS
    CURR_PARS = torch.where(UPDATE_MASK > 0, PARS_OPTI, CURR_PARS)

    #impose symmetries
    
    CURR_PARS[SYMM_LIST_SYMM] = CURR_PARS[SYMM_LIST]
    
    #impose continuity
    
    #delta
    CURR_PARS[3] = torch.square( CURR_PARS[2] )
    
    #f1
    EXP1 = torch.exp( -CURR_PARS[f1_A_ID]*(CURR_PARS[f1_RL_ID]-CURR_PARS[f1_R0_ID]) )
    EXP2 = torch.exp( -CURR_PARS[f1_A_ID]*(CURR_PARS[f1_RC_ID]-CURR_PARS[f1_R0_ID]) )
    EXP3 = torch.exp( -CURR_PARS[f1_A_ID]*(CURR_PARS[f1_RH_ID]-CURR_PARS[f1_R0_ID]) )

    CURR_PARS[f1_BL_ID] = torch.square( CURR_PARS[f1_A_ID] )*torch.square( EXP1*(1-EXP1) )/( torch.square(1-EXP1) - torch.square(1-EXP2) )
    CURR_PARS[f1_BH_ID] = torch.square( CURR_PARS[f1_A_ID] )*torch.square( EXP3*(1-EXP3) )/( torch.square(1-EXP3) - torch.square(1-EXP2) )
    
    CURR_PARS[f1_RCL_ID] = CURR_PARS[f1_RL_ID] - CURR_PARS[f1_A_ID]/CURR_PARS[f1_BL_ID]*( EXP1*(1-EXP1) )
    CURR_PARS[f1_RCH_ID] = CURR_PARS[f1_RH_ID] - CURR_PARS[f1_A_ID]/CURR_PARS[f1_BH_ID]*( EXP3*(1-EXP3) )
    
    #f4
    CURR_PARS[f4_TS_ID] = torch.sqrt(0.81225/CURR_PARS[f4_A_ID])
    CURR_PARS[f4_TC_ID] = 1./CURR_PARS[f4_A_ID]/CURR_PARS[f4_TS_ID]
    CURR_PARS[f4_B_ID] = CURR_PARS[f4_A_ID]*CURR_PARS[f4_TS_ID]/(CURR_PARS[f4_TC_ID]-CURR_PARS[f4_TS_ID])
    
    #REWEIGHTING
    #compute Boltzmann weights (reweighting)
    EN_FENE_REW,  STCK_MOD_REW, HYDR_MOD_REW, CRST_33_MOD_REW, CRST_55_MOD_REW = compute_rew_factor(CURR_PARS)

    WEIGHTS = torch.exp( 10*torch.sum( EN_FENE_REW + EN_STCK_IN*(1-STCK_MOD_REW/STCK_MOD_IN) + EN_HYDR_IN*(1-HYDR_MOD_REW/HYDR_MOD_IN) + EN_CRST_33_IN*(1-CRST_33_MOD_REW/CRST_33_MOD_IN) +\
                                     EN_CRST_55_IN*(1-CRST_55_MOD_REW/CRST_55_MOD_IN),dim=2 ))
    #normailise weights
    NORM = 1/WEIGHTS.sum(dim=1)
    WEIGHTS = NORM.unsqueeze(1)*WEIGHTS
    
    #compute reweighted average
    REW_MU = torch.sum(WEIGHTS.unsqueeze(2)*INTERNAL_COORDS,dim=1)  #sum because weights are normalised
    
    #compute reweighted covariance
    DIFF = INTERNAL_COORDS - REW_MU.unsqueeze(1)
    
    sx,sy,sz = INTERNAL_COORDS.size()
    DIFF = DIFF.reshape(sx*sy,sz)
    
    PROD = WEIGHTS.unsqueeze(2).unsqueeze(2)*torch.bmm(DIFF.unsqueeze(2), DIFF.unsqueeze(1)).reshape(sx,sy,sz,sz) #just magic
    REW_COV = PROD.sum(dim=1)
    
    #reduce reweighted covariances
    #XXX todo masks
    COV_RED_REW_COV = REW_COV[COV_RED_COV_MASK]
    MU_RED_REW_COV = REW_COV[MU_RED_COV_MASK]
    MU_RED_REW_MU = REW_MU[MU_RED_MU_MASK]
    
    DELTA = MU_RED_REW_MU - MU_RED_TARGET_MU
    #DELTA_T = torch.transpose(DELTA)
    
    #AVERAGE
    S = 0.5*torch.tensordot(DELTA, torch.tensordot(MU_RED_TARGET_COV_m1, DELTA, dims = ([1,1])), dims = ([1,1])).sum(dim=0)
    #inv_ex operations computes the inverse of a batch of square matrices. Does NOT syncronize with cpu
    S += 0.5*torch.tensordot(DELTA, torch.tensordot(torch.linalg.inv_ex(MU_RED_REW_COV), DELTA, dims = ([1,1])), dims = ([1,1])).sum(dim=0)
    
    #COVARIANCE
    
    #xxxtodo apply mask for diagonal covariance
    
    #diagonal(...) computes the trace of a batch of matrices
    S += 0.5*( torch.tensordot(COV_RED_TARGET_COV_m1,COV_RED_REW_COV, dims = ([1,2])).diagonal(offset=0, dim1=-1, dim2=-2).sum(-1) ).sum(dim=0)
    
    S += 0.5*( torch.tensordot(torch.linalg.inv_ex(COV_RED_REW_COV),COV_RED_TARGET_COV, dims = ([1,2])).diagonal(offset=0, dim1=-1, dim2=-2).sum(-1) ).sum(dim=0)
    
    return S

