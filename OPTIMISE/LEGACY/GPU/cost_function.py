#	#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:00:47 2024

@author: yqb22156
"""

import parameters_list as parl
import torch
import functions as fun
import config as cg
import numpy as np


#target_lb = 45.0 #target average long range bending persistence length (in nm).
#target_lt = 220.0 #target average long range torsional persistence length (in nm).

#target_C = 140 #target average long range C modulus
#target_Ar = 40 #target average long range Ar (A2) modulus

PARS_LIST = parl.PARS_LIST
par_index = parl.par_index

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

FENE_R = None
STCK_R = None
TH4_BN = None
TH5 = None
TH6 = None
COSPHI1 = None
COSPHI2 = None

TYPES_BN = None


HYDR_R = None
TH1 = None
TH2 = None
TH3 = None
TH4_UNBN = None
TH7 = None
TH8 = None

TYPES_UNBN_33 = None
TYPES_UNBN_55 = None

COS_THETAB = None
COS_OMEGAT = None


SHIFT_STCK = None
SHIFT_HYDR = None

PARS_IN = None
CURR_PARS = None
AVE_DELTA = None
ZEROS_BN = None
ONES_BN = None
ZEROS_UNBN = None


EN_FENE_IN = None
EN_STCK_IN = None
STCK_MOD_IN = None
STCK_MOD_FIX = None
EN_HYDR_IN = None
HYDR_MOD_IN = None
HYDR_MOD_FIX = None
EN_CRST_33_IN = None
CRST_33_MOD_IN = None
CRST_33_MOD_FIX = None
EN_CRST_55_IN = None
CRST_55_MOD_IN = None
CRST_55_MOD_FIX = None

UPDATE_MAP = None
SYMM_LIST = None #list with parameters to symmetrise - dim 1 is par index, dim 2 type
SYMM_LIST_SYMM = None #list of symmetric parameters - SYMM_LIST_SYMM[i][j] is the symmetric counterpart of SYMM_LIST[i][j]
CRST_TETRA_TYPE_33 = None
CRST_TETRA_TYPE_55 = None


internal_coords = []
INTERNAL_COORDS = None

save_mu = []
SAVE_MU = None

target_mu = []
target_cov = []
TARGET_MU = None
TARGET_COV = None

COV_RED_COV_MASK = None
MU_RED_COV_MASK = None
MU_RED_MU_MASK = None
COV_RED_COV_MASK_0 = None
MU_RED_COV_MASK_0 = None
MU_RED_MU_MASK_0 = None

COV_RED_TARGET_COV = None
COV_RED_TARGET_COV_m1 = None
MU_RED_TARGET_MU = None
MU_RED_TARGET_COV = None
MU_RED_TARGET_COV_m1 = None

COV_RED_REW_COV = None
COV_RED_REW_COV_m1 = None
AVE_COV_RED_TARGET_COV = None
MU_RED_REW_COV = None
MU_RED_REW_MU = None
MU_RED_REW_COV_m1 = None

COV_RED_IDS = None
MU_RED_IDS = None
PROP_IDS = None
LP_RED_IDS = None

LP_JK_IDS = None
LP_SUM_IDS_ROWS = None
LP_SUM_IDS_COLS = None

#initailise tensors with coordinates and
def init_tensors(dev, fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3,\
                  th4_unbn, th7, th8, types_unbn_33, types_unbn_55, costb, cosot, shifts, OXPS_zero) :

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
    global TH7
    global TH8
    global TYPES_UNBN_33
    global TYPES_UNBN_55

    global COS_THETAB
    global COS_OMEGAT

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

    #global UPDATE_MASK

    global internal_coords
    global INTERNAL_COORDS

    global save_mu
    global SAVE_MU

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

    TYPES_UNBN_33 = torch.tensor(types_unbn_33,device=device)
    TYPES_UNBN_55 = torch.tensor(types_unbn_55,device=device)

    COS_THETAB = torch.tensor(costb,device=device)
    COS_OMEGAT = torch.tensor(cosot,device=device)

    SHIFT_STCK = torch.tensor(shifts[1], device=device)
    SHIFT_HYDR = torch.tensor(shifts[0], device=device)

    PARS_IN = torch.tensor(OXPS_zero,device=device)
    CURR_PARS = torch.tensor(OXPS_zero,device=device)
    ZEROS_BN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    ONES_BN = torch.ones(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    ZEROS_UNBN = torch.zeros(len(types_unbn_33),len(types_unbn_33[0]),len(types_unbn_33[0][0]),device=device)

    EN_FENE_IN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    EN_STCK_IN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    STCK_MOD_IN = torch.zeros(len(types_bn),len(types_bn[0]),len(types_bn[0][0]),device=device)
    EN_HYDR_IN = torch.zeros(len(types_unbn_33),len(types_unbn_33[0]),len(types_unbn_33[0][0]),device=device)
    HYDR_MOD_IN = torch.zeros(len(types_unbn_33),len(types_unbn_33[0]),len(types_unbn_33[0][0]),device=device)
    EN_CRST_33_IN = torch.zeros(len(types_unbn_33),len(types_unbn_33[0]),len(types_unbn_33[0][0]),device=device)
    CRST_33_MOD_IN = torch.zeros(len(types_unbn_33),len(types_unbn_33[0]),len(types_unbn_33[0][0]),device=device)
    EN_CRST_55_IN = torch.zeros(len(types_unbn_33),len(types_unbn_33[0]),len(types_unbn_33[0][0]),device=device)
    CRST_55_MOD_IN = torch.zeros(len(types_unbn_33),len(types_unbn_33[0]),len(types_unbn_33[0][0]),device=device)

    #UPDATE_MASK = torch.zeros(len(PARS_LIST),256,device=device) # 1 if parameter is optimised, 0 otherwise

    INTERNAL_COORDS = torch.tensor(internal_coords,device=device)
    SAVE_MU = torch.tensor(save_mu,device=device)

    TARGET_MU = torch.tensor(target_mu,dtype=float,device=device)
    TARGET_COV = torch.tensor(target_cov,dtype=float,device=device)

def MORSE(R,EPS,R0,A) :
    return EPS*torch.square(1.-torch.exp(-(R-R0)*A))


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


f1_R0_ID = None
f1_A_ID = None
f1_RC_ID = None
f1_BL_ID = None
f1_BH_ID = None
f1_RL_ID = None
f1_RH_ID = None
f1_RCL_ID = None
f1_RCH_ID = None

OFFSET_f1_RC = None
OFFSET_f1_RL = None
OFFSET_f1_RH = None
OFFSET_f1_ID = None

f2_R0_ID = None
f2_A_ID = None
f2_RC_ID = None
f2_BL_ID = None
f2_BH_ID = None
f2_RL_ID = None
f2_RH_ID = None
f2_RCL_ID = None
f2_RCH_ID = None

OFFSET_f2_RC = None
OFFSET_f2_RL = None
OFFSET_f2_RH = None
OFFSET_f2_ID = None


f4_A_ID = None
f4_B_ID = None
f4_TS_ID = None
f4_TC_ID = None

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

    global OFFSET_f1_RC
    global OFFSET_f1_RL
    global OFFSET_f1_RH
    global OFFSET_f1_ID

    global f2_R0_ID
    global f2_RC_ID
    global f2_BL_ID
    global f2_BH_ID
    global f2_RL_ID
    global f2_RH_ID
    global f2_RCL_ID
    global f2_RCH_ID

    global OFFSET_f2_RC
    global OFFSET_f2_RL
    global OFFSET_f2_RH
    global OFFSET_f2_ID

    global f4_A_ID
    global f4_B_ID
    global f4_TS_ID
    global f4_TC_ID

    #f1
    f1_r0_id = [5,45]
    f1_a_id = [6,46]
    f1_rc_id = [7,47]
    f1_bl_id = [8,48]
    f1_bh_id = [9,49]
    f1_rl_id = [10,50]
    f1_rh_id = [11,51]
    f1_rcl_id = [12,52]
    f1_rch_id = [13,53]

    OFFSET_f1_RC = torch.tensor([[0.35]*256,[0.5]*256],device=device)
    OFFSET_f1_RL = torch.tensor([[-0.06]*256,[-0.08]*256],device=device)
    OFFSET_f1_RH = torch.tensor([[0.3]*256,[0.35]*256],device=device)
    OFFSET_f1_ID = torch.tensor([0,1],device=device)

    #f1
    f2_r0_id = [78,117]
    f2_rc_id = [79,118]
    f2_bl_id = [80,119]
    f2_bh_id = [81,120]
    f2_rl_id = [82,121]
    f2_rh_id = [83,122]
    f2_rcl_id = [84,123]
    f2_rch_id = [85,124]

    OFFSET_f2_RC = torch.tensor([[0.1]*256,[0.1]*256],device=device)
    OFFSET_f2_RL = torch.tensor([[-0.08]*256,[-0.08]*256],device=device)
    OFFSET_f2_RH = torch.tensor([[0.08]*256,[0.08]*256],device=device)
    OFFSET_f2_ID = torch.tensor([0,1],device=device)

    #f4
    f4_a_id = [15,20,25,30,35,40,55,60,65,87,92,97,102,107,112,126,131,136,141,146,151]
    f4_b_id = [16,21,26,31,36,41,56,61,66,88,93,98,103,108,113,127,132,137,142,147,152]
    f4_ts_id = [17,22,27,32,37,42,57,62,67,89,94,99,104,109,114,128,133,138,143,148,153]
    f4_tc_id = [18,23,28,33,38,43,58,63,68,90,95,100,105,110,115,129,134,139,144,149,154]

    f1_R0_ID = torch.tensor(f1_r0_id,device=device)
    f1_A_ID = torch.tensor(f1_a_id,device=device)
    f1_RC_ID = torch.tensor(f1_rc_id,device=device)
    f1_BL_ID = torch.tensor(f1_bl_id,device=device)
    f1_BH_ID = torch.tensor(f1_bh_id,device=device)
    f1_RL_ID = torch.tensor(f1_rl_id,device=device)
    f1_RH_ID = torch.tensor(f1_rh_id,device=device)
    f1_RCL_ID = torch.tensor(f1_rcl_id,device=device)
    f1_RCH_ID = torch.tensor(f1_rch_id,device=device)

    f2_R0_ID = torch.tensor(f2_r0_id,device=device)
    f2_RC_ID = torch.tensor(f2_rc_id,device=device)
    f2_BL_ID = torch.tensor(f2_bl_id,device=device)
    f2_BH_ID = torch.tensor(f2_bh_id,device=device)
    f2_RL_ID = torch.tensor(f2_rl_id,device=device)
    f2_RH_ID = torch.tensor(f2_rh_id,device=device)
    f2_RCL_ID = torch.tensor(f2_rcl_id,device=device)
    f2_RCH_ID = torch.tensor(f2_rch_id,device=device)

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
    elif vals_cn[0] == "HYDR" or (vals_cn[0] == "CRST" and vals_cn[1] == "K") :
        par_name = vals_cn[0]
        for i in range(1,len(vals_cn)-2):
            par_name += "_"+vals_cn[i]
        ty=fun.base_to_id(vals_cn[len(vals_cn)-2])*4+fun.base_to_id(vals_cn[len(vals_cn)-1])*4*4  #from base 4 to base 10
        OPT_PAR_LIST.append( [PARS_LIST.index(par_name),ty] )

    elif vals_cn[0] == "STCK" or vals_cn[0] == "FENE" or vals_cn[0] == "CRST":
        par_name = vals_cn[0]
        for i in range(1,len(vals_cn)-4):
            par_name += "_"+vals_cn[i]
        ty=fun.base_to_id(vals_cn[len(vals_cn)-4])+fun.base_to_id(vals_cn[len(vals_cn)-3])*4+fun.base_to_id(vals_cn[len(vals_cn)-2])*4*4+fun.base_to_id(vals_cn[len(vals_cn)-1])*4*4*4    #from base 4 to base 10
        OPT_PAR_LIST.append( [PARS_LIST.index(par_name),ty] )

    return


#UPDATE_MAP is the the tensor which tells which parameters are optimised
#SYMM_LIST and SYMM_LIST_SYMM are used to impose symmetries. i.e. SYMM_LIST_SYMM[i,j] = SYMM_LIST [i,j]
def build_masks_and_symm_tensors() :

    global OPT_PAR_LIST
    global UPDATE_MAP
    global device
    global SYMM_LIST
    global SYMM_LIST_SYMM
    global CRST_TETRA_TYPE_33
    global CRST_TETRA_TYPE_55

    global COV_RED_TARGET_COV
    global COV_RED_TARGET_COV_m1
    global AVE_COV_RED_TARGET_COV
    global MU_RED_TARGET_MU
    global MU_RED_TARGET_COV
    global MU_RED_TARGET_COV_m1

    global COV_RED_REW_COV
    global COV_RED_REW_COV_m1
    global MU_RED_REW_MU
    global MU_RED_REW_COV
    global MU_RED_REW_COV_m1

    #reduction masks (ids)
    global COV_RED_IDS
    global MU_RED_IDS
    global PROP_IDS

    #persistence length calculation masks
    global LP_RED_IDS
    global LP_JK_IDS
    global LP_SUM_IDS_ROWS
    global LP_SUM_IDS_COLS

    #reduction

    l = 0

    red_n = (cg.fin_j[l]-cg.in_j[l]+1)*(len(cg.ids)-len(cg.ids_cov))

    COV_RED_TARGET_COV = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)
    COV_RED_TARGET_COV_m1 = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)

    AVE_COV_RED_TARGET_COV = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)

    COV_RED_REW_COV = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)
    COV_RED_REW_COV_m1 = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)

    crids = []
    mrids = []
    lpids = []
    propids = []

    print(cg.ids)
    print(cg.ids_gs)
    print(cg.ids_cov)
    print(cg.ids_lp)


    for i in range(cg.dimension[0]):
        if cg.ids[i%len(cg.ids)] in cg.ids_cov: crids.append(i)
        if cg.ids[i%len(cg.ids)] in cg.ids_gs: mrids.append(i)
        if cg.ids[i%len(cg.ids)] in cg.ids_lp: lpids.append(i)
        if cg.ids[i%len(cg.ids)] == 1: propids.append(i)

    COV_RED_IDS = torch.tensor(crids,device=device)
    MU_RED_IDS = torch.tensor(mrids,device=device)
    LP_RED_IDS = torch.tensor(lpids,device=device)
    if 1 in cg.ids_gs: PROP_IDS = torch.tensor(propids,device=device)

    red_n = (cg.fin_j[l]-cg.in_j[l]+1)*(len(cg.ids)-len(cg.ids_gs))

    MU_RED_TARGET_MU = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,dtype=float,device=device)
    MU_RED_TARGET_COV = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)
    MU_RED_TARGET_COV_m1 = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)

    MU_RED_REW_MU = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,dtype=float,device=device)
    MU_RED_REW_COV = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)
    MU_RED_REW_COV_m1 = torch.zeros(cg.Nseq,cg.dimension[l]-red_n,cg.dimension[l]-red_n,dtype=float,device=device)


    #lp computation

    #N = (cg.fin_j[l]-cg.in_j[l]+1)*len(cg.ids_lp)

    #cg.lp_m = 3
    cg.lp_comp_range = 1- cg.in_j[l] + cg.fin_j[l] - cg.lp_m

    N = cg.lp_comp_range*len(cg.ids_lp)

    ids_jk = np.zeros((cg.lp_m,N),dtype=int)

    ids_rows = np.zeros((cg.Nseq,cg.lp_m,N,N),dtype=int)
    ids_cols = np.zeros((cg.Nseq,cg.lp_m,N,len(cg.ids_lp)),dtype=int)

    for m in range(cg.lp_m):
        for i in range(N) :
            ids_jk[m,i] = i+m*len(cg.ids_lp)

    for l in range(cg.Nseq) :
        for m in range(cg.lp_m):
            for i in range(N) :
                for j in range(len(cg.ids_lp)) :
                    ids_cols[l,m,i,j] = i%len(cg.ids_lp)
                for j in range(N) :
                    ids_rows[l,m,i,j] = j%len(cg.ids_lp)

    LP_JK_IDS = torch.tensor(ids_jk,device=device)
    LP_SUM_IDS_ROWS = torch.tensor(ids_rows,device=device)
    LP_SUM_IDS_COLS = torch.tensor(ids_cols,device=device)

    #parameters update

    sl = []
    sls = []
    umap = []


    for i in range(len(OPT_PAR_LIST)) :

        ID = OPT_PAR_LIST[i][0]
        TY = OPT_PAR_LIST[i][1]
        umap.append(ID*256+TY)

        if ID in parl.compl_symm_ids_2d :
            ty0 = TY%4
            ty1 = (TY//4)%4
            ty2 = (TY//4//4)%4
            ty3 = (TY//4//4//4)%4
            for l in range(4) :
                for m in range(4) :
                    TY_S = m+(3-ty2)*4+(3-ty1)*4*4+l*4*4*4
                    if TY != TY_S:
                        #sl.append(ID*256+TY)
                        sl.append(i)
                        sls.append(ID*256+TY_S)
        if ID in parl.compl_symm_ids_2d_no_TT_AA :
            ty0 = TY%4
            ty1 = (TY//4)%4
            ty2 = (TY//4//4)%4
            ty3 = (TY//4//4//4)%4
            if (ty1 == 0 and ty3 == 0) or (ty1 == 3 and ty3 == 3):
                continue
            for l in range(4) :
                for m in range(4) :
                    TY_S = m+(3-ty2)*4+(3-ty1)*4*4+l*4*4*4
                    if TY != TY_S:
                        #sl.append(ID*256+TY)
                        sl.append(i)
                        sls.append(ID*256+TY_S)

        if ID in parl.compl_symm_ids :
            ty0 = TY%4
            ty1 = (TY//4)%4
            ty2 = (TY//4//4)%4
            ty3 = (TY//4//4//4)%4
            TY_S = (3-ty3)+(3-ty2)*4+(3-ty1)*4*4+(3-ty0)*4*4*4
            if TY != TY_S:
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(ID*256+TY_S)
            if ID in parl.is_th2 :
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY)
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY_S)
            if ID in parl.is_th5 :
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY)
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY_S)
            if ID in parl.is_th7 :
                #sl.append(256*ID+TY)
                sl.append(i)
                sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY)
                #sl.append(256*ID+TY)
                sl.append(i)
                sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY_S)

        if ID in parl.perm_symm_ids_2d :
            ty0 = TY%4
            ty1 = (TY//4)%4
            ty2 = (TY//4//4)%4
            ty3 = (TY//4//4//4)%4
            for l in range(4) :
                for m in range(4) :
                    TY_S = m+(ty2)*4+(ty1)*4*4+l*4*4*4
                    if TY != TY_S:
                        #sl.append(ID*256+TY)
                        sl.append(i)
                        sls.append(ID*256+TY_S)
                    if ID in parl.is_th2 :
                        #sl.append(ID*256+TY)
                        sl.append(i)
                        sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY)
                        #sl.append(ID*256+TY)
                        sl.append(i)
                        sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY_S)
                    if ID in parl.is_th5 :
                        #sl.append(ID*256+TY)
                        sl.append(i)
                        sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY)
                        #sl.append(ID*256+TY)
                        sl.append(i)
                        sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY_S)
                    if ID in parl.is_th7 :
                        #sl.append(256*ID+TY)
                        sl.append(i)
                        sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY)
                        #sl.append(256*ID+TY)
                        sl.append(i)
                        sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY_S)

        if ID in parl.perm_symm_ids_4d :
            ty0 = TY%4
            ty1 = (TY//4)%4
            ty2 = (TY//4//4)%4
            ty3 = (TY//4//4//4)%4
            TY_S = ty3+ty2*4+ty1*4*4+ty0*4*4*4
            if TY != TY_S:
                sl.append(i)
                sls.append(ID*256+TY_S)
            if ID in parl.is_th2 :
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY)
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY_S)
            if ID in parl.is_th5 :
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY)
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY_S)
            if ID in parl.is_th7 :
                #sl.append(256*ID+TY)
                sl.append(i)
                sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY)
                #sl.append(256*ID+TY)
                sl.append(i)
                sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY_S)

        if ID in parl.perm_symm_ids :
            ty0 = TY%4
            ty1 = (TY//4)%4
            ty2 = (TY//4//4)%4
            ty3 = (TY//4//4//4)%4
            TY_S = ty3+ty2*4+ty1*4*4+ty0*4*4*4
            if TY != TY_S:
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(ID*256+TY_S)
            if ID in parl.is_th2 :
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY)
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th3[parl.is_th2.index(ID)]+TY_S)
            if ID in parl.is_th5 :
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY)
                #sl.append(ID*256+TY)
                sl.append(i)
                sls.append(256*parl.is_th6[parl.is_th5.index(ID)]+TY_S)
            if ID in parl.is_th7 :
                #sl.append(256*ID+TY)
                sl.append(i)
                sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY)
                #sl.append(256*ID+TY)
                sl.append(i)
                sls.append(256*parl.is_th8[parl.is_th7.index(ID)]+TY_S)

    UPDATE_MAP = torch.tensor(umap,device=device)
    SYMM_LIST = torch.tensor(sl,device=device)
    SYMM_LIST_SYMM = torch.tensor(sls,device=device)


    # Nota:
    # stacking - qn3, q, p, pn5
    # cross_33 - qn5, q, p, pn5
    # corss_55 - qn3, q, p, pn3
    # p = chosen particle, q = interacting particle

    sl_33 = []
    sl_55 = []
    for m in range(4):
        for n in range(4):
             for p in range(4):
                 for q in range(4):
                     sl_33.append((3-m)+(3-n)*4+p*16+q*64)
                     sl_55.append(q+p*4+(3-n)*16+(3-m)*64)

    CRST_TETRA_TYPE_33 = torch.tensor(sl_33,device=device)
    CRST_TETRA_TYPE_55 = torch.tensor(sl_55,device=device)

    return



#PREPARE REDUCED TARGETS AND COV^{-1}

IDS_AVE_COV_SUM = None
IDS_AVE_COV_EXPAND = None

def reduce_targets() :

    global COV_RED_TARGET_COV
    global COV_RED_TARGET_COV_m1
    global IDS_AVE_COV_SUM
    global IDS_AVE_COV_EXPAND
    global AVE_COV_RED_TARGET_COV
    global MU_RED_TARGET_MU
    global MU_RED_TARGET_COV
    global MU_RED_TARGET_COV_m1

    TMP = TARGET_COV[:,:,COV_RED_IDS]
    COV_RED_TARGET_COV = TMP[:,COV_RED_IDS]

    #TMP = torch.zeros(cg.Nseq,cg.ids_cov,device=device)

    nentries = int(len(COV_RED_TARGET_COV[0])/len(cg.ids_cov))

    ids_s = np.zeros((len(cg.ids_cov), nentries), dtype=int)

    for l in range(len(cg.ids_cov)):
       for i in range(nentries) :
            ids_s[l,i] = l+i*len(cg.ids_cov)

    IDS_AVE_COV_SUM = torch.tensor(ids_s,device=device)

    TMP = torch.diagonal(COV_RED_TARGET_COV, 0, 1, 2)[:,IDS_AVE_COV_SUM].mean(dim=2)

    ids_e = []

    for l in range(len(COV_RED_TARGET_COV[0])) :
        ids_e.append(l%len(cg.ids_cov))

    IDS_AVE_COV_EXPAND = torch.tensor(ids_e,device=device)

    AVE_COV_RED_TARGET_COV = TMP[:,IDS_AVE_COV_EXPAND]

    TMP = TARGET_COV[:,:,MU_RED_IDS]
    MU_RED_TARGET_COV = TMP[:,MU_RED_IDS]
    MU_RED_TARGET_MU = TARGET_MU[:,MU_RED_IDS]

    #print(len(TARGET_COV), len(TARGET_COV[0]), len(TARGET_COV[0][0]))
    #print(len(COV_RED_TARGET_COV), len(COV_RED_TARGET_COV[0]), len(COV_RED_TARGET_COV[0][0]))
    #print(len(TARGET_COV), len(TARGET_COV[0]), len(TARGET_COV[0][0]))
    #print(len(MU_RED_TARGET_COV), len(MU_RED_TARGET_COV[0]), len(MU_RED_TARGET_COV[0][0]))

    #print(COV_RED_TARGET_COV[0])
    #print(MU_RED_TARGET_COV[0])

    sx,sy,sz = COV_RED_TARGET_COV.size()
    COV_RED_TARGET_COV_m1 = torch.zeros(sx,sy,sz,dtype=float,device=device)

    sx,sy,sz = MU_RED_TARGET_COV.size()
    MU_RED_TARGET_COV_m1 = torch.zeros(sx,sy,sz,dtype=float,device=device)

    for l in range(cg.Nseq) :
        for i in range(len(COV_RED_TARGET_COV[l])):
            COV_RED_TARGET_COV_m1[l][i][i] = 1/COV_RED_TARGET_COV[l][i][i]
        for i in range(len(MU_RED_TARGET_COV[l])):
            MU_RED_TARGET_COV_m1[l][i][i] = 1/MU_RED_TARGET_COV[l][i][i]

        #COV_RED_TARGET_COV_m1[l] = torch.linalg.inv( COV_RED_TARGET_COV[l] )
        #MU_RED_TARGET_COV_m1[l] = torch.linalg.inv( MU_RED_TARGET_COV[l] )

    return

#COMPUTE INITIAL ENERGY AND MODULATIONS

def compute_initial_energy() :

    global EN_FENE_IN
    global EN_STCK_IN
    global STCK_MOD_IN
    global STCK_MOD_FIX
    global EN_HYDR_IN
    global HYDR_MOD_IN
    global HYDR_MOD_FIX
    global EN_CRST_33_IN
    global CRST_33_MOD_IN
    global CRST_33_MOD_FIX
    global EN_CRST_55_IN
    global CRST_55_MOD_IN
    global CRST_55_MOD_FIX

    #FENE
    EN_FENE_IN = -PARS_IN[par_index[0]][TYPES_BN]/2.*torch.log( 1.-torch.square( FENE_R-PARS_IN[par_index[1]][TYPES_BN] )/PARS_IN[par_index[3]][TYPES_BN])

    #STACKING
    #radial part
    #piecewise energy: compute all possibilities, and filter with where

    f1 = F1(STCK_R, PARS_IN[par_index[44]][TYPES_BN], PARS_IN[par_index[45]][TYPES_BN], PARS_IN[par_index[46]][TYPES_BN], PARS_IN[par_index[48]][TYPES_BN],\
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
    STCK_MOD_FIX = f5_phi1*f5_phi2

    #HYDROGEN
    #TYPES_UNBN_33 and TYPES_UNBN_55 are the same here: hydr is a 2d interaction, flanking types are irrelevant
    f1 = F1(HYDR_R, PARS_IN[par_index[4]][TYPES_UNBN_33], PARS_IN[par_index[5]][TYPES_UNBN_33], PARS_IN[par_index[6]][TYPES_UNBN_33],PARS_IN[par_index[8]][TYPES_UNBN_33],\
            PARS_IN[par_index[9]][TYPES_UNBN_33], PARS_IN[par_index[10]][TYPES_UNBN_33], PARS_IN[par_index[11]][TYPES_UNBN_33], PARS_IN[par_index[12]][TYPES_UNBN_33], PARS_IN[par_index[13]][TYPES_UNBN_33], SHIFT_HYDR[TYPES_UNBN_33], ZEROS_UNBN)

    f4_th1 = F4(TH1,PARS_IN[par_index[14]][TYPES_UNBN_33],PARS_IN[par_index[15]][TYPES_UNBN_33],PARS_IN[par_index[16]][TYPES_UNBN_33],\
                PARS_IN[par_index[17]][TYPES_UNBN_33],PARS_IN[par_index[18]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th2 = F4(TH2,PARS_IN[par_index[19]][TYPES_UNBN_33],PARS_IN[par_index[20]][TYPES_UNBN_33],PARS_IN[par_index[21]][TYPES_UNBN_33],\
                PARS_IN[par_index[22]][TYPES_UNBN_33],PARS_IN[par_index[23]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th3 = F4(TH3,PARS_IN[par_index[24]][TYPES_UNBN_33],PARS_IN[par_index[25]][TYPES_UNBN_33],PARS_IN[par_index[26]][TYPES_UNBN_33],\
                PARS_IN[par_index[27]][TYPES_UNBN_33],PARS_IN[par_index[28]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th4 = F4(TH4_UNBN,PARS_IN[par_index[29]][TYPES_UNBN_33],PARS_IN[par_index[30]][TYPES_UNBN_33],PARS_IN[par_index[31]][TYPES_UNBN_33],\
                PARS_IN[par_index[32]][TYPES_UNBN_33],PARS_IN[par_index[33]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th7 = F4(TH7,PARS_IN[par_index[34]][TYPES_UNBN_33],PARS_IN[par_index[35]][TYPES_UNBN_33],PARS_IN[par_index[36]][TYPES_UNBN_33],\
                PARS_IN[par_index[37]][TYPES_UNBN_33],PARS_IN[par_index[38]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th8 = F4(TH8,PARS_IN[par_index[39]][TYPES_UNBN_33],PARS_IN[par_index[40]][TYPES_UNBN_33],PARS_IN[par_index[41]][TYPES_UNBN_33],\
                PARS_IN[par_index[42]][TYPES_UNBN_33],PARS_IN[par_index[43]][TYPES_UNBN_33],ZEROS_UNBN)

    EN_HYDR_IN = f1*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    HYDR_MOD_IN  = f4_th4
    HYDR_MOD_FIX = f1*f4_th1*f4_th2*f4_th3*f4_th7*f4_th8

    #CROSS STACKING
    #33
    f2 = F2(HYDR_R, PARS_IN[par_index[77]][TYPES_UNBN_33], PARS_IN[par_index[78]][TYPES_UNBN_33], PARS_IN[par_index[79]][TYPES_UNBN_33], PARS_IN[par_index[80]][TYPES_UNBN_33], PARS_IN[par_index[81]][TYPES_UNBN_33],\
            PARS_IN[par_index[82]][TYPES_UNBN_33], PARS_IN[par_index[83]][TYPES_UNBN_33], PARS_IN[par_index[84]][TYPES_UNBN_33], PARS_IN[par_index[85]][TYPES_UNBN_33], ZEROS_UNBN)

    f4_th1 = F4(TH1,PARS_IN[par_index[86]][TYPES_UNBN_33],PARS_IN[par_index[87]][TYPES_UNBN_33],PARS_IN[par_index[88]][TYPES_UNBN_33],\
                PARS_IN[par_index[89]][TYPES_UNBN_33],PARS_IN[par_index[90]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th2 = F4(TH2,PARS_IN[par_index[91]][TYPES_UNBN_33],PARS_IN[par_index[92]][TYPES_UNBN_33],PARS_IN[par_index[93]][TYPES_UNBN_33],\
                PARS_IN[par_index[94]][TYPES_UNBN_33],PARS_IN[par_index[95]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th3 = F4(TH3,PARS_IN[par_index[96]][TYPES_UNBN_33],PARS_IN[par_index[97]][TYPES_UNBN_33],PARS_IN[par_index[98]][TYPES_UNBN_33],\
                PARS_IN[par_index[99]][TYPES_UNBN_33],PARS_IN[par_index[100]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th4 = F4(TH4_UNBN,PARS_IN[par_index[101]][TYPES_UNBN_33],PARS_IN[par_index[102]][TYPES_UNBN_33],PARS_IN[par_index[103]][TYPES_UNBN_33],\
                PARS_IN[par_index[104]][TYPES_UNBN_33],PARS_IN[par_index[105]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th7 = F4(TH7,PARS_IN[par_index[106]][TYPES_UNBN_33],PARS_IN[par_index[107]][TYPES_UNBN_33],PARS_IN[par_index[108]][TYPES_UNBN_33],\
                PARS_IN[par_index[109]][TYPES_UNBN_33],PARS_IN[par_index[110]][TYPES_UNBN_33],ZEROS_UNBN)

    f4_th8 = F4(TH8,PARS_IN[par_index[111]][TYPES_UNBN_33],PARS_IN[par_index[112]][TYPES_UNBN_33],PARS_IN[par_index[113]][TYPES_UNBN_33],\
                PARS_IN[par_index[114]][TYPES_UNBN_33],PARS_IN[par_index[115]][TYPES_UNBN_33],ZEROS_UNBN)

    EN_CRST_33_IN = f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    CRST_33_MOD_IN  = f2*f4_th4
    CRST_33_MOD_FIX  = f4_th1*f4_th2*f4_th3*f4_th7*f4_th8


    #55
    f2 = F2(HYDR_R, PARS_IN[par_index[116]][TYPES_UNBN_55], PARS_IN[par_index[117]][TYPES_UNBN_55], PARS_IN[par_index[118]][TYPES_UNBN_55], PARS_IN[par_index[119]][TYPES_UNBN_55], PARS_IN[par_index[120]][TYPES_UNBN_55],\
            PARS_IN[par_index[121]][TYPES_UNBN_55], PARS_IN[par_index[122]][TYPES_UNBN_55], PARS_IN[par_index[123]][TYPES_UNBN_55], PARS_IN[par_index[124]][TYPES_UNBN_55], ZEROS_UNBN)

    f4_th1 = F4(TH1,PARS_IN[par_index[125]][TYPES_UNBN_55],PARS_IN[par_index[126]][TYPES_UNBN_55],PARS_IN[par_index[127]][TYPES_UNBN_55],\
                PARS_IN[par_index[128]][TYPES_UNBN_55],PARS_IN[par_index[129]][TYPES_UNBN_55],ZEROS_UNBN)

    f4_th2 = F4(TH2,PARS_IN[par_index[130]][TYPES_UNBN_55],PARS_IN[par_index[131]][TYPES_UNBN_55],PARS_IN[par_index[132]][TYPES_UNBN_55],\
                PARS_IN[par_index[133]][TYPES_UNBN_55],PARS_IN[par_index[134]][TYPES_UNBN_55],ZEROS_UNBN)

    f4_th3 = F4(TH3,PARS_IN[par_index[135]][TYPES_UNBN_55],PARS_IN[par_index[136]][TYPES_UNBN_55],PARS_IN[par_index[137]][TYPES_UNBN_55],\
                PARS_IN[par_index[138]][TYPES_UNBN_55],PARS_IN[par_index[139]][TYPES_UNBN_55],ZEROS_UNBN)

    f4_th4 = F4(TH4_UNBN,PARS_IN[par_index[140]][TYPES_UNBN_55],PARS_IN[par_index[141]][TYPES_UNBN_55],PARS_IN[par_index[142]][TYPES_UNBN_55],\
                PARS_IN[par_index[143]][TYPES_UNBN_55],PARS_IN[par_index[144]][TYPES_UNBN_55],ZEROS_UNBN)

    f4_th7 = F4(TH7,PARS_IN[par_index[145]][TYPES_UNBN_55],PARS_IN[par_index[146]][TYPES_UNBN_55],PARS_IN[par_index[147]][TYPES_UNBN_55],\
                PARS_IN[par_index[148]][TYPES_UNBN_55],PARS_IN[par_index[149]][TYPES_UNBN_55],ZEROS_UNBN)

    f4_th8 = F4(TH8,PARS_IN[par_index[150]][TYPES_UNBN_55],PARS_IN[par_index[151]][TYPES_UNBN_55],PARS_IN[par_index[152]][TYPES_UNBN_55],\
                PARS_IN[par_index[153]][TYPES_UNBN_55],PARS_IN[par_index[154]][TYPES_UNBN_55],ZEROS_UNBN)

    EN_CRST_55_IN = f2*f4_th1*f4_th2*f4_th3*f4_th4*f4_th7*f4_th8
    CRST_55_MOD_IN = f2*f4_th4
    CRST_55_MOD_FIX  = f4_th1*f4_th2*f4_th3*f4_th7*f4_th8

    return



def compute_rew_factor(PARS,SH_ST) :

    #FENE
    EN_FENE_REW = -PARS[par_index[0]][TYPES_BN]/2.*torch.log( 1.-torch.square( FENE_R-PARS[par_index[1]][TYPES_BN] )/PARS[par_index[2]][TYPES_BN]/PARS[par_index[2]][TYPES_BN])

    #FENE_ENERGY+CONTINUITY-TO AVOID DIVERGENT ENERGY

    h = 0.01
    LOWA = -PARS[par_index[0]][TYPES_BN]*2*(PARS[par_index[2]][TYPES_BN]-h)/h*(2*PARS[par_index[2]][TYPES_BN]-h)
    LOWB = -PARS[par_index[0]][TYPES_BN]*torch.log( 1.-torch.square( PARS[par_index[2]][TYPES_BN]-h )/PARS[par_index[2]][TYPES_BN]/PARS[par_index[2]][TYPES_BN])\
          -LOWA*(PARS[par_index[1]][TYPES_BN]-PARS[par_index[2]][TYPES_BN]+h)

    UPA = -LOWA
    UPB = -PARS[par_index[0]][TYPES_BN]*torch.log( 1.-torch.square( PARS[par_index[2]][TYPES_BN]-h )/PARS[par_index[2]][TYPES_BN]/PARS[par_index[2]][TYPES_BN])\
          -UPA*(PARS[par_index[1]][TYPES_BN]+PARS[par_index[2]][TYPES_BN]-h)


    EN_FENE_TMP = -PARS[par_index[0]][TYPES_BN]/2.*torch.log( 1.-torch.square( FENE_R-PARS[par_index[1]][TYPES_BN] )/PARS[par_index[2]][TYPES_BN]/PARS[par_index[2]][TYPES_BN])
    EN_FENE_REW = torch.where(FENE_R > PARS[par_index[1]][TYPES_BN]-PARS[par_index[2]][TYPES_BN]+h, EN_FENE_TMP, LOWA*FENE_R+LOWB)
    EN_FENE_REW = torch.where(FENE_R > PARS[par_index[1]][TYPES_BN]+PARS[par_index[2]][TYPES_BN]-h,UPA*FENE_R+UPB, EN_FENE_REW)

    #STACKING

    f1 = F1(STCK_R, PARS[par_index[44]][TYPES_BN], PARS[par_index[45]][TYPES_BN], PARS[par_index[46]][TYPES_BN], PARS[par_index[48]][TYPES_BN],\
            PARS[par_index[49]][TYPES_BN], PARS[par_index[50]][TYPES_BN], PARS[par_index[51]][TYPES_BN], PARS[par_index[52]][TYPES_BN], PARS[par_index[53]][TYPES_BN], SH_ST[TYPES_BN], ZEROS_BN)

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
    f4_th4 = F4(TH4_UNBN,PARS[par_index[29]][TYPES_UNBN_33],PARS[par_index[30]][TYPES_UNBN_33],PARS[par_index[31]][TYPES_UNBN_33],\
                PARS[par_index[32]][TYPES_UNBN_33],PARS[par_index[33]][TYPES_UNBN_33],ZEROS_UNBN)

    HYDR_MOD_REW = f4_th4

    #CROSS STACKING
    #33

    f2 = F2(HYDR_R, PARS_IN[par_index[77]][TYPES_UNBN_33], PARS_IN[par_index[78]][TYPES_UNBN_33], PARS_IN[par_index[79]][TYPES_UNBN_33], PARS_IN[par_index[80]][TYPES_UNBN_33], PARS_IN[par_index[81]][TYPES_UNBN_33],\
            PARS_IN[par_index[82]][TYPES_UNBN_33], PARS_IN[par_index[83]][TYPES_UNBN_33], PARS_IN[par_index[84]][TYPES_UNBN_33], PARS_IN[par_index[85]][TYPES_UNBN_33], ZEROS_UNBN)

    f4_th4 = F4(TH4_UNBN,PARS[par_index[101]][TYPES_UNBN_33],PARS[par_index[102]][TYPES_UNBN_33],PARS[par_index[103]][TYPES_UNBN_33],\
                PARS[par_index[104]][TYPES_UNBN_33],PARS[par_index[105]][TYPES_UNBN_33],ZEROS_UNBN)

    CRST_33_MOD_REW = f2*f4_th4


    f2 = F2(HYDR_R, PARS_IN[par_index[116]][TYPES_UNBN_55], PARS_IN[par_index[117]][TYPES_UNBN_55], PARS_IN[par_index[118]][TYPES_UNBN_55], PARS_IN[par_index[119]][TYPES_UNBN_55], PARS_IN[par_index[120]][TYPES_UNBN_55],\
            PARS_IN[par_index[121]][TYPES_UNBN_55], PARS_IN[par_index[122]][TYPES_UNBN_55], PARS_IN[par_index[123]][TYPES_UNBN_55], PARS_IN[par_index[124]][TYPES_UNBN_55], ZEROS_UNBN)

    f4_th4 = F4(TH4_UNBN,PARS[par_index[140]][TYPES_UNBN_55],PARS[par_index[141]][TYPES_UNBN_55],PARS[par_index[142]][TYPES_UNBN_55],\
                PARS[par_index[143]][TYPES_UNBN_55],PARS[par_index[144]][TYPES_UNBN_55],ZEROS_UNBN)

    CRST_55_MOD_REW = f2*f4_th4

    return EN_FENE_REW,  STCK_MOD_REW, HYDR_MOD_REW, CRST_33_MOD_REW, CRST_55_MOD_REW



def compute_stiffness(cov) :

    #print("Evaluating stiffness")
    #torch.set_printoptions(profile="full")

    #print(cov[0])

    M = torch.zeros(3,3,device=device)

    M_i = torch.zeros((cg.Nseq,cg.lp_m,len(cg.ids_lp),len(cg.ids_lp)),device=device)

    n = cg.lp_comp_range*len(cg.ids_lp)
    M_jk = torch.zeros(cg.Nseq,cg.lp_m,n,n,device=device)
    #TMP = cov[LP_JK_IDS]
    TMP = torch.zeros((cg.Nseq,cg.lp_m,n,len(cov[0,0])),device=device)
    for l in range(cg.lp_m) : #can we avoid a for loop here? It's not too bad, though: lp_m = 3
        TMP[:,l,:,:] = cov[:,LP_JK_IDS[l],:]
        M_jk[:,l,:,:] = TMP[:,l,:,LP_JK_IDS[l]]

    #print(M_jk[0])

    TMP1 = torch.zeros((cg.Nseq,cg.lp_m,n,len(cg.ids_lp)),device=device)

    TMP1.scatter_add_(3,LP_SUM_IDS_ROWS,M_jk)
    M_i.scatter_add_(2,LP_SUM_IDS_COLS,TMP1)

    #print(M_i[0])
    M_i /= (25.*0.34*0.34)

    M = (torch.linalg.inv(M_i)).sum(dim=1)/0.34/cg.lp_m*cg.lp_comp_range

    return M

PAR0 = []
LB = None
LT = None
CURR_AVE_DELTA = 0.2
ave_prop_flag = True

def COST(PARS) :


    PARS_OPTI = torch.tensor(PARS,device=device)
    #PARS_OPTI = torch.clone(PARS)


    CURR_PARS = torch.clone(PAR0)
    SH_ST = torch.clone(SHIFT_STCK)
    #SH_HY = torch.clone(SHIFT_HYDR)


    #UPDATE PARAMETERS
    CURR_PARS.put_(UPDATE_MAP, PARS_OPTI)

    #impose symmetries
    VALS = torch.gather( torch.reshape(PARS_OPTI,(-1,)),0,SYMM_LIST )
    CURR_PARS.put_(SYMM_LIST_SYMM,VALS)

    #crst r0
    #crst_r0 = sqrt( stck_r0^2+hydr_r0^2/2*(1+cos(2asin(sqrt(fene_ro^2-stck_r0^2)))) )

    #Constraints - no continuity
    """
    fene_r02_crst = torch.square(CURR_PARS[1][CRST_TETRA_TYPE_33])
    stck_r02_crst = torch.square(CURR_PARS[45][CRST_TETRA_TYPE_33])

    #0.08 <- if hydr_r0 = 0.4. otherwise this has to be changed
    CURR_PARS[78] = torch.sqrt( stck_r02_crst+0.08*(1+torch.cos(2*torch.arcsin(0.5*torch.sqrt(fene_r02_crst-stck_r02_crst)))) )

    fene_r02_crst = torch.square(CURR_PARS[1][CRST_TETRA_TYPE_55])
    stck_r02_crst = torch.square(CURR_PARS[45][CRST_TETRA_TYPE_55])

    CURR_PARS[117] = torch.sqrt( stck_r02_crst+0.08*(1+torch.cos(2*torch.arcsin(0.5*torch.sqrt(fene_r02_crst-stck_r02_crst)))) )
    """

    #Fix delta average

    global CURR_AVE_DELTA
    aveD = torch.mean(CURR_PARS[2])
    CURR_PARS[2] = CURR_PARS[2] - aveD + AVE_DELTA
    CURR_AVE_DELTA = torch.mean(CURR_PARS[2])

    #Make STCK_THETA5A average

    aveTH5A = torch.mean(CURR_PARS[60])
    CURR_PARS[60] = aveTH5A

    #delta2
    #delta2 = delta^2
    CURR_PARS[3] = torch.square( CURR_PARS[2] )

    #Constraints - continuity
    #f1

    CURR_PARS[f1_RL_ID] = OFFSET_f1_RL[OFFSET_f1_ID] + CURR_PARS[f1_R0_ID]
    CURR_PARS[f1_RH_ID] = OFFSET_f1_RH[OFFSET_f1_ID] + CURR_PARS[f1_R0_ID]
    CURR_PARS[f1_RC_ID] = OFFSET_f1_RC[OFFSET_f1_ID] + CURR_PARS[f1_R0_ID]

    EXP1 = torch.exp( -CURR_PARS[f1_A_ID]*(CURR_PARS[f1_RL_ID]-CURR_PARS[f1_R0_ID]) )
    EXP2 = torch.exp( -CURR_PARS[f1_A_ID]*(CURR_PARS[f1_RC_ID]-CURR_PARS[f1_R0_ID]) )
    EXP3 = torch.exp( -CURR_PARS[f1_A_ID]*(CURR_PARS[f1_RH_ID]-CURR_PARS[f1_R0_ID]) )

    #print("EXP ok")

    CURR_PARS[f1_BL_ID] = torch.square( CURR_PARS[f1_A_ID] )*torch.square( EXP1*(1-EXP1) )/( torch.square(1-EXP1) - torch.square(1-EXP2) )
    CURR_PARS[f1_BH_ID] = torch.square( CURR_PARS[f1_A_ID] )*torch.square( EXP3*(1-EXP3) )/( torch.square(1-EXP3) - torch.square(1-EXP2) )

    CURR_PARS[f1_RCL_ID] = CURR_PARS[f1_RL_ID] - CURR_PARS[f1_A_ID]/CURR_PARS[f1_BL_ID]*( EXP1*(1-EXP1) )
    CURR_PARS[f1_RCH_ID] = CURR_PARS[f1_RH_ID] - CURR_PARS[f1_A_ID]/CURR_PARS[f1_BH_ID]*( EXP3*(1-EXP3) )

    #update shift

    #SH_HY = MORSE(CURR_PARS[7],CURR_PARS[4],CURR_PARS[5],CURR_PARS[6])
    SH_ST = CURR_PARS[44]*torch.square(1.-torch.exp(-(CURR_PARS[47]-CURR_PARS[45])*CURR_PARS[46]))

    #f2
    CURR_PARS[f2_RL_ID] = OFFSET_f2_RL[OFFSET_f2_ID] + CURR_PARS[f2_R0_ID]
    CURR_PARS[f2_RH_ID] = OFFSET_f2_RH[OFFSET_f2_ID] + CURR_PARS[f2_R0_ID]
    CURR_PARS[f2_RC_ID] = OFFSET_f2_RC[OFFSET_f2_ID] + CURR_PARS[f2_R0_ID]

    TERM1 = CURR_PARS[f2_RL_ID]-CURR_PARS[f2_R0_ID]
    TERM2 = CURR_PARS[f2_RH_ID]-CURR_PARS[f2_R0_ID]
    TERM3 = CURR_PARS[f2_RC_ID]-CURR_PARS[f2_R0_ID]

    CURR_PARS[f2_RCL_ID] = CURR_PARS[f2_RL_ID] - (torch.square(TERM1)-torch.square(TERM3))/TERM1
    CURR_PARS[f2_RCH_ID] = CURR_PARS[f2_RH_ID] - (torch.square(TERM2)-torch.square(TERM3))/TERM2

    CURR_PARS[f2_BL_ID] = -0.5*TERM1/( CURR_PARS[f2_RCL_ID]-CURR_PARS[f2_RL_ID] )
    CURR_PARS[f2_BH_ID] = -0.5*TERM2/( CURR_PARS[f2_RCH_ID]-CURR_PARS[f2_RH_ID] )

    #f4
    CURR_PARS[f4_TS_ID] = torch.sqrt(0.81225/CURR_PARS[f4_A_ID])
    CURR_PARS[f4_TC_ID] = 1./CURR_PARS[f4_A_ID]/CURR_PARS[f4_TS_ID]
    CURR_PARS[f4_B_ID] = CURR_PARS[f4_A_ID]*CURR_PARS[f4_TS_ID]/(CURR_PARS[f4_TC_ID]-CURR_PARS[f4_TS_ID])

    #REWEIGHTING
    #compute Boltzmann weights (reweighting)
    EN_FENE_REW,  STCK_MOD_REW, HYDR_MOD_REW, CRST_33_MOD_REW, CRST_55_MOD_REW = compute_rew_factor(CURR_PARS,SH_ST)
    D1 = EN_FENE_IN - EN_FENE_REW + EN_STCK_IN-STCK_MOD_FIX*STCK_MOD_REW
    D2 = EN_HYDR_IN-HYDR_MOD_FIX*HYDR_MOD_REW
    D3 = EN_CRST_33_IN-CRST_33_MOD_FIX*CRST_33_MOD_REW
    D4 = EN_CRST_55_IN-CRST_55_MOD_FIX*CRST_55_MOD_REW

    #D2 = torch.where(HYDR_MOD_IN > 0, EN_HYDR_IN*(1-HYDR_MOD_REW/(HYDR_MOD_IN)), ZEROS_UNBN)
    #D3 = torch.where(CRST_33_MOD_IN > 0, EN_CRST_33_IN*(1-CRST_33_MOD_REW/(CRST_33_MOD_IN+1e-12)), ZEROS_UNBN)
    #D4 = torch.where(CRST_55_MOD_IN > 0, EN_CRST_55_IN*(1-CRST_55_MOD_REW/(CRST_55_MOD_IN+1e-12)), ZEROS_UNBN)


    WEIGHTS = torch.exp( 10*torch.sum( D1,dim=2 ) + 10*torch.sum( D2+D3+D4,dim=2 ) )

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

    #global MU_RED_TARGET_MU
    #print(REW_MU)
    #print(MU_RED_TARGET_MU)
    #AVERAGE PROPELLER!
    if 1 in cg.ids_gs:
        AVE_REW_PROP = torch.mean(REW_MU[:,PROP_IDS])
        REW_MU[:,PROP_IDS] = AVE_REW_PROP

        global ave_prop_flag
        if ave_prop_flag :
            global TARGET_MU
            global MU_RED_TARGET_MU
            AVE_TARGET_PROP = torch.mean(TARGET_MU[:,PROP_IDS])
            TARGET_MU[:,PROP_IDS] = AVE_TARGET_PROP
            MU_RED_TARGET_MU = TARGET_MU[:,MU_RED_IDS]
            ave_prop_flag = False

    #print(REW_MU)
    #print(MU_RED_TARGET_MU)

    #reduce reweighted covariances
    TMP = REW_COV[:,:,COV_RED_IDS]
    COV_RED_REW_COV = TMP[:,COV_RED_IDS]
    TMP = REW_COV[:,:,MU_RED_IDS]
    MU_RED_REW_COV = TMP[:,MU_RED_IDS]
    MU_RED_REW_MU = REW_MU[:,MU_RED_IDS]

    DELTA = MU_RED_REW_MU - MU_RED_TARGET_MU

    #GROUND STATE
    S = 0.5*torch.einsum('ij,ij->i',DELTA, torch.einsum('ijk,ij->ik',MU_RED_TARGET_COV_m1, DELTA))

    #If we optimise the covariance, we add this bit to the cost of the GS

    DIA_MRRC = torch.diagonal(MU_RED_REW_COV, 0, 1, 2)
    DIA_MRRC_m1 = 1/DIA_MRRC
    #MRRC_m1 = torch.linalg.inv(MU_RED_REW_COV)

    S += 0.5*torch.einsum('ij,ij->i',DELTA, DIA_MRRC_m1*DELTA)

    S = S.sum(dim=0)

    #COVARIANCE

    #NOTE: we are optimising the STD only (i.e. the diagonal of the covariance matrix).
    #This is why we use the diagonal only


    DIA_CRTC = torch.diagonal(COV_RED_TARGET_COV, 0, 1, 2)
    DIA_CRTC_m1 = torch.diagonal(COV_RED_TARGET_COV_m1, 0, 1, 2)
    DIA_CRRC = torch.diagonal(COV_RED_REW_COV, 0, 1, 2)
    #DIA_CRRC_m1 = torch.diagonal(torch.linalg.inv(COV_RED_REW_COV), 0, 1, 2)

    #offset covariance so that Target and reweighted have the same average. In this way we are optimising the shape and not the average (i.e. persistence length)
    TMP = DIA_CRRC[:,IDS_AVE_COV_SUM].mean(dim=2)

    AVE_DIA_CRRC = TMP[:,IDS_AVE_COV_EXPAND]

    DIA_CRRC = DIA_CRRC - AVE_DIA_CRRC + AVE_COV_RED_TARGET_COV

    #print("offset: ", - AVE_DIA_CRRC + AVE_COV_RED_TARGET_COV)
    #print("dia crrc: ", DIA_CRRC)

    DIA_CRRC_m1 = 1/(DIA_CRRC+1e-12)

    #S += 0.5*(torch.einsum('ij,ij->i',DIA_CRTC_m1,DIA_CRRC)+torch.einsum('ij,ij->i',DIA_CRRC_m1,DIA_CRTC)-2).sum(dim=0)

    S += 0.5*(torch.square(DIA_CRRC-DIA_CRTC)/DIA_CRRC/DIA_CRTC).sum(dim=1).sum(dim=0)

    #print("S2: ", torch.einsum('ij,ij->i',DIA_CRTC_m1,DIA_CRRC)+torch.einsum('ij,ij->i',DIA_CRRC_m1,DIA_CRTC))
    #print("S2 - 2: ", torch.einsum('ij,ij->i',DIA_CRTC_m1,DIA_CRRC)+torch.einsum('ij,ij->i',DIA_CRRC_m1,DIA_CRTC)-2)
    #print("S2 - 2 v2: ", 0.5*(torch.square(DIA_CRRC-DIA_CRTC)/DIA_CRRC/DIA_CRTC).sum(dim=1))

    #print("S: "+str(float(S)))

    #diagonal(...) computes the trace of a batch of matrices
    #S += 0.5*( torch.tensordot(COV_RED_TARGET_COV_m1,COV_RED_REW_COV, dims = ([1,2])).diagonal(offset=0, dim1=-1, dim2=-2).sum(-1) ).sum(dim=0)
    #S += 0.5*( torch.tensordot(torch.linalg.inv_ex(COV_RED_REW_COV),COV_RED_TARGET_COV, dims = ([1,2])).diagonal(offset=0, dim1=-1, dim2=-2).sum(-1) ).sum(dim=0)


    #PERSISTENCE LENGTH

    global LB
    global LT

    if cg.opti_lp :

        """
        TMP = REW_COV[:,:,LP_RED_IDS]
        LP_RED_REW_COV = TMP[:,LP_RED_IDS]

        K = compute_stiffness(LP_RED_REW_COV)

        A1 = K[:,0,0]
        A2 = K[:,1,1]
        C = K[:,2,2]
        G = K[:,1,2]

        lb = 2*A1*(A2-G*G/C)/(A1+A2-G*G/C)
        lt = 2*C*(1-G*G/A2/C)

        LT = lt
        LB = lb

        s0, s1, s2 = LP_RED_REW_COV.size()
        #C = C.sum(dim=0)
        #A2 = A2.sum(dim=0)
        """

        #reweight costb and cosot

        COS_THETAB_REW = WEIGHTS*COS_THETAB
        COS_OMEGAT_REW = WEIGHTS*COS_OMEGAT

        ave_rise = SAVE_MU[:,cg.ids.index(11)]/10.

        lb = -cg.lp_m*ave_rise/torch.log(COS_THETAB_REW.sum(dim=1)) *cg.lb_corr
        lt = -cg.lp_m*ave_rise/torch.log(COS_OMEGAT_REW.sum(dim=1)) *cg.lt_corr

        #print(lb)
        #print(lt)

        LT = lt.clone()
        LB = lb.clone()

        lt = lt.sum(dim=0)/cg.Nseq
        lb = lb.sum(dim=0)/cg.Nseq

        N = cg.fin_j[0] - cg.in_j[0] + 1

        #S += 0.5*cg.Nseq*(C/140+140/C+40/A2+A2/40-4) #*(s1/3)
        #S += 0.5*cg.Nseq*((lt-210)*(lt-210)/lt/210+(lb-45)*(lb-45)/45/lb)*N
        S += 0.5*cg.Nseq*((lb-45)*(lb-45)/45/lb)*N

    Scpu = float(S)

    return Scpu
