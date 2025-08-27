#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 14:50:24 2024

@author: yqb22156
"""

#list of oxdna parameters used in energy computation
PARS_LIST = [
            'FENE_EPS', 'FENE_R0', 'FENE_DELTA', 'FENE_DELTA2',
             'HYDR_EPS', 'HYDR_R0', 'HYDR_A', 'HYDR_RC', 'HYDR_BLOW', 'HYDR_BHIGH', 'HYDR_RLOW', 'HYDR_RHIGH', 'HYDR_RCLOW', 'HYDR_RCHIGH',
             'HYDR_THETA1_T0', 'HYDR_THETA1_A', 'HYDR_THETA1_B', 'HYDR_THETA1_TS', 'HYDR_THETA1_TC',
             'HYDR_THETA2_T0', 'HYDR_THETA2_A', 'HYDR_THETA2_B', 'HYDR_THETA2_TS', 'HYDR_THETA2_TC',
             'HYDR_THETA3_T0', 'HYDR_THETA3_A', 'HYDR_THETA3_B', 'HYDR_THETA3_TS', 'HYDR_THETA3_TC',
             'HYDR_THETA4_T0', 'HYDR_THETA4_A', 'HYDR_THETA4_B', 'HYDR_THETA4_TS', 'HYDR_THETA4_TC',
             'HYDR_THETA7_T0', 'HYDR_THETA7_A', 'HYDR_THETA7_B', 'HYDR_THETA7_TS', 'HYDR_THETA7_TC',
             'HYDR_THETA8_T0', 'HYDR_THETA8_A', 'HYDR_THETA8_B', 'HYDR_THETA8_TS', 'HYDR_THETA8_TC',
             'STCK_EPS', 'STCK_R0', 'STCK_A', 'STCK_RC', 'STCK_BLOW', 'STCK_BHIGH', 'STCK_RLOW', 'STCK_RHIGH', 'STCK_RCLOW', 'STCK_RCHIGH',
             'STCK_THETA4_T0', 'STCK_THETA4_A', 'STCK_THETA4_B', 'STCK_THETA4_TS', 'STCK_THETA4_TC',
             'STCK_THETA5_T0', 'STCK_THETA5_A', 'STCK_THETA5_B', 'STCK_THETA5_TS', 'STCK_THETA5_TC',
             'STCK_THETA6_T0', 'STCK_THETA6_A', 'STCK_THETA6_B', 'STCK_THETA6_TS', 'STCK_THETA6_TC',
             'STCK_PHI1_A', 'STCK_PHI1_B', 'STCK_PHI1_XC', 'STCK_PHI1_XS',
             'STCK_PHI2_A', 'STCK_PHI2_B', 'STCK_PHI2_XC', 'STCK_PHI2_XS',
             'CRST_K_33', 'CRST_R0_33', 'CRST_RC_33', 'CRST_BLOW_33', 'CRST_BHIGH_33', 'CRST_RLOW_33', 'CRST_RHIGH_33', 'CRST_RCLOW_33', 'CRST_RCHIGH_33',
             'CRST_THETA1_T0_33', 'CRST_THETA1_A_33', 'CRST_THETA1_B_33', 'CRST_THETA1_TS_33', 'CRST_THETA1_TC_33',
             'CRST_THETA2_T0_33', 'CRST_THETA2_A_33', 'CRST_THETA2_B_33', 'CRST_THETA2_TS_33', 'CRST_THETA2_TC_33',
             'CRST_THETA3_T0_33', 'CRST_THETA3_A_33', 'CRST_THETA3_B_33', 'CRST_THETA3_TS_33', 'CRST_THETA3_TC_33',
             'CRST_THETA4_T0_33', 'CRST_THETA4_A_33', 'CRST_THETA4_B_33', 'CRST_THETA4_TS_33', 'CRST_THETA4_TC_33',
             'CRST_THETA7_T0_33', 'CRST_THETA7_A_33', 'CRST_THETA7_B_33', 'CRST_THETA7_TS_33', 'CRST_THETA7_TC_33',
             'CRST_THETA8_T0_33', 'CRST_THETA8_A_33', 'CRST_THETA8_B_33', 'CRST_THETA8_TS_33', 'CRST_THETA8_TC_33',
             'CRST_K_55', 'CRST_R0_55', 'CRST_RC_55', 'CRST_BLOW_55', 'CRST_BHIGH_55', 'CRST_RLOW_55', 'CRST_RHIGH_55', 'CRST_RCLOW_55', 'CRST_RCHIGH_55',
             'CRST_THETA1_T0_55', 'CRST_THETA1_A_55', 'CRST_THETA1_B_55', 'CRST_THETA1_TS_55', 'CRST_THETA1_TC_55',
             'CRST_THETA2_T0_55', 'CRST_THETA2_A_55', 'CRST_THETA2_B_55', 'CRST_THETA2_TS_55', 'CRST_THETA2_TC_55',
             'CRST_THETA3_T0_55', 'CRST_THETA3_A_55', 'CRST_THETA3_B_55', 'CRST_THETA3_TS_55', 'CRST_THETA3_TC_55',
             'CRST_THETA4_T0_55', 'CRST_THETA4_A_55', 'CRST_THETA4_B_55', 'CRST_THETA4_TS_55', 'CRST_THETA4_TC_55',
             'CRST_THETA7_T0_55', 'CRST_THETA7_A_55', 'CRST_THETA7_B_55', 'CRST_THETA7_TS_55', 'CRST_THETA7_TC_55',
             'CRST_THETA8_T0_55', 'CRST_THETA8_A_55', 'CRST_THETA8_B_55', 'CRST_THETA8_TS_55', 'CRST_THETA8_TC_55',
             'EXCL_S1', 'EXCL_R1', 'EXCL_B1', 'EXCL_RC1',
             'EXCL_S2', 'EXCL_R2', 'EXCL_B2', 'EXCL_RC2',
             'EXCL_S3', 'EXCL_R3', 'EXCL_B3', 'EXCL_RC3',
             'EXCL_S4', 'EXCL_R4', 'EXCL_B4', 'EXCL_RC4',
             'EXCL_S5', 'EXCL_R5', 'EXCL_B5', 'EXCL_RC5',
             'EXCL_S6', 'EXCL_R6', 'EXCL_B6', 'EXCL_RC6',
             'EXCL_S7', 'EXCL_R7', 'EXCL_B7', 'EXCL_RC7'
             ]



#create list with indices of parameters used in compute_energy.
#The purpuse of this is that even if PARS_LIST gets changed, compute_energy doesn't need to be updated. This should be faster than searching within a map
#ids are reported below
par_index = []


#FENE

par_index.append(PARS_LIST.index("FENE_EPS"))
par_index.append(PARS_LIST.index("FENE_R0"))
par_index.append(PARS_LIST.index("FENE_DELTA"))
par_index.append(PARS_LIST.index("FENE_DELTA2"))


#HYDROGEN

par_index.append(PARS_LIST.index("HYDR_EPS"))

par_index.append(PARS_LIST.index("HYDR_R0"))
par_index.append(PARS_LIST.index("HYDR_A"))
par_index.append(PARS_LIST.index("HYDR_RC"))
par_index.append(PARS_LIST.index("HYDR_BLOW"))
par_index.append(PARS_LIST.index("HYDR_BHIGH"))
par_index.append(PARS_LIST.index("HYDR_RLOW"))
par_index.append(PARS_LIST.index("HYDR_RHIGH"))
par_index.append(PARS_LIST.index("HYDR_RCLOW"))
par_index.append(PARS_LIST.index("HYDR_RCHIGH"))

par_index.append(PARS_LIST.index("HYDR_THETA1_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA1_A"))
par_index.append(PARS_LIST.index("HYDR_THETA1_B"))
par_index.append(PARS_LIST.index("HYDR_THETA1_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA1_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA2_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA2_A"))
par_index.append(PARS_LIST.index("HYDR_THETA2_B"))
par_index.append(PARS_LIST.index("HYDR_THETA2_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA2_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA3_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA3_A"))
par_index.append(PARS_LIST.index("HYDR_THETA3_B"))
par_index.append(PARS_LIST.index("HYDR_THETA3_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA3_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA4_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA4_A"))
par_index.append(PARS_LIST.index("HYDR_THETA4_B"))
par_index.append(PARS_LIST.index("HYDR_THETA4_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA4_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA7_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA7_A"))
par_index.append(PARS_LIST.index("HYDR_THETA7_B"))
par_index.append(PARS_LIST.index("HYDR_THETA7_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA7_TC"))

par_index.append(PARS_LIST.index("HYDR_THETA8_T0"))
par_index.append(PARS_LIST.index("HYDR_THETA8_A"))
par_index.append(PARS_LIST.index("HYDR_THETA8_B"))
par_index.append(PARS_LIST.index("HYDR_THETA8_TS"))
par_index.append(PARS_LIST.index("HYDR_THETA8_TC"))


#STACKING

par_index.append(PARS_LIST.index("STCK_EPS"))

par_index.append(PARS_LIST.index("STCK_R0"))
par_index.append(PARS_LIST.index("STCK_A"))
par_index.append(PARS_LIST.index("STCK_RC"))
par_index.append(PARS_LIST.index("STCK_BLOW"))
par_index.append(PARS_LIST.index("STCK_BHIGH"))
par_index.append(PARS_LIST.index("STCK_RLOW"))
par_index.append(PARS_LIST.index("STCK_RHIGH"))
par_index.append(PARS_LIST.index("STCK_RCLOW"))
par_index.append(PARS_LIST.index("STCK_RCHIGH"))

par_index.append(PARS_LIST.index("STCK_THETA4_T0"))
par_index.append(PARS_LIST.index("STCK_THETA4_A"))
par_index.append(PARS_LIST.index("STCK_THETA4_B"))
par_index.append(PARS_LIST.index("STCK_THETA4_TS"))
par_index.append(PARS_LIST.index("STCK_THETA4_TC"))

par_index.append(PARS_LIST.index("STCK_THETA5_T0"))
par_index.append(PARS_LIST.index("STCK_THETA5_A"))
par_index.append(PARS_LIST.index("STCK_THETA5_B"))
par_index.append(PARS_LIST.index("STCK_THETA5_TS"))
par_index.append(PARS_LIST.index("STCK_THETA5_TC"))

par_index.append(PARS_LIST.index("STCK_THETA6_T0"))
par_index.append(PARS_LIST.index("STCK_THETA6_A"))
par_index.append(PARS_LIST.index("STCK_THETA6_B"))
par_index.append(PARS_LIST.index("STCK_THETA6_TS"))
par_index.append(PARS_LIST.index("STCK_THETA6_TC"))

par_index.append(PARS_LIST.index("STCK_PHI1_A"))
par_index.append(PARS_LIST.index("STCK_PHI1_B"))
par_index.append(PARS_LIST.index("STCK_PHI1_XC"))
par_index.append(PARS_LIST.index("STCK_PHI1_XS"))

par_index.append(PARS_LIST.index("STCK_PHI2_A"))
par_index.append(PARS_LIST.index("STCK_PHI2_B"))
par_index.append(PARS_LIST.index("STCK_PHI2_XC"))
par_index.append(PARS_LIST.index("STCK_PHI2_XS"))


#CROSS STACKING

par_index.append(PARS_LIST.index("CRST_K_33"))
par_index.append(PARS_LIST.index("CRST_R0_33"))
par_index.append(PARS_LIST.index("CRST_RC_33"))
par_index.append(PARS_LIST.index("CRST_BLOW_33"))
par_index.append(PARS_LIST.index("CRST_BHIGH_33"))
par_index.append(PARS_LIST.index("CRST_RLOW_33"))
par_index.append(PARS_LIST.index("CRST_RHIGH_33"))
par_index.append(PARS_LIST.index("CRST_RCLOW_33"))
par_index.append(PARS_LIST.index("CRST_RCHIGH_33"))

par_index.append(PARS_LIST.index("CRST_THETA1_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA1_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA2_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA2_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA3_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA3_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA4_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA4_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA7_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA7_TC_33"))

par_index.append(PARS_LIST.index("CRST_THETA8_T0_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_A_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_B_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_TS_33"))
par_index.append(PARS_LIST.index("CRST_THETA8_TC_33"))

par_index.append(PARS_LIST.index("CRST_K_55"))
par_index.append(PARS_LIST.index("CRST_R0_55"))
par_index.append(PARS_LIST.index("CRST_RC_55"))
par_index.append(PARS_LIST.index("CRST_BLOW_55"))
par_index.append(PARS_LIST.index("CRST_BHIGH_55"))
par_index.append(PARS_LIST.index("CRST_RLOW_55"))
par_index.append(PARS_LIST.index("CRST_RHIGH_55"))
par_index.append(PARS_LIST.index("CRST_RCLOW_55"))
par_index.append(PARS_LIST.index("CRST_RCHIGH_55"))

par_index.append(PARS_LIST.index("CRST_THETA1_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA1_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA2_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA2_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA3_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA3_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA4_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA4_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA7_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA7_TC_55"))

par_index.append(PARS_LIST.index("CRST_THETA8_T0_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_A_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_B_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_TS_55"))
par_index.append(PARS_LIST.index("CRST_THETA8_TC_55"))

#EXCLUDED VOLUME
par_index.append(PARS_LIST.index("EXCL_S1"))
par_index.append(PARS_LIST.index("EXCL_R1"))
par_index.append(PARS_LIST.index("EXCL_B1"))
par_index.append(PARS_LIST.index("EXCL_RC1"))

par_index.append(PARS_LIST.index("EXCL_S2"))
par_index.append(PARS_LIST.index("EXCL_R2"))
par_index.append(PARS_LIST.index("EXCL_B2"))
par_index.append(PARS_LIST.index("EXCL_RC2"))

par_index.append(PARS_LIST.index("EXCL_S3"))
par_index.append(PARS_LIST.index("EXCL_R3"))
par_index.append(PARS_LIST.index("EXCL_B3"))
par_index.append(PARS_LIST.index("EXCL_RC3"))

par_index.append(PARS_LIST.index("EXCL_S4"))
par_index.append(PARS_LIST.index("EXCL_R4"))
par_index.append(PARS_LIST.index("EXCL_B4"))
par_index.append(PARS_LIST.index("EXCL_RC4"))

par_index.append(PARS_LIST.index("EXCL_S5"))
par_index.append(PARS_LIST.index("EXCL_R5"))
par_index.append(PARS_LIST.index("EXCL_B5"))
par_index.append(PARS_LIST.index("EXCL_RC5"))

par_index.append(PARS_LIST.index("EXCL_S6"))
par_index.append(PARS_LIST.index("EXCL_R6"))
par_index.append(PARS_LIST.index("EXCL_B6"))
par_index.append(PARS_LIST.index("EXCL_RC6"))

par_index.append(PARS_LIST.index("EXCL_S7"))
par_index.append(PARS_LIST.index("EXCL_R7"))
par_index.append(PARS_LIST.index("EXCL_B7"))
par_index.append(PARS_LIST.index("EXCL_RC7"))

compl_symm_ids_2d = [0]
compl_symm_ids_2d_no_TT_AA = []
compl_symm_ids = [1,2]
"""
for i in range(45,77) :
    compl_symm_ids.append(i)
compl_symm_ids_no_TT_AA = []
perm_symm_ids_4d = []
"""
perm_symm_ids_2d = [77,116]

for i in range(4,44) :
    perm_symm_ids_2d.append(i)

perm_symm_ids_4d = []

for i in range(78,116) :
    perm_symm_ids_4d.append(i)

for i in range(117,155) :
    perm_symm_ids_4d.append(i)

#nosymm_ids_2d = [44]
nosymm_ids_2d = []
compl_symm_ids_no_TT_AA = []

is_th2 = [19,20,21,22,23,91,92,93,94,95,130,131,132,133,134]
is_th3 = [24,25,26,27,28,96,97,98,99,100,135,136,137,138,139]
is_th5 = [59,60,61,62,63]
is_th6 = [64,65,66,67,68]
is_th7 = [34,35,36,37,38,106,107,108,109,110,145,146,147,148,149]
is_th8 = [39,40,41,42,43,111,112,113,114,115,150,151,152,153,154]

