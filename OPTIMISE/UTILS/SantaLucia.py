#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:19:00 2024

@author: yqb22156
"""

import numpy as np
import math

def base_to_id(base) :
    
    if base == 'A':
        return  0
    elif base == 'G':
        return 1
    elif base == 'C':
        return 2
    elif base == 'T':
        return 3
    
    return -1


def id_to_base(base) :
    
    if base == 0:
        return  'A'
    elif base == 1:
        return 'G'
    elif base == 2:
        return 'C'
    elif base == 3:
        return 'T'
    
    return 'Z'

DH0in = 0.2
DS0in = -5.7

DH0fin_AT = 2.2
DS0fin_AT = 6.9

DS0Symm = -1.4

R = 1.9872

def self_complementary(seq_num) :    
    
    ls = len(seq_num)
    
    #print(seq_num)
    
    if ls%2 == 1:
        return False
    
    for i in range(int(ls/2)) :
        if (seq_num[i]+seq_num[ls-1-i]) != 3:
            return False    
    
    return True


DH0_step = np.array( [ [-7.6,-7.8,-8.4,-7.2], 
                  [-8.2,-8.0,-9.8,-8.4],        
                  [-8.5,-10.6,-8.0,-7.8],        
                  [-7.2,-8.5,-8.2,-7.6] ] )


DS0_step = np.array( [ [-21.3,-21.0,-22.4,-20.4],
                  [-22.2,-19.9,-24.4,-22.4],
                  [-22.7,-27.2,-19.9,-21.0],
                  [-21.3,-22.7,-22.2,-21.3] ] )


DH0in_ave = 0.2
DS0in_ave = -5.7

DH0_step_ave = 0
DS0_step_ave = 0

for i in range(len(DH0_step)) :
    for j in range(len(DH0_step[i])) :
        DH0_step_ave += DH0_step[i][j]/16.
        DS0_step_ave += DS0_step[i][j]/16.
        
DH0fin_ave = 2.2
DS0fin_ave = 6.9


def melting_temperature(seq,Ct,Cs) :
    
    #Santa Lucia (and cgna) reads 5'->3', oxdna reads 3'->5'.
    #To get the right sequence, must reverse oxdna seq!
    SL_seq = "".join(reversed(seq))
    #SL_seq = seq
    #print(SL_seq)
    
    seq_num = []    #from letters to numbers
    
    for i in range(len(SL_seq)) :
        seq_num.append(base_to_id(SL_seq[i]))   
    
    
    DS0 = 0.
    DH0 = 0.
    
    DS0 += DS0in
    DH0 += DH0in
    for i in range(len(seq_num)-1):
        DS0 += DS0_step[seq_num[i]][seq_num[i+1]]
        DH0 += DH0_step[seq_num[i]][seq_num[i+1]]
        
    if seq_num[0] == 0 or seq_num[0] == 3 :
        #print("In AT")
        DS0 += DS0fin_AT
        DH0 += DH0fin_AT
        
    if seq_num[len(seq_num)-1] == 0 or seq_num[len(seq_num)-1] == 3 :
        #print("Fin AT")
        DS0 += DS0fin_AT
        DH0 += DH0fin_AT
        
        
    #apply salt correction
    
    DS0 = DS0 + 0.368*(len(seq_num)-1)*math.log(Cs)
    
    #print(DS0)
    #print(DH0)
    
    if self_complementary(seq_num) :
        #print("complementary")
        DS0 += DS0Symm
        
    x = 4.
    
    if self_complementary(seq_num) :
        x = 1.
    
    Tm = DH0*1000/(DS0 + R*math.log(Ct/x))-273.15
    
    return Tm






def melting_temperature_ave(seq,Ct,Cs) :
        
    #Santa Lucia (and cgna) reads 5'->3', oxdna reads 3'->5'.
    #To get the right sequence, must reverse oxdna seq!
    SL_seq = "".join(reversed(seq))
    #SL_seq = seq
    #print(SL_seq)
    
    seq_num = []    #from letters to numbers
    
    for i in range(len(SL_seq)) :
        seq_num.append(base_to_id(SL_seq[i]))   
    
    
    DS0 = 0.
    DH0 = 0.
    
    DS0 += DS0in_ave + (len(seq_num)-1)*DS0_step_ave + DS0fin_ave
    DH0 += DH0in_ave + (len(seq_num)-1)*DH0_step_ave + DH0fin_ave
 
        
    #apply salt correction
    
    if len(seq_num)%2 == 0 :
        DS0 += DS0Symm*pow(4,-len(seq_num)/2)
    
    
    DS0 = DS0 + 0.368*(len(seq_num)-1)*math.log(Cs)
    
    #print(DS0)
    #print(DH0)
    
    """
    
    #lower order correction. Can be neglected
    
    if self_complementary(seq_num) :
        #print("complementary")
        DS0 += DS0Symm
        
    if self_complementary(seq_num) :
        x = 1.
    """
    x = 4.
    
    Tm = DH0*1000/(DS0 + R*math.log(Ct/x))-273.15
    
    return Tm