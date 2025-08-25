#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 11:22:36 2023

@author: yqb22156
"""

import sys, random
import numpy as np
from get_cgdna_pars import constructSeqParms
from get_cgdna_pars import reduce



seq1 = "ACACACACACACACACACAC"  # randomly chosen sequence


gs,stiff = constructSeqParms(seq1)

Inv_cov = stiff.A
cov = np.linalg.inv(Inv_cov)

#0,1,2,3,4,5 phosphate
#6,7,8,9,10,11 intra
#12,13,14,15,16,17 phosphate
#18,19,20,21,22,23 inter

ids = [6,7,8,9,10,11,18,19,20,21,22,23]

jx_in = 1
jx_from_fin = 1

nbp = int((len(Inv_cov)+18)/24)
nj = nbp-1

N = int(len(ids)*((nj)-jx_in-jx_from_fin))

ave_red_cov = np.zeros((N,N),dtype=float)
ave_red_gs = np.zeros(N,dtype=float)

#produce average covariance

Nseq = 10000
Lseq = 20

for i in range(Nseq) :
    
    print("seq " + str(i))
    seq = ""
    for j in range(Lseq):
        
        x = random.randint(0, 3)        
        
        if x == 0:
            seq += 'A'  # randomly chosen sequence
        if x == 1:
            seq += 'C'  # randomly chosen sequence
        if x == 2:
            seq += 'G'  # randomly chosen sequence
        if x == 3:
            seq += 'T'  # randomly chosen sequence

    print(seq)

    gs,stiff = constructSeqParms(seq)

    Inv_cov = stiff.A
    
    red_gs, red_stiff, red_cov = reduce(gs,Inv_cov,ids,jx_in,jx_from_fin)
    
    ave_red_cov += red_cov*(1./Nseq)
    ave_red_gs += red_gs*(1./Nseq)
 


ave_SEQ_gs = np.load("cgna_SEQave_gs_20bp_17jx_noends.npy")
ave_SEQ_cov = np.load("cgna_SEQave_cov_20bp_17jx_noends.npy")

ave_gs = np.zeros(12,dtype=float)
ave_cov = np.zeros((12,12),dtype=float)


jx_in = 6
jx_fin = 12


counts_i = 0
for i in range(jx_in,jx_fin) :
    print(i)
    for n in range(12) :
        ave_gs[n] += ave_SEQ_gs[i*12+n]/(jx_fin-jx_in)
        for m in range(12) :
            ave_cov[n][m] += ave_SEQ_cov[i*12+n][i*12+m]/(jx_fin-jx_in)
                              