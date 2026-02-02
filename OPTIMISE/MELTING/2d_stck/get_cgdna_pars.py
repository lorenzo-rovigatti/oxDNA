#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:35:31 2023

@author: andrea bonato - constructSeqParms function and parameters (cgDNA+ps1.mat) copyed from rahul sharma (cgna+)
"""

import sys
import numpy as np
import scipy, scipy.io
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
import os
#path = os.getcwd()
path=os.path.dirname(os.path.realpath(__file__))

par_file = path+'/Param/cgDNA+ps1.mat'
path_averages = path+'/Param/'


def finder(seq):
	istart = []  
	end = {}
	start = []
	for i, c in enumerate(seq):
		if c == '[':
			istart.append(i)
			start.append(i)
		if c == ']':
			try:
				end[istart.pop()] = i
			except IndexError:
				print('Too many closing parentheses')
	if istart:  # check if stack is empty afterwards
		print('Too many opening parentheses')
	return end, start

def mult(seq):
	i =seq.rfind('_') 
	if seq[i+1].isdigit():
		a = seq[i+1]
		if seq[i+2].isdigit():
			a = a + seq[i+2]
			if seq[i+3].isdigit():
				a = a + seq[i+3]
				if seq[i+4].isdigit():
					a = a + seq[i+4]
					if seq[i+5].isdigit():
						a = a + seq[i+5]


def seq_edit(seq):
	s = seq.upper()
	while s.rfind('_')>0:
		if s[s.rfind('_')-1].isdigit():
			print("Please write the input sequence correctly. Two or more _ can't be put consequently. You can use the brackets. i.e. A_2_2 can be written as [A_2]_2")
			exit()
		if s[s.rfind('_')-1] != ']':
			a = int(mult(s))
			s = s[:s.rfind('_')-1]+ s[s.rfind('_')-1]*a +  s[s.rfind('_')+1+len(str((a))):]
		if s[s.rfind('_')-1] == ']':
			end,start = finder(s)
			ka=(2,len(start))
			h=np.zeros(ka)
			for i in range(len(start)):
				h[0][i] = start[i]
				h[1][i] = end[start[i]]	
			ss=  int(max(h[1]))
			ee=  int(h[0][np.argmax(h[1])])
			a = int(mult(s))
			s =  s[0:ee] + s[ee+1:ss]*a + s[ss+2+len(str((a))):] 
	return s






def constructSeqParms(sequence):

    
    ps = scipy.io.loadmat(par_file)

	#### Following loop take every input sequence and construct shape and stiff matrix ###
    s_seq = seq_edit(sequence)
    nbp = len(s_seq.strip())
    N = 24*nbp-18

	#### Initialise the sigma vector ###		
    s = np.zeros((N,1))

    #### Error report if sequence provided is less than 2 bp #### 

    if nbp <= 3:
        print("sequence length must be greater than or equal to 2")
        sys.exit() 

    elif nbp > 3 :
    
        data,row,col = {},{},{}
        
        ### 5' end #### 
        tmp_ind = np.nonzero(ps['stiff_end5'][s_seq[0:2]][0][0][0:36,0:36])
        row[0],col[0] = tmp_ind[0][:],tmp_ind[1][:]
        data[0] = ps['stiff_end5'][s_seq[0:2]][0][0][row[0],col[0]]
        
        s[0:36] = ps['sigma_end5'][s_seq[0:2]][0][0][0:36]
        #### interior blocks  ###
        for i in range(2,nbp-1):
            tmp_ind = np.nonzero(ps['stiff_int'][s_seq[i-1:i+1]][0][0][0:42, 0:42])
            data[i-1] = ps['stiff_int'][s_seq[i-1:i+1]][0][0][tmp_ind[0][:], tmp_ind[1][:]]
            
            di = 24*(i-2)+18
            row[i-1] = tmp_ind[0][:]+np.ones((1,np.size(tmp_ind[0][:])))*di
            col[i-1] = tmp_ind[1][:]+np.ones((1,np.size(tmp_ind[1][:])))*di
            
            s[di:di+42] = np.add(s[di:di+42],ps['sigma_int'][s_seq[i-1:i+1]][0][0][0:42])
			
		#### 3' end ####
        tmp_ind = np.nonzero(ps['stiff_end3'][s_seq[nbp-2:nbp]][0][0][0:36, 0:36])
        data[nbp-1] = ps['stiff_end3'][s_seq[nbp-2:nbp]][0][0][tmp_ind[0][:], tmp_ind[1][:]]
        
        di = 24*(nbp-3)+18
        row[nbp-1] = tmp_ind[0][:]+np.ones((1,np.size(tmp_ind[0][:])))*di
        col[nbp-1] = tmp_ind[1][:]+np.ones((1,np.size(tmp_ind[1][:])))*di
        s[N-36:N] = s[N-36:N] + ps['sigma_end3'][s_seq[nbp-2:nbp]][0][0][0:36]
       
        tmp = list(row.values())
        row = np.concatenate(tmp,axis=None)
        
        tmp = list(col.values())
        col = np.concatenate(tmp,axis=None)
   
        tmp = list(data.values())
        data = np.concatenate(tmp,axis=None)
        
    
    #### Create the sparse Stiffness matrix from data,row_ind,col_ind  ###
        stiff =  csc_matrix((data, (row,col)), shape =(N,N))	

	#### Groudstate calculation ####
        ground_state = spsolve(stiff, s) 

    return ground_state,stiff


def reduce(gs,stiff,ids,jx_in,jx_from_fin) :
    
    nbp = int((len(stiff)+18)/24)
    nj = nbp-1
    
    N = int(len(ids)*((nj)-jx_in-jx_from_fin))
    
    index_in = 0
    if jx_in >= 1 :
        index_in+=18
        if jx_in > 1 :
            index_in+=24*(jx_in-1)
    
    index_fin = len(stiff)
    if jx_from_fin >= 1 :
        index_fin-=36
        if jx_in > 1 :
            index_fin-=24*(jx_from_fin-1)
            
    print(index_in)
    print(index_fin)
    
    cov = np.linalg.inv(stiff) #invert stiffness, reduce covariance (multivaraiate dist)     
    red_cov = np.zeros((N,N),dtype=float)
    red_gs = np.zeros(N,dtype=float)
    
    nrow = -1
    ncol = -1
    
    for l in range(index_in,index_fin) :

        if (l-index_in)%24 in ids :
            nrow += 1  
            #print("norw: "+str(nrow))              
            ncol = -1
            red_gs[nrow] = gs[l]
            for m in range(index_in,index_fin) :
                if (m-index_in)%24 in ids :
                    
                    ncol+=1
                    #print("m: "+str(m))
                    #print("ncol: "+str(ncol))
                    red_cov[nrow,ncol]=cov[l,m]
            
    
    
    red_stiff = np.linalg.inv(red_cov) 
    
    return red_gs, red_stiff, red_cov  


def get_target_mean_and_covariance(seq, ids, jx_in, jx_from_end) :
    #map optimisation ids to cgna+ids (i.e.: remove phosphate)
    cgDNA_ids = []
    for i in range(len(ids)) :
        if ids[i] >= 0 and ids[i] <= 5: 
            cgDNA_ids.append(ids[i]+6)
        else :
            cgDNA_ids.append(ids[i]+12)
            
            
    #cgna reads 5'->3', oxdna reads 3'->5'.
    #To get the right sequence, must reverse oxdna seq!
    cgna_seq = "".join(reversed(seq))            
            
    gs,stiff = constructSeqParms(cgna_seq)
    Inv_cov = stiff.A
    red_gs, red_stiff, red_cov = reduce(gs,Inv_cov,cgDNA_ids,jx_in,jx_from_end)
    
    return red_gs.tolist(), red_cov.tolist()
    

def get_target_mean_and_covariance_diag_ave(N, ids) :

    gs = np.load(path_averages+"cgna_ave_gs.npy")
    cov = np.load(path_averages+"cgna_ave_cov.npy")
    
    
    red_gs = np.zeros(len(ids),dtype = float)
    red_cov = np.zeros((len(ids),len(ids)),dtype = float)
    
    #reduce gs and cov

    nrow = -1
    
    for i in range(len(gs)):        
        if i%12 in ids :
            nrow+=1
            red_gs[nrow] = gs[i]
            ncol = -1
            for j in range(len(gs)):
                if j%12 in ids:
                    ncol +=1
                    red_cov[nrow,ncol]=cov[i,j]
                    
    
    #build Nbps cov and gs
    
    tot_gs = np.zeros(N*len(red_gs),dtype=float)
    tot_cov = np.zeros((N*len(red_gs),N*len(red_gs)),dtype=float)
    
    
    for i in range(N) :
        for j in range(len(red_gs)):
            tot_gs[i*len(red_gs)+j] = red_gs[j]
            tot_cov[i*len(red_gs)+j,i*len(red_gs)+j] = red_cov[j,j]
    
    return tot_gs, tot_cov

