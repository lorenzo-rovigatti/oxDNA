#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 17:31:11 2023

@author: andrea

####################################
reweight program: 
after a first VMMC run with some educately guessed weight, 
compute weights such that the biased order parameter distribution is approximately uniform

this is done by using:

    n_i/w_i = n_i'/w_i'
    n_i'= n_j' for every i,j
    
    => w_i'=w_i*n_1*w_1'/w_1*n_i = n_1^u*w_1'/n_i^u
    
where n_i is the number of biased sampled states with order parameter = i, and n_i^u the corresponding unbiased frequence

we take n_1^u = max(n_i^u) and set w_1'=1. 

input:
-original weigth file
-last_hist file (VMMC oxDNA simulations)

output:
-new weight file

"""
import sys
import numpy as np

if len(sys.argv) != 3:
    print ('Usage: %s wfile histfile' % (sys.argv[0]))
    exit(1)


wfile = sys.argv[1]
histfile = sys.argv[2]

def reweight() :
        
    weights = np.loadtxt(wfile)
    
    #print(weights)
    histo = np.loadtxt(histfile)
   # print(histo)
    
    Nop = len(weights[:,0])
    print(Nop)
    w = np.zeros(Nop)   #new weights
    h_ub = np.zeros(Nop) #unbiased histo (we have to do this because if one state is never reached, it is not printed in last_hist)
    
    new_wfile = open("wfile_new.txt", 'w')
    
    for i in range(0,len(histo[:,2])) :
        #print(i)
        h_ub[int(histo[i,0])] = histo[i,2]   #in the histo file, the third column is the unbiased frequence (distribution)
        
    i_max = -1
    f_max = -1 #max unbiase frequence
    
    for i in range(0,Nop) : 
        if h_ub[i] > f_max :
            f_max = histo[i,2]
            i_max = i
    
    for i in range(0,Nop) :
        if h_ub[i] > 0 : #if hist[i,2] == 0 (can happen if we do not sample for long enough), w[i] is 0 for now
            w[i] = h_ub[i_max]/h_ub[i] #always >= 1
        
    w_max = -1

    for i in range(0,Nop) :
        if w[i] > w_max :
            w_max = w[i]        
        
    for i in range(0,Nop) :
        
        if w[i] < 0.1 :  #here we set the weigth for i | h_ub[i] == 0 to 2* max weight
            w[i] = w_max*2
        
        line = str(int(weights[i,0])) + " " +str(w[i])
        print(line, file=new_wfile)
        
    return    
        
    
reweight()