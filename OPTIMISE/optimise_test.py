#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:12:36 2023

@author: yqb22156
"""

import sys
from scipy import optimize
import oxpy
import numpy as np

import functions_test as functions
from oxdna_to_internal_wflip import read_oxdna_trajectory_standard_order
import config as cg


# READ CONFIG FILE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 optimise.py config_file")
    sys.exit()
    

config_file = sys.argv[1]

#READ PARAMETERS
if functions.read_config(config_file) == False :
    sys.exit()


#Compute and store internal coordinates from sampled trajectories. Compute gs and cov
for i in range(cg.Nreps) :

    iname = './Rep'+str(i)+'/trajectory.dat'
    tname = './Rep'+str(i)+'/generated.top'
    
    ifile = open(iname,'r')
    tfile = open(tname,'r')
    
    traj = read_oxdna_trajectory_standard_order(ifile, tfile)
    
    ifile.close()
    tfile.close()
    
    # True = overwrite, False = append
    # append is for averaging over multiple trajectories
    
    if i == 0 :
        functions.store_internal_coord(traj,cg.ids,cg.in_j,cg.fin_j,cg.in_snap,True)
    else  :
        functions.store_internal_coord(traj,cg.ids,cg.in_j,cg.fin_j,cg.in_snap,False)
 
#average coordinates and energy
functions.ave_and_cov_stored()

print("mu sampled:")
print(cg.mu_sampled)

print("cov sampled:")
print(cg.cov_sampled)


inp = []
backend = []
obs = []

#compute energy of trajectory by using oxpy
for i in range(cg.Nreps) :
    with oxpy.Context():
        
        #read input script specifying sequence dependent file
        inp.append(oxpy.InputFile())
        
        inp[i].init_from_filename("./Rep"+str(i)+"/input1.an")
        #create analysys backend
        backend.append(oxpy.analysis.AnalysisBackend(inp[i]))
    
        obs.append(backend[i].config_info().observables)
        
        read = False
        counts = -1
        """
        while backend[i].read_next_configuration() :
            
            counts+=1
            
            if(counts < in_snap) :
                continue

            a = float(obs[i][0].get_output_string(backend[i].conf_step).split()[0])
            print(a)
            energy_sampled.append((Njuns+1)*10*a)
        """
        while 1==1 : 
            try:
                read =  backend[i].read_next_configuration()
            except:
                counts+=1
                functions.energy_sampled.append(999)
                print("Warning: exception in oxpy energy computation;. Rep "+str(i)+",conf" +str(counts))
                continue
            if read == False :
                break
            counts+=1
            
            if(counts < cg.in_snap) :
                continue
            a = float(obs[i][0].get_output_string(backend[i].conf_step).split()[0])
            cg.energy_sampled.append((cg.Njuns+1)*20*a)
            
       

#reduce covariance for stiff part: remove from the covariance all coordinates we want to exclude from cov optimisation


mu = cg.mu_sampled
cov = cg.cov_sampled

ave_target_mu = np.zeros(len(cg.ids),dtype=float)
ave_target_cov = np.zeros((len(cg.ids),len(cg.ids)),dtype=float)
ave_cov = np.zeros((len(cg.ids),len(cg.ids)),dtype=float)


for i in range(len(cov[0])) :
    ave_target_mu[i%len(cg.ids)] += cg.target_mu[i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
    if cg.diag == True:
        ave_cov[i%len(cg.ids),i%len(cg.ids)] += cov[i,i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
        ave_target_cov[i%len(cg.ids),i%len(cg.ids)] += cg.target_cov[i,i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))

red_n = (cg.fin_j-cg.in_j+1)*(len(cg.ids)-len(cg.ids_cov))

reduced_cov = np.zeros((cg.dimension-red_n,cg.dimension-red_n),dtype=float)
reduced_target_cov = np.zeros((cg.dimension-red_n,cg.dimension-red_n),dtype=float)

ave_reduced_cov = np.zeros((len(cg.ids_cov),len(cg.ids_cov)),dtype=float)
ave_reduced_target_cov = np.zeros((len(cg.ids_cov),len(cg.ids_cov)),dtype=float)

nrow=-1

for i in range(cg.dimension):
    if cg.ids[i%len(cg.ids)]  in cg.ids_cov: #remove rise from covariance

        nrow+=1
        ncol=-1
        for j in range(cg.dimension):
            if cg.ids[j%len(cg.ids)]  in cg.ids_cov: #remove rise from covariance
                ncol+=1
                reduced_cov[nrow,ncol] = cov[i,j]
                reduced_target_cov[nrow,ncol] = cg.target_cov[i,j]
nrow=-1          
for i in range(len(cg.ids)) :
    if cg.ids[i] in cg.ids_cov :
        nrow+=1
        ncol=-1
        for j in range(len(cg.ids)):
            if cg.ids[j]  in cg.ids_cov: #remove rise from covariance
                ncol+=1
                ave_reduced_cov[nrow,ncol] = ave_cov[i,j]
                ave_reduced_target_cov[nrow,ncol] = ave_target_cov[i,j]
                
                
                
#PERSISTENCE LENGTHS            
opti_lp = True
    
for z in range(len(cg.ids_inter_rot)) :
    if cg.ids_inter_rot[z] in cg.ids_cov :
        continue
    else :
        print("Warning: Not all inter rotations are used for covariance optimisation. Cannot tune persistence length.")
        opti_lp = False
    
    
if opti_lp :
    
    M = functions.Stiff(reduced_cov,cg.Njuns-10)
    lb,lt = functions.lps(M)
    
    lt_half = lt/2
    
    print("Long range (m="+str(cg.Njuns-10)+"):")
    print(M)
    print("lb: "+str(lb))     
    print("lt/2: "+str(lt_half))


"""
#read initial values of the optimisation parameters (from initial seq dep file)

ifile = open("oxDNA_sequence_dependent_parameters_in.txt",'r')

par0 = []
order = []

for line in ifile.readlines() :
    vals = line.split()
    
    #average    
    if cg.ave :
        if len(vals) == 0 :
                continue
        for i in range(len(cg.par_codename)) :
            if vals[0] == cg.par_codename[i] and cg.par_codename[i] == 'FENE_DELTA' :
                par0.append(float(vals[2]))
                order.append(i)
            elif vals[0] == cg.par_codename[i]+'_A_A' :
                par0.append(float(vals[2]))    
                order.append(i)
    #SD           
    else :
        if len(vals) == 0 :
                continue
        for i in range(len(cg.par_codename)) :
            if vals[0] == cg.par_codename[i] :
                par0.append(float(vals[2]))
                order.append(i)

ifile.close()

#sort par0 so that the values match the parameters in par_codenames 
par0 = [par0[i] for i in order]

print("Initial values of the optimisation parameters: ")
print(par0)

par = par0

print("RUNNING MINIMISATION")

sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(par0),method='nelder-mead',options={'maxiter':cg.miter})

print(sol)

#print(sol.x)

if cg.ave :
    functions.update_rew_seq_dep_file_ave(sol.x)
else :
    functions.update_rew_seq_dep_file(sol.x)

print("STEP DONE!")

"""