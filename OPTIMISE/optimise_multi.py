#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:12:36 2023

@author: yqb22156
"""

import sys
from scipy import optimize
import oxpy
import copy

import functions_multi as functions
from oxdna_to_internal_wflip import read_oxdna_trajectory_standard_order
import config_multi as cg


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
for l in range(cg.Nseq) : 
    for i in range(cg.Nreps) :
    
        iname = './Seq'+str(l)+'/Rep'+str(i)+'/trajectory.dat'
        tname = './Seq'+str(l)+'/Rep'+str(i)+'/generated.top'
        
        ifile = open(iname,'r')
        tfile = open(tname,'r')
        
        traj = read_oxdna_trajectory_standard_order(ifile, tfile)
        
        ifile.close()
        tfile.close()
        
        # True = overwrite, False = append
        # append is for averaging over multiple trajectories
        
        print("READING COORDINATES SEQ "+str(l)+" REP" +str(i))
        
        if i == 0 :
            functions.store_internal_coord(traj,l,cg.ids,cg.in_j[l],cg.fin_j[l],cg.in_snap,True)
        else  :
            functions.store_internal_coord(traj,l,cg.ids,cg.in_j[l],cg.fin_j[l],cg.in_snap,False)
 
#average coordinates and energy
functions.ave_and_cov_stored()
for l in range(cg.Nseq) :
    
    print("SEQUENCE "+str(l))
    
    print("mu sampled:")
    print(cg.mu_sampled[l])
    
    print("cov sampled:")
    #functions.print_matrix(cg.cov_sampled)
    print(cg.cov_sampled[l])


inp = []
backend = []
obs = []

#compute energy of trajectory by using oxpy
for l in range(cg.Nseq) :
    for i in range(cg.Nreps) :
        with oxpy.Context():
            
            nrep = l*cg.Nreps +i
            
            #read input script specifying sequence dependent file
            inp.append(oxpy.InputFile())
            
            inp[nrep].init_from_filename('./Seq'+str(l)+'/Rep'+str(i)+'/input1.an')
            #create analysys backend
            backend.append(oxpy.analysis.AnalysisBackend(inp[nrep]))
        
            obs.append(backend[nrep].config_info().observables)
            
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
                    read =  backend[nrep].read_next_configuration()
                except:
                    counts+=1
                    cg.energy_sampled[l].append(999)
                    print("Warning: exception in oxpy energy computation;. Rep "+str(i)+",conf" +str(counts))
                    continue
                if read == False :
                    break
                counts+=1
                
                if(counts < cg.in_snap) :
                    continue
                a = float(obs[nrep][0].get_output_string(backend[nrep].conf_step).split()[0])
                cg.energy_sampled[l].append((cg.Njuns[l]+1)*20*a)
            
       
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

par0_c = copy.deepcopy(par0)

par0 = [par0_c[i] for i in order]

print("Initial values of the optimisation parameters: ")
print(par0)

par = par0

print("RUNNING MINIMISATION")

sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(par0),method='nelder-mead',options={'maxiter':cg.miter, 'eps':0.1})

print(sol)

#print(sol.x)

if cg.ave :
    functions.update_rew_seq_dep_file_ave(sol.x)
else :
    functions.update_rew_seq_dep_file(sol.x)

print("STEP DONE!")

