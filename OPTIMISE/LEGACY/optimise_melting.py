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
import time

import functions_multi as functions
import functions_multi_melting as functions_melting
import config_multi as cg
import Utils


# READ CONFIG FILE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 optimise.py config_file")
    sys.exit()
    

start_time = time.time()
    

config_file = sys.argv[1]

#READ PARAMETERS
if functions_melting.read_config(config_file) == False :
    sys.exit()


inp = []
backend = []
obs = []

#compute energy and hbs (order parameter) of sampled trajectory by using oxpy
for l in range(cg.Nseq) :
    
    energy_ratio = 300.0/(cg.simTs[l]+273.15) #300K/simT in K 
    print("Energy ratio sampled: "+str(energy_ratio))
    
    for i in range(cg.Nreps) :
        with oxpy.Context():
            
            nrep = l*cg.Nreps +i
            
            file_name = './Seq'+str(l)+'/Rep'+str(i)+'/input1_melting.an'
            functions_melting.update_T_input_file(cg.simTs[l],file_name)
            
            #read input script specifying sequence dependent file
            inp.append(oxpy.InputFile())
            
            inp[nrep].init_from_filename(file_name)
            
            print(l,i)
            
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
                a = float(obs[nrep][0].get_output_string(backend[nrep].conf_step).split()[0])   #total energy per nucleotide
                #stk = float(obs[nrep][1].get_output_string(backend[nrep].conf_step).split()[2])   #total stacking energy per nucleotide
                b = int(obs[nrep][1].get_output_string(backend[nrep].conf_step).split()[0])   #number of hbonds (value of the order parameter)
                
                cg.energy_sampled[l].append((cg.Njuns[l]+1)*20*a*energy_ratio)
                #cg.energy_stk_sampled[l].append((cg.Njuns[l]+1)*20*stk)
                cg.hbs_sampled[l].append(b)
                
                            
#Utils.plot_gs_sampled()
#Utils.plot_mt_sampled() #XXXtodo plot comparison with Santa Lucia
            
       
#read initial values of the optimisation parameters (from initial seq dep file)
ifile = open("oxDNA_sequence_dependent_parameters_in.txt",'r')

par0 = []
order = []

up_bond = []    #upper parameter bond
low_bond = []   #lower parameter bond

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
                
                print(i)
                print(cg.par_codename[i])
                
                vals1 = cg.par_codename[i].split('_')
                
                if vals1[0] == "FENE" and vals1[1] == "R0" :    #stricter for FENE_R0 (+-3%), to avoid problems with the FENE potential
                    
                    up_bond.append(float(vals[2])*1.03)
                    low_bond.append(float(vals[2])*0.97)
                
                elif vals1[0] == "FENE" and vals1[1] == "DELTA" :    #stricter for FENE_DELTA (+-5%), to avoid problems with the FENE potential
                    
                    up_bond.append(float(vals[2])*1.05)
                    low_bond.append(float(vals[2])*0.95)
                    
                elif vals1[0] == "STCK" and vals1[1] == "R0" :
                    
                    up_bond.append(float(vals[2])*1.05)
                    low_bond.append(float(vals[2])*0.95)
                    
                elif (vals1[0] == "STCK" or vals1[0] == "HYDR") and len(vals) == 3 :
                    
                    up_bond.append(float(vals[2])*3)
                    low_bond.append(0.)
                    
                    
                else : 
                    
                    up_bond.append(float(vals[2])*1.5)
                    low_bond.append(float(vals[2])*0.5)
                    
                
                order.append(i)

ifile.close()

#sort par0 so that the values match the parameters in par_codenames 

par0_c = copy.deepcopy(par0)
up_bond_c = copy.deepcopy(up_bond)
low_bond_c = copy.deepcopy(low_bond)

up_bond = [up_bond_c[i] for i in order]
low_bond = [low_bond_c[i] for i in order]


print(order)
print(par0_c)

for i in range(len(order)) :
    par0[order[i]] = par0_c[i]
    up_bond[order[i]] = up_bond_c[i]
    low_bond[order[i]] = low_bond_c[i]



print(par0)

print("Initial values of the optimisation parameters: ")
print(par0)

par = par0

bnd = optimize.Bounds(low_bond,up_bond)


print("RUNNING MINIMISATION")

#sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(par0),method='nelder-mead',options={'maxiter':cg.miter, 'eps':0.1})
if cg.algo == "L-BFGS-B" :
    sol = optimize.minimize(functions_melting.Cost_function_mT,par,args=(par0),method='L-BFGS-B', callback=functions_melting.callbackF, bounds=bnd ,options={'maxfun':cg.neva, 'eps':cg.LBFGSB_eps, 'iprint':cg.LBFGSB_iprint})
    
elif cg.algo == "nelder-mead"  :
    sol = optimize.minimize(functions_melting.Cost_function_mT,par,args=(par0),method='nelder-mead', callback=functions_melting.callbackF, bounds=bnd ,options={'maxfev':cg.neva})

else :
    print("UNKNOWN ALGORITHM!")
    sys.exit()


print(sol)

mTs, mTs_w = functions_melting.reweight_melting_temperature(sol.x)
"""
for l in range(len(mTs)) :
    for i in range(len(mTs[l])) :
        mTs[l][i] = mTs[l][i] - 273.15 #from K to C
"""

functions_melting.print_final_melting_temperatures(mTs)

#print(sol.x)

if cg.ave :
    functions.update_rew_seq_dep_file_ave(sol.x)
else :
    functions.update_rew_seq_dep_file(sol.x)
    
    
timefile = open("runtime.txt",'w')
print("Run time [s]: " + str(time.time() - start_time),file=timefile)
timefile.close()

print("STEP DONE!")

