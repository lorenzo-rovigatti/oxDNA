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
import functions_multi_melting_mpi as functions_melting
import config_multi as cg
#import Utils

cg.size = cg.comm.Get_size()
cg.rank = cg.comm.Get_rank()

# READ CONFIG FILE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 optimise.py config_file")
    sys.exit()
    
    
if cg.rank == 0 :
    start_time = time.time()
        

config_file = sys.argv[1]

#READ PARAMETERS
if functions_melting.read_config(config_file) == False :
    sys.exit()

#SET UP MPI COMMUNICATORS

cg.seq_id = int(cg.rank / cg.Nreps) 
cg.rep_id = cg.rank % cg.Nreps

cg.comm_seq = cg.comm.Split(cg.seq_id, cg.rank) #create second communicator. Communication between all reps of a given sequence
#set communication between leaders of each seq group

cg.rank_seq = cg.comm_seq.Get_rank()

cg.comm_leaders = cg.comm.Split(cg.rep_id, cg.rank) #create third communicator. Colour 0 set communication between leaders of each sequence (i.e. rep 0)

cg.rank_leaders = cg.comm_leaders.Get_rank()

"""
for i in range(cg.Nseq) :
    cg.leaders.append(i*cg.Nreps)
"""

if cg.rank_seq == 0:
    cg.leaders.append(cg.rank)
    
print("rank: " + str(cg.rank) + ", seq: " + str(cg.seq_id) + "rep: " + str(cg.rep_id) )
if cg.rank in cg.leaders :
    print("rank " + str(cg.rank) + " is a comm leader" )


#compute energy and hbs (order parameter) of sampled trajectory by using oxpy


with oxpy.Context():
    
    l = cg.seq_id 
    i = cg.rep_id
    
    energy_ratio = 300.0/(cg.simTs[l]+273.15) #300K/simT in K 
    print("Energy ratio sampled: "+str(energy_ratio))
        
    file_name = './Seq'+str(l)+'/Rep'+str(i)+'/input1_melting.an'
    functions_melting.update_T_input_file(cg.simTs[l],file_name)
    
    #read input script specifying sequence dependent file
    inp = oxpy.InputFile()
    
    inp.init_from_filename(file_name)
    #create analysys backend
    backend = oxpy.analysis.AnalysisBackend(inp)

    obs = backend.config_info().observables
    
    read = False
    counts = -1

    while 1==1 : 
        try:
            read =  backend.read_next_configuration()
        except:
            counts+=1
            cg.energy_sampled.append(999)
            print("Warning: exception in oxpy energy computation;. Rep "+str(i)+",conf" +str(counts))
            continue
        if read == False :
            break
        counts+=1
        
        if(counts < cg.in_snap) :
            continue
        a = float(obs[0].get_output_string(backend.conf_step).split()[0])   #total energy per nucleotide
        #stk = float(obs[nrep][1].get_output_string(backend[nrep].conf_step).split()[2])   #total stacking energy per nucleotide
        b = int(obs[1].get_output_string(backend.conf_step).split()[0])   #number of hbonds (value of the order parameter)
        
        cg.energy_sampled.append((cg.Njuns[l]+1)*20*a*energy_ratio)
        #cg.energy_stk_sampled[l].append((cg.Njuns[l]+1)*20*stk)
        cg.hbs_sampled.append(b)
                
                            
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
    
    
if cg.rank == 0:
    
    stop = [0]
    
    #sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(par0),method='nelder-mead',options={'maxiter':cg.miter, 'eps':0.1})
    if cg.algo == "L-BFGS-B" :
        sol = optimize.minimize(functions_melting.Cost_function_mT,par,args=(stop,par0),method='L-BFGS-B', callback=functions_melting.callbackF, bounds=bnd ,options={'maxfun':cg.neva, 'eps':cg.LBFGSB_eps, 'iprint':cg.LBFGSB_iprint})
        
    elif cg.algo == "nelder-mead"  :
        sol = optimize.minimize(functions_melting.Cost_function_mT,par,args=(stop,par0),method='nelder-mead', callback=functions_melting.callbackF, bounds=bnd ,options={'maxfev':cg.neva})
    
    else :
        print("UNKNOWN ALGORITHM!")
        sys.exit()
    
    #OPTI DONE
    print("OPTI DONE")
    print("TERMINATING")
    stop = [1] #this stops the while on the other processors
    functions_melting.Cost_function_mT(par,stop,par0)   #extra call to stop other processors
        
        
else :
    #ite = 0
    stop = [0]
    while stop[0]==0:   #keep computing local term of cost function untill stop[0] != 0
        #ite += 1
        C = functions_melting.Cost_function_mT(par,stop,par0)
        #print("ite",cg.rank,ite)


#compute and gather final melting temperatures
#update final Ts file (needed for next iteration of optimisation)

solution = None

if cg.rank == 0 :
    solution = sol.x

solution = cg.comm.bcast(solution, root=0)
print(cg.rank,solution)

cg.comm.Barrier()

stop = [0]

C = functions_melting.Cost_function_mT(solution,stop,par0)

#print("final C: ", C)

mTs = cg.comm.gather(cg.current_mT,root=0)

if cg.rank == 0 :
    
    print("mTs")
    print(mTs)

    print("gathered")

    print(sol)
    functions_melting.print_final_melting_temperatures(mTs)

    #print(sol.x)
    
    if cg.ave :
        functions.update_rew_seq_dep_file_ave(sol.x)
    else :
        functions.update_rew_seq_dep_file(sol.x)
        
    timefile = open("runtime_mpi.txt",'w')
    print("Run time [s]: " + str(time.time() - start_time),file=timefile)
    timefile.close()
    
    print("STEP DONE!")

