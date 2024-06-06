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
from oxdna_to_internal_wflip import read_oxdna_trajectory_standard_order
import config_multi as cg
import Utils_mpi
from mpi4py import MPI

cg.size = cg.comm.Get_size()
cg.rank = cg.comm.Get_rank()

# READ CONFIG FILE
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 optimise.py config_file")
    sys.exit(1)

start_time = 0.

if cg.rank == 0 :
    start_time = time.time()

config_file = sys.argv[1]

#READ PARAMETERS
if functions.read_config(config_file) == False :
    sys.exit(1)
    
    
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


#Compute and store internal coordinates from sampled trajectories. Compute gs and cov
l = cg.seq_id 
i = cg.rep_id
    
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


functions.store_internal_coord(traj,cg.ids,cg.in_j[l],cg.fin_j[l],cg.in_snap,True)


#check that the simulation went well:
    
lsamp = len(cg.internal_coords)

if lsamp < cg.in_snap + 100:
    print("FATAL:")
    print("Not enough sampled configurations.")
    print("At least 100+IN_SNAP sampled configurations are needed." )
    print("Did something go wrong with the simulation?")
    sys.exit(1)    


#compute energy of trajectory by using oxpy

#accepted minumn value of hb bonds.
#if the number of hb of a configuration is less then min_hb, than it is discarded
#since configurations melt from the ends, this we can make it so the protion of DNA used for optimisation is always a double helix
"""
min_rj = cg.inj
if cg.jfe < min_rj :
    min_rj = cg.jfe
    
min_hb = cg.Njuns[cg.seq_id ]+1
if min_rj  > 0 :
    min_hb = min_hb - min_rj + 1
"""
    
min_hb = cg.Njuns[cg.seq_id ]+1 -  cg.inj - cg.jfe
    
    
discarded = 0

with oxpy.Context():
    
    l = cg.seq_id 
    i = cg.rep_id
    
    
    #read input script specifying sequence dependent file
    inp = oxpy.InputFile()
    
    inp.init_from_filename('./Seq'+str(l)+'/Rep'+str(i)+'/input1.an')
    #create analysys backend
    backend = oxpy.analysis.AnalysisBackend(inp)

    obs = backend.config_info().observables
    
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
        a = float(obs[0].get_output_string(backend.conf_step).split()[0])
        b = int(obs[1].get_output_string(backend.conf_step).split()[0])   #number of hbonds (value of the order parameter)
        
        if b < min_hb :
            cg.energy_sampled.append(999)
            print("rank " + str(cg.rank) + ", seq"+ str(cg.seq_id)+ " rep"+str(cg.rep_id)+", discarded one conf. Hb = "+str(b))
            discarded += 1
        else:        
            cg.energy_sampled.append((cg.Njuns[l]+1)*20*a)
            
print("rank " + str(cg.rank) + ", seq"+ str(cg.seq_id)+ " rep"+str(cg.rep_id)+ " total discarded confs: "+str(discarded))


good_confs = lsamp - cg.in_snap - discarded

good_confs = cg.comm_seq.reduce(good_confs,op=MPI.SUM, root=0)

if cg.rank_seq == 0: #same as cg.rank in cg.leaders
    if lsamp - cg.in_snap - discarded <  50:
        print("FATAL:")
        print("Discarded too many configurations.")
        print("Not enough sampled configurations left.")
        print("At least 100 configurations are needed for acceptable statistics." )
        print("Did something go wrong with the simulation?")
        sys.exit(1)   
            
#average coordinates and energy
functions.ave_and_cov_stored()
            
#compute mu and cov for sequence (summing over replicas)      
cg.mu_sampled = cg.comm_seq.reduce(cg.mu_sampled,op=MPI.SUM, root=0)
cg.cov0_sampled = cg.comm_seq.reduce(cg.cov0_sampled,op=MPI.SUM, root=0)

if cg.rank_seq == 0: #same as cg.rank in cg.leaders

    cg.mu_sampled /= cg.Nreps
    cg.cov0_sampled /= cg.Nreps  
    
    for p in range(len(cg.cov0_sampled)) :
        for q in range(len(cg.cov0_sampled[p])) :
            cg.cov_sampled[p][q] =  cg.cov0_sampled[p][q] - cg.mu_sampled[p]*cg.mu_sampled[q]    

    Utils_mpi.plot_gs_sampled()
    Utils_mpi.plot_cov_sampled_diag()
    
    print("SEQUENCE "+str(cg.seq_id))

    print("mu sampled:")
    print(cg.mu_sampled)

    print("cov sampled:")
    #functions.print_matrix(cg.cov_sampled)
    print(cg.cov_sampled)
    
    cg.mu_curr = cg.mu_sampled
    
    print("Computing Deltas "+str(cg.seq_id))
    
    functions.compute_deltas(True)
    
    print("Communicating Deltas "+str(cg.seq_id))
    
    cg.Deltas = cg.comm_leaders.reduce(cg.Deltas,op=MPI.SUM, root=0)
    cg.bsteps_counts = cg.comm_leaders.reduce(cg.bsteps_counts,op=MPI.SUM, root=0)
    
    
    if cg.rank == 0:
        for i in range(len(cg.ids)):
            for j in range(4):
                for k in range(4):
                    if cg.bsteps_counts[j][k] > 0:
                        cg.Deltas[i][j][k]/=cg.bsteps_counts[j][k]
            print(cg.Deltas[i])
        print(cg.bsteps_counts)
                    
    
       
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
                    
                    up_bond.append(float(vals[2])*1.02)
                    low_bond.append(float(vals[2])*0.98)
                
                elif vals1[0] == "FENE" and vals1[1] == "DELTA" :    #stricter for FENE_DELTA (+-5%), to avoid problems with the FENE potential
                    
                    up_bond.append(float(vals[2])*1.02)
                    low_bond.append(float(vals[2])*0.98)
                    
                elif vals1[0] == "STCK" and vals1[1] == "R0" :
                    
                    up_bond.append(float(vals[2])*1.1)
                    low_bond.append(float(vals[2])*0.9) 
                    
                elif vals1[0] == "STCK" and vals1[2] == "T0" and abs(float(vals[2])) < 0.01:   #stricter for FENE_R0 (+-3%), to avoid problems with the FENE potential
                    
                    up_bond.append(0.2)
                    low_bond.append(-0.2)
                    
                elif vals1[0] == "CRST" and vals1[1] == "THETA4" and vals1[2] == "T0" :   #stricter for FENE_R0 (+-3%), to avoid problems with the FENE potential
                    
                    up_bond.append(3.14159+0.1745) #pi +- 20deg to avoid weird things
                    low_bond.append(3.14159-0.1745)
                    
                else :                    
                    up_bond.append(float(vals[2])*6.0)
                    low_bond.append(float(vals[2])*0.1)
                    
                
                order.append(i)

ifile.close()

par0_c = copy.deepcopy(par0)


print(order)
print(par0_c)

for i in range(len(order)) :
    par0[order[i]] = par0_c[i]



#print(par0)
print(["{0:0.3f}".format(i) for i in par0])

#sort par0 so that the values match the parameters in par_codenames 


up_bond_c = copy.deepcopy(up_bond)
low_bond_c = copy.deepcopy(low_bond)

up_bond = [up_bond_c[i] for i in order]
low_bond = [low_bond_c[i] for i in order]



for i in range(len(order)) :
    up_bond[order[i]] = up_bond_c[i]
    low_bond[order[i]] = low_bond_c[i]




print("Initial values of the optimisation parameters: ")

par = copy.deepcopy(par0)


#print(par)

bnd = optimize.Bounds(low_bond,up_bond)


cg.curr_feva = 0

print("RUNNING MINIMISATION")
    
for n in range(0,4) :
    
    S = functions.Relative_entropy_wRew(par,[0],par0)
    
    if n > 0:
        if cg.rank_seq == 0: #same as cg.rank in cg.leaders
            functions.compute_deltas(False)
            
            cg.Deltas = cg.comm_leaders.reduce(cg.Deltas,op=MPI.SUM, root=0)
            
            if cg.rank == 0:
                for i in range(len(cg.ids)):
                    for j in range(4):
                        for k in range(4):
                            if cg.bsteps_counts[j][k] > 0:
                                cg.Deltas[i][j][k]/=cg.bsteps_counts[j][k]    
    
    if cg.rank == 0:
                
        stop = [0]
        
        
        #sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(par0),method='nelder-mead',options={'maxiter':cg.miter, 'eps':0.1})
        if cg.algo == "L-BFGS-B" :
            sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(stop,par0),method='L-BFGS-B', callback=functions.callbackF, bounds=bnd ,options={'maxfun':cg.neva, 'eps':cg.LBFGSB_eps, 'iprint':cg.LBFGSB_iprint})
            
        elif cg.algo == "nelder-mead"  :
            
            in_simplex = functions.build_initial_simplex_for_nm(par,up_bond,low_bond)
            
            
            #sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(stop,par0),method='nelder-mead', callback=functions.callbackF, bounds=bnd ,options={'maxfev': cg.neva, 'adaptive': True})
            #sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(stop,par0),method='nelder-mead', callback=functions.callbackF, bounds=bnd ,options={'maxfev':cg.neva, 'initial_simplex' : in_simplex})
            #sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(stop,par0),method='nelder-mead', callback=functions.callbackF, bounds=bnd ,options={'maxiter': cg.miter, 'adaptive' : True, 'initial_simplex' : in_simplex})
            sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(stop,par0),method='nelder-mead', callback=functions.callbackF, bounds=bnd ,options={'maxiter': cg.miter, 'initial_simplex' : in_simplex})
            #sol = optimize.minimize(functions.Relative_entropy_wRew,par,args=(stop,par0),method='nelder-mead', callback=functions.callbackF, bounds=bnd ,options={'maxfev':2, 'initial_simplex' : in_simplex})
        
        else :
            print("UNKNOWN ALGORITHM!")
            sys.exit(1)
        
        #OPTI DONE
        print("OPTI DONE")
        print("TERMINATING")
        stop = [1] #this stops the while on the other processors
        functions.Relative_entropy_wRew(par,stop,par0)   #extra call to stop other processors
        par = sol.x
            
            
    else :
        #ite = 0
        stop = [0]
        while stop[0]==0:   #keep computing local term of cost function untill stop[0] != 0
            #ite += 1
            C = functions.Relative_entropy_wRew(par,stop,par0)
            #print("ite",cg.rank,ite)



if cg.rank == 0 :
    
    print("gathered")

    print(sol)

    #print(sol.x)
    
    if cg.ave :
        functions.update_rew_seq_dep_file_ave(sol.x)
    else :
        functions.update_rew_seq_dep_file(sol.x)
        
        
    timefile = open("runtime_mpi.txt",'w')
    print("Run time [s]: " + str(time.time() - start_time),file=timefile)
    timefile.close()
    
    print("STEP DONE!")