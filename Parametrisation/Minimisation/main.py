#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:12:36 2023

@author: yqb22156
"""

import numpy as np
import math
import copy
from scipy import optimize
import oxpy
import matplotlib.pyplot as plt


from oxdna_to_internal_wflip import read_oxdna_trajectory_standard_order
import continuity_constraints

bases = ['A','C','G','T']

internal_coords = [] #stores arrays of internal coordinates for sampled configurations
energy_sampled = [] #stores the enegy of the sampled configurations

dimension = 28 #number of coordinates to optimise
Njuns = 19 #total number of junctions
mu_sampled = np.zeros(dimension, dtype = float)
cov_sampled = np.zeros((dimension,dimension), dtype = float)

# IDS:
# 0,1,2 = intra tran
# 3,4,5 = intra rot
# 6,7,8 = inter tran
# 9,10,11 = inter rot
ids = [4,11]
in_j = 3 #ignore ends
fin_j = Njuns -3 #ignore ends
in_snap = 100 #ignore first in_snap snapshots (equilibration)

Nreps = 4

#associate each entry of par to a specific oxdna parameter
codename_file = open('codenames.txt','r')
par_codename = []
continuity_par_codenames = []
continuity_par_values = []

for line in codename_file.readlines() :
    vals = line.split()
    if vals[0][0] == '#':
        continue
    #print(vals)
    if vals[1] == 'OPTIMISE' :
        par_codename.append(vals[0])
    else :
        continuity_par_codenames.append(vals[0])
        continuity_par_values.append(float(vals[2]))
    
codename_file.close()

print("Optimising:")
print(par_codename)
print("Auxiliary parameters (continuity):")
print(continuity_par_codenames)
print(continuity_par_values)

par_dimension = len(par_codename) #number of parameters to use

used = [] #track used parameters when imposing continuity (avoid printing duplication)
for i in range(par_dimension) :
    used.append(False)

#read file containing ausiliar parameter when changing parameters requiring continuity

#Mean and covariance of the target distribution
target_cov = np.zeros([dimension,dimension], dtype=float)
target_mu = np.zeros(dimension, dtype=float)

#read target distribution
"""
ifile = open('target_mean_and_covariance_twist.txt','r')

ln = 0
for line in ifile.readlines() :   
    vals = line.split()
    
    if len(vals) == 0 :
        continue
    
    if vals[0][0] == '#' :
        continue
    
    if len(vals) != dimension :
        print("Warning: mismatch between dimension and target distribution dimension!")
    if ln == 0 :            
        for i in range(len(vals)) :
            target_mu[i] = float(vals[i])
    else :
        for i in range(len(vals)) :
            target_cov[ln-1][i] = float(vals[i])
    ln+=1
    
ifile.close()

target_inv_cov = np.linalg.inv(target_cov)
det_target_inv_cov = np.linalg.det(target_inv_cov)
"""



#given a junction trajectory (read_oxdna_trajectory_standard_order), store specific internal coordinates in
#global variable internal_coords
def store_internal_coord(traj,ids,in_j,fin_j,in_snap,overwrite=True) :
    
    global internal_coords
     
    coords = []
    
    Nsnaps = len(traj)
    Njuns = len(traj[0])
    
    for i in range(in_snap,Nsnaps) :

        coord = []
        for j in range(Njuns) :
            #print("snap:"+str(i)+" junn:"+str(j))
            if j < in_j or j > fin_j :
                continue
            for z in range(len(ids)) :
                if ids[z] == 0 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[0])
                elif ids[z] == 1 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[1])
                elif ids[z] == 2 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[2])
                elif ids[z] == 3 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[0])
                elif ids[z] == 4 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[1])
                elif ids[z] == 5 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[2])
                    
                elif ids[z] == 6 :
                    coord.append(traj[i][j].inter_coord.tran[0])
                elif ids[z] == 7 :
                    coord.append(traj[i][j].inter_coord.tran[1])
                elif ids[z] == 8 :
                    coord.append(traj[i][j].inter_coord.tran[2])
                elif ids[z] == 9 :
                    coord.append(traj[i][j].inter_coord.rot[0])
                elif ids[z] == 10 :
                    coord.append(traj[i][j].inter_coord.rot[1])
                elif ids[z] == 11 :
                    coord.append(traj[i][j].inter_coord.rot[2])
                    
        if overwrite == False :
            internal_coords.append(coord)
        else :
            coords.append(coord)
           
    if overwrite == True :
        internal_coords = coords
                    
    return

def ave_and_cov_stored() :
    
    global mu_sampled
    global cov_sampled
    
    Nsnaps = len(internal_coords)
    Ncoords = len(internal_coords[0])
    
    for i in range(Ncoords) :
        mu_sampled[i] = 0.
        for j in range(Ncoords) :
            cov_sampled[i][j] = 0.
    
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            mu_sampled[j] += internal_coords[i][j]/Nsnaps
    
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            for z in range(j+1,Ncoords) :
                cov_sampled[j][z] += (internal_coords[i][j] - mu_sampled[j])*(internal_coords[i][z] - mu_sampled[z])/Nsnaps
    
    for j in range(Ncoords) :
        for z in range(j+1,Ncoords) :
            cov_sampled[z][j] = cov_sampled[j][z]    
    
    return

#update sequence dependent file for computing energy with new value of parameters
def update_rew_seq_dep_file(par) :
    
    ofile = open('oxDNA2_sequence_dependent_parameters_tmp.txt','w')
    ifile = open('oxDNA2_sequence_dependent_parameters.txt','r')    #file with fixed sequence dependent parameters
    
    if len(par) != len(par_codename) :
        print("Something is not right with the parameters! Check codename file.")
    
    else :
        for line in ifile.readlines() :
            print(line.strip('\n'), file=ofile)
            
        print('\n', file=ofile)
        
        for i in range(len(par_codename)) :
            print(par_codename[i]+" = "+str(par[i]),file=ofile)
            
    ofile.close()
    ifile.close()
    
    return

#update sequence dependent file for computing energy with new value of parameters
#version with average parameters

def impose_continuity(par_cname,p_id,pars) :
    vals = par_cname.split('_')
    auxiliars = []
    aux_values = []
    output = []
    
    global used
    
    f1 = False
    f4 = False
    
    #f1   
    r0 = 0.
    rc = 0.
    a = 0.

    if vals[1] == 'R0' :
        f1 = True
        r0 = pars[p_id]
        auxiliars.append('A')
        auxiliars.append('RC')
    elif vals[1] == 'RC' :
        f1 = True
        rc = pars[p_id]
        auxiliars.append('R0')
        auxiliars.append('A')
    elif vals[1] == 'A' :
        f1 = True
        a = pars[p_id]
        auxiliars.append('R0')
        auxiliars.append('RC')
    if f1 :
        for i in range(len(auxiliars)) : 
            found = False
            aux_pname = vals[0]+'_'+auxiliars[i]
            for j in range(2,len(vals)) :
                aux_pname = aux_pname + '_' + vals[j]
            for j in range(len(continuity_par_codenames)) :
                if aux_pname == continuity_par_codenames[j] :
                    aux_values.append(continuity_par_values[j])
                    found = True
            if found == False :
                for j in range(len(par_codename)) :
                    if aux_pname == par_codename[j] :
                        aux_values.append(pars[j])
                        used[j] = True
                        found = True
            if found == False :
                print("WARNNG: Can't impose contnuity, parameter(s) missing!")
            
      
        for i in range(len(auxiliars)) :
            if auxiliars[i] == 'R0' :
                r0 = aux_values[i]
            elif auxiliars[i] == 'RC' :
                rc = aux_values[i]
            elif auxiliars[i] == 'A' :
                a = aux_values[i]
            
        rl,rh,bl,bh,rcl,rch,continuous = continuity_constraints.continuity_f1(r0,a,rc)
        
        output.append('f1')
        output.append(rl)
        output.append(rh)
        output.append(bl)
        output.append(bh)
        output.append(rcl)
        output.append(rch)
        
        return output

    #f4
    a = 0.

    if vals[2] == 'A' and (vals[1] != 'HYDR' or vals[1] != 'STCK'):
        f4 = True
        a = pars[p_id]
    """
        auxiliars.append('T0')
    elif vals[2] == 'T0' and (vals[1] != 'HYDR' or vals[1] != 'STCK'):
        f4 = True
        th0 = pars[p_id]
        auxiliars.append('A')
    """
    if f4 :
        """
        for i in len(auxiliars) : 
            found = False
            aux_pname = vals[0]+'_'+vals[1]+'_'+auxiliars[i]
            for j in range(2,len(vals)) :
                aux_pname = aux_pname + '_' + vals[j]
            for j in range(len(continuity_par_codenames)) :
                if aux_pname == continuity_par_codenames[j] :
                    aux_values[i].insert(continuity_par_values[j])
                    found = True
            if found == False :
                for j in range(len(par_codename)) :
                    if aux_pname == par_codename[j] :
                        aux_values[i].insert(pars[j])
                        used[j] = True
                        found = True
            if found == False :
                print("WARNNG: Can't impose contnuity, parameter(s) missing!")
            
      
        for i in range(len(auxiliars)) :
            if auxiliars[i] == 'T0' :
                th0 = aux_values[i]
            elif auxiliars[i] == 'A' :
                a = aux_values[i]
        """
        dts,dtc,b= continuity_constraints.continuity_f4(a)
        
        output.append('f4')
        output.append(dts)
        output.append(dtc)
        output.append(b)
        
        return output
    
    if f1 == False and f4 == False:
        output.append('No')
        return output
                
        
    

def update_rew_seq_dep_file_ave(par) :
    
    ofile = open('oxDNA2_sequence_dependent_parameters_tmp.txt','w')
    ifile = open('oxDNA2_sequence_dependent_parameters.txt','r')    #file with fixed sequence dependent parameters
    
    #continuity for f1
    
    if len(par) != len(par_codename) :
        print("Something is not right with the parameters! Check codename file.")
    
    else :
        for line in ifile.readlines() :
            print(line.strip('\n'), file=ofile)
            
        print('\n', file=ofile)
        
        for i in range(len(par_codename)) :
            for l in range(4):
                for m in range(4) :
                    print(par_codename[i]+"_"+bases[l]+"_"+bases[m]+" = "+str(par[i]),file=ofile)
            #impose continuity!
            if used[i] == False :
                output = impose_continuity(par_codename[i],i,par)
                vals = par_codename[i].split('_')
                names = []
                if output[0] == 'f1':                    
                    string = vals[0]+"_RLOW"
                    if len(vals) >= 3 :
                        if vals[2] == '33' or vals[2] == '55' :
                            string = string + '_' + vals[2]
                    names.append(string)
                    string = vals[0]+"_RHIGH"
                    if len(vals) >= 3 :
                        if vals[2] == '33' or vals[2] == '55' :
                            string = string + '_' + vals[2]
                    names.append(string)
                    string = vals[0]+"_BLOW"
                    if len(vals) >= 3 :
                        if vals[2] == '33' or vals[2] == '55' :
                            string = string + '_' + vals[2]
                    names.append(string)
                    string = vals[0]+"_BHIGH"
                    if len(vals) >= 3 :
                        if vals[2] == '33' or vals[2] == '55' :
                            string = string + '_' + vals[2]
                    names.append(string)
                    string = vals[0]+"_RCLOW"
                    if len(vals) >= 3 :
                        if vals[2] == '33' or vals[2] == '55' :
                            string = string + '_' + vals[2]
                    names.append(string)
                    string = vals[0]+"_RCHIGH"
                    if len(vals) >= 3 :
                        if vals[2] == '33' or vals[2] == '55' :
                            string = string + '_' + vals[2]
                    names.append(string)
                          
                elif output[0] == 'f4':
                    string = vals[0]+"_"+vals[1]+"_TS"
                    if len(vals) >= 4 :
                        if vals[3] == '33' or vals[3] == '55' :
                            string = string + '_' + vals[3]
                    names.append(string)
                    string = vals[0]+"_"+vals[1]+"_TC"
                    if len(vals) >= 4 :
                        if vals[3] == '33' or vals[3] == '55' :
                            string = string + '_' + vals[3]
                    names.append(string)
                    string = vals[0]+"_"+vals[1]+"_B"
                    if len(vals) >= 4 :
                        if vals[3] == '33' or vals[3] == '55' :
                            string = string + '_' + vals[3]
                    names.append(string)
                    
                elif output[0] == 'No':
                    continue
                    
                for k in range(len(names)) :
                    for l in range(4):
                        for m in range(4) :
                            print(names[k]+"_"+bases[l]+"_"+bases[m]+" = "+str(output[k+1]),file=ofile)                
                    
    ofile.close()
    ifile.close()
    
    return



#Compute mean and covariance with parameters par, by reweighting mean and covariance with parameters par0.
#The data used to compute ensamble avarages is stored inside the global list data.
#Data must be sampled with parameters par0 before reweighting.
#The reweighting procedure is for multivariate normal distributions.
def reweight_cov_and_mu(par,par0) :
    
    cov = np.zeros([dimension,dimension], dtype=float)    
    mu = np.zeros(dimension, dtype=float)    
    
    #<e^-DH>
    av_e_to_deltaH = np.zeros(dimension, dtype=float) 

    #data are sampled with par0 before calling reweighting, and stored globally.
    #This is to avoid sampling multiple times when unnecessary.
    #Fits well with what we want to do with oxDNA
    
    energy1 = []
    
    #update_rew_seq_dep_file(par)
    update_rew_seq_dep_file_ave(par)
    
    inp = []
    backend = []
    obs = []
    
    for i in range(Nreps) :
        with oxpy.Context():
            
            #read input script specifying sequence dependent file
            inp.append(oxpy.InputFile())
            
            inp[i].init_from_filename("./Rep"+str(i)+"/input2.an")
            #create analysys backend
            backend.append(oxpy.analysis.AnalysisBackend(inp[i]))
        
            obs.append(backend[i].config_info().observables)
            
            counts = -1
            
            while backend[i].read_next_configuration() :
                
                counts+=1
                
                if(counts < in_snap) :
                    continue
                a = float(obs[i][0].get_output_string(backend[i].conf_step).split()[0])
                energy1.append((Njuns+1)*20*a)
                
   # print(energy1)    
        
    ave_mu = []
    for i in range(len(ids)) :
        ave_mu.append(0.)
    
    
    #reweight mean
    for i in range(len(internal_coords)) :

        deltaH = 2*(energy1[i] - energy_sampled[i])
        
        print(-deltaH)
        
        for j in range(len(internal_coords[i])) :
            #if internal_coords[i][j] > -30 or internal_coords[i][j] < -10 : 
                mu[j] += internal_coords[i][j]*math.exp(-deltaH)
            
                av_e_to_deltaH[j] += math.exp(-deltaH)
        
        """
        ave_mu_i = 0.
        for j in range(len(internal_coords[i])) :
            ave_mu_i += internal_coords[i][j]/len(internal_coords[i])
        
        ave_mu += ave_mu_i*math.exp(-deltaH)
        """
    
    #print(av_e_to_deltaH)
    for i in range(len(mu)) :
        mu[i] = mu[i]*(1./av_e_to_deltaH[i])
        #ave_mu = ave_mu*(1./av_e_to_deltaH)
        ave_mu[i%len(ids)] += mu[i]/(1.*len(internal_coords[0])/(1.*len(ids)))
    
    """
    #reweight covariance
    for i in range(len(internal_coords)) :
        for j in range(len(internal_coords[i])) :
            for k in range(j+1,len(internal_coords[i])) :
     
                deltaH =  energy1[i] - energy_sampled[i]
                
                cov[j,k]+=(internal_coords[i][j]-mu[j])*(internal_coords[i][k]-mu[k])*math.exp(-deltaH)/(1.0*av_e_to_deltaH)
                
    for j in range(dimension) :
        for k in range(j+1,dimension) :
                cov[k,j] = cov[j,k]
    """
    
    return cov,mu,ave_mu



"""
#Compute Relative Entropy.
#It is a function of some parameters par;
#args are extra parameters:
#the mean and covariance of a given set of parameters par0, and par0 itself.
#If par = par0, the function uses the mean and covariance in args to compute the relative entropy,
#otherwise it estimates mean and covariance for par by reweighting par0.
#This function is built this way so that optimize.minimize uses reweighting to estimate Rel_entropy(par0+delta_par)
#e.g. when estimating the gradient and Hessian of the Relative entropy in par0
def Relative_entropy_wRew(par,par0):
    
    #cov0,mu0,par0 = decode_args(args)
    rew = False
    
    #check if par == par0 (1/1000 tolerance)
    for i in range(par_dimension) :
        if par[i] < par0[i]*(999./1000.) or par[i] > par0[i]*(1001./1000.) :
            rew = True
            break
    
    #if par != par0, then reweight mean and covariance from mean and covariance with par0
    if rew :
        cov,mu = reweight_cov_and_mu(par,par0)
    else :
        cov,mu = cov_sampled,mu_sampled
    

    #print(cov)
    #print(mu)
    
    #compute relative entropy
    
    S = 0.
    
    S += np.dot(cov,target_inv_cov).trace()
    
    S -= math.log(det_target_inv_cov*np.linalg.det(cov))+dimension
    
    delta_mu = mu - target_mu
    
    S += np.dot(np.dot(delta_mu.transpose(),target_inv_cov),delta_mu)
    
    S = S*0.5
    
    return S
"""

#Compute Relative Entropy. Mean only version
#The covariance of the target distribution is updated so that it coincides with the sampled covariance
def Relative_entropy_wRew_meanOnly(par,par0):
    
    #cov0,mu0,par0 = decode_args(args)
    rew = False
    
    #check if par == par0 (1/1000 tolerance)
    """
    for i in range(par_dimension) :
        if par[i] < par0[i]*(9999./10000.) or par[i] > par0[i]*(10001./10000.) :
            rew = True
            break

    
    #if par != par0, then reweight mean and covariance from mean and covariance with par0
    if rew :
        cov,mu = reweight_cov_and_mu(par,par0)
    else :
        cov,mu = cov_sampled,mu_sampled
    """
    

    cov,mu,ave_mu = reweight_cov_and_mu(par,par0)

    #print(cov)
    #print(mu)
    
    #compute relative entropy
    
    S = 0.
    
    #These 2 pieces are zero if cov = target_cov!
    #S += np.dot(cov,target_inv_cov).trace()
    #S -= math.log(det_target_inv_cov*np.linalg.det(cov))+dimension
    
    delta_mu = ave_mu
    delta_mu[0] += 12
    delta_mu[1] -= 36
    
    #Here we replace inv_target_cov with inv_cov
    
    #S = np.linalg.norm(delta_mu)
    S = np.linalg.norm(delta_mu)
    print(S)
    
    #S += np.dot(np.dot(delta_mu.transpose(),np.linalg.inv(cov)),delta_mu)
    
    #S = S*0.5
    
    return S
                


for i in range(Nreps) :

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
        store_internal_coord(traj,ids,in_j,fin_j,in_snap,True)
    else  :
        store_internal_coord(traj,ids,in_j,fin_j,in_snap,False)
 
#average coordinates and energy
ave_and_cov_stored()

print(mu_sampled)


inp = []
backend = []
obs = []

for i in range(Nreps) :
    with oxpy.Context():
        
        #read input script specifying sequence dependent file
        inp.append(oxpy.InputFile())
        
        inp[i].init_from_filename("./Rep"+str(i)+"/input1.an")
        #create analysys backend
        backend.append(oxpy.analysis.AnalysisBackend(inp[i]))
    
        obs.append(backend[i].config_info().observables)
        
        counts = -1
        
        while backend[i].read_next_configuration() :
            
            counts+=1
            
            if(counts < in_snap) :
                continue

            a = float(obs[i][0].get_output_string(backend[i].conf_step).split()[0])
            print(a)
            energy_sampled.append((Njuns+1)*20*a)
            
       
#set par0



ifile = open("oxDNA2_sequence_dependent_parameters_in.txt",'r')

par0 = []
#par0 = [2.8,0.5,0.5,1.0,0.7,1.8]

for line in ifile.readlines() :
    vals = line.split()
    
    #SD
    #for i in range(len(par_codename)) :
    #    if vals[0] == par_codename[i] :
    #        par0.append(float(vals[2]))
    
    #average
    
  
    if len(vals) == 0 :
            continue
    for i in range(len(par_codename)) :
        if vals[0] == par_codename[i]+'_A_A' :
            par0.append(float(vals[2]))            
            

ifile.close()


print(par0)



par = par0

sol = optimize.minimize(Relative_entropy_wRew_meanOnly,par,args=(par0),method='nelder-mead',options={'maxiter':5, 'eps':0.02})

print(sol)

#print(sol.x)

update_rew_seq_dep_file_ave(sol.x)

print("Step done!")

