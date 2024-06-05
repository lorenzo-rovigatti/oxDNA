#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:00:12 2023

@author: yqb22156
"""

import numpy as np
import math
import oxpy
import copy
import continuity_constraints
import config_multi as cg
import get_cgdna_pars
from mpi4py import MPI

#Parse config file (read parameters)
#Returns False if mandatory parameters are missing 
def read_config(cfile_name) :
    
    cfile = open(cfile_name,'r')
    
    checklist = np.zeros(18, dtype = int) #check if some mandatory parameters are missing
    
    for line in cfile.readlines() :
        vals = line.split()
        if len(vals) == 0:
            continue
        if vals[0][0] == '#':
            continue        
        #print(vals)
        
        #read sequence
        if(vals[0] == 'SEQ'):
            cg.seq.append(vals[1])
            cg.Njuns.append(len(vals[1])-1)
            #cg.internal_coords.append([])
            #cg.energy_sampled.append([])
            checklist[0] = 1
        
        #read initial and final junctions ids
        if(vals[0] == "IN_J") :
            cg.inj = int(vals[1])
            checklist[1] = 1
            
        if(vals[0] == "J_FROM_END") :
            cg.jfe = int(vals[1])
            checklist[2] = 1
        
        #read which coordinates to optimise (ground state)
        if(vals[0] == 'IDS_GS') :
            for i in range(1,len(vals)) :
                cg.ids_gs.append(int(vals[i]))
            cg.ids_gs.sort()
            
            checklist[3] = 1
            
        #read which coordinates to optimise (covariance)
        if(vals[0] == 'IDS_COV') :
            for i in range(1,len(vals)) :
                cg.ids_cov.append(int(vals[i]))
            cg.ids_cov.sort()
            
            checklist[4] = 1
            
        #read which coordinates to optimise (large m elstic moduli). 
        #Currently we optimise the persistence length, not the moduli.
        if(vals[0] == 'IDS_LRS') :
            for i in range(1,len(vals)) :
                cg.ids_lrs.append(int(vals[i]))
            cg.ids_lrs.sort()
            
            checklist[5] = 1
        
        #read snapshots to discard (equilibration)
        if(vals[0] == "IN_SNAP") :
            cg.in_snap = int(vals[1]) 
            checklist[6] = 1
            
        #read how many replicas
        if(vals[0] == "REPS") :
            cg.Nreps = int(vals[1])
            checklist[7] = 1

        #read parameters to use for optimisation
        if(vals[0] == "cname") :
            #print(vals)
            if vals[2] == 'OPTIMISE' :
                cg.par_codename.append(vals[1])
                checklist[8] = 1
            else :
                cg.continuity_par_codenames.append(vals[1])
                cg.continuity_par_values.append(float(vals[3]))
            
        #options
            
        if(vals[0] == "MODE") :
            if(vals[1]=="ave") :
                cg.ave = True
                checklist[9] = 1
            elif(vals[1]=="sd") :
                cg.ave = False
                checklist[9] = 1
                
        if(vals[0] == "ALGO") :
            cg.algo = vals[1]
            checklist[10] = 1
            
        if(vals[0] == "MAXITER") :
            cg.miter = int(vals[1])
            checklist[11] = 1
            
        if(vals[0] == "NEVA") :
            cg.neva = int(vals[1])
            checklist[12] = 1
            
        if(vals[0] == "LBFGSB_EPS") :
            cg.LBFGSB_eps = float(vals[1])
            checklist[13] = 1            
            
        if(vals[0] == "LBFGSB_IPRINT") :
            cg.LBFGSB_iprint = int(vals[1])
            checklist[14] = 1
            
        if(vals[0] == "WEIGHT_GS") :
            cg.weight_gs = int(vals[1])
            checklist[15] = 1
            
        #optimise persistence lengths 0==true
        if(vals[0] == "LPS") :
            if int(vals[1]) != 0 :
                cg.opti_lp = True                
                checklist[16] = 1
                
        #symm stacking default = True
        if(vals[0] == "SYMM_STCK") :
            if int(vals[1]) != 0 :
                cg.symm_stck = True
            else :
                cg.symm_stck = False
                checklist[17] = 1
        

    #CHECK AND PRINT
    
    if checklist[0] == 1:
        
        cg.Nseq = len(cg.seq)
        
        for i in range(cg.Nseq):
                print("SEQUENCE " + str(i) + ": "+cg.seq[i])
                print("Njunctions Seq" + str(i) + ": "+str(cg.Njuns[i]))
    else:
        print("MANDATORY. Sequences missing from input file.")
        print("Usage:")
        print("SEQ seq")
        return False
    
    
    
    for i in range(cg.Nseq):
        cg.in_j.append(cg.inj)
    
    if checklist[1] == 1:
        print("IN JUNCTION: "+str(cg.inj))
        print("Ignoring all junctions < "+ str(cg.inj) +" in the optimisation.")
    else :
        print("OPTION. No in junction specified. Using default 0")
        print("Usage:")
        print("IN_J in_junction")
        
        
    for i in range(cg.Nseq):
        cg.fin_j.append(cg.Njuns[i]-cg.jfe-1)
        
    if checklist[2] == 1:
        for i in range(len(cg.fin_j)) :
            print("END JUNCTION Seq " + str(i) +": "+str(cg.fin_j[i]))
            print("Ignoring all junctions > "+ str(cg.fin_j[i]) +" for Seq "+str(i)+" in the optimisation.")
    else :
        print("OPTION. No junction from end specified. Using all junctions")
        print("Usage:")
        print("J_FROM_END junctions from end")
        
    if checklist[3] == 1:
        print("IDS GROUND STATE:")        
        print(cg.ids_gs)
        print("Optimising: ")
        for i in range(len(cg.ids_gs)) :
            if cg.ids_gs[i] == 0:
                print("Optimising buckle")
            if cg.ids_gs[i] == 1:
                print("Optimising propeller")
            if cg.ids_gs[i] == 2:
                print("Optimising opening")
            if cg.ids_gs[i] == 3:
                print("Optimising shear")
            if cg.ids_gs[i] == 4:
                print("Optimising stretch")
            if cg.ids_gs[i] == 5:
                print("Optimising stagger")
            if cg.ids_gs[i] == 6:
                print("Optimising tilt")
            if cg.ids_gs[i] == 7:
                print("Optimising roll")
            if cg.ids_gs[i] == 8:
                print("Optimising twist")
            if cg.ids_gs[i] == 9:
                print("Optimising shift")
            if cg.ids_gs[i] == 10:
                print("Optimising slide")
            if cg.ids_gs[i] == 11:
                print("Optimising rise")
    else:
        print("MANDATORY. Ids ground state missing from config file.")
        print("Usage:")
        print("IDS_GS id1 id2 id3 ...")
        return False
    
    
    if checklist[4] == 1:
        print("IDS COVARIANCE:")
        print(cg.ids_cov)
        print("Optimising: ")
        for i in range(len(cg.ids_cov)) :
            if cg.ids_cov[i] == 0:
                print("Optimising buckle")
            if cg.ids_cov[i] == 1:
                print("Optimising propeller")
            if cg.ids_cov[i] == 2:
                print("Optimising opening")
            if cg.ids_cov[i] == 3:
                print("Optimising shear")
            if cg.ids_cov[i] == 4:
                print("Optimising stretch")
            if cg.ids_cov[i] == 5:
                print("Optimising stagger")
            if cg.ids_cov[i] == 6:
                print("Optimising tilt")
            if cg.ids_cov[i] == 7:
                print("Optimising roll")
            if cg.ids_cov[i] == 8:
                print("Optimising twist")
            if cg.ids_cov[i] == 9:
                print("Optimising shift")
            if cg.ids_cov[i] == 10:
                print("Optimising slide")
            if cg.ids_cov[i] == 11:
                print("Optimising rise")
    else:
        print("OPTION. Ids covariance missing from config file.")
        print("Optimising ground state only")
        print("Usage:")
        print("ids_cov id1 id2 id3 ...")
        
        
    if checklist[5] == 1:
        print("IDS LONG RANGE STIFFNESS (q=0):")
        print(cg.ids_lrs)
        print("Optimising: ")
        for i in range(len(cg.ids_cov)) :
            if cg.ids_lrs[i] == 0:
                print("Optimising buckle")
            if cg.ids_lrs[i] == 1:
                print("Optimising propeller")
            if cg.ids_lrs[i] == 2:
                print("Optimising opening")
            if cg.ids_lrs[i] == 3:
                print("Optimising shear")
            if cg.ids_lrs[i] == 4:
                print("Optimising stretch")
            if cg.ids_lrs[i] == 5:
                print("Optimising stagger")
            if cg.ids_lrs[i] == 6:
                print("Optimising tilt")
            if cg.ids_lrs[i] == 7:
                print("Optimising roll")
            if cg.ids_lrs[i] == 8:
                print("Optimising twist")
            if cg.ids_lrs[i] == 9:
                print("Optimising shift")
            if cg.ids_lrs[i] == 10:
                print("Optimising slide")
            if cg.ids_lrs[i] == 11:
                print("Optimising rise")
                
        

    if cg.opti_lp == True :
        print("Optimising persistence lengths (i.e. Ar and C)")
    else:
        print("Neglecting persistence lengths (i.e. Ar and C)")
        print("Usage:")
        print("LPS i")
        print("i = 0 for optimising lps.")
        
    
    
    #collect all ids:
        
    if checklist[3] == 1 :    
        for i in range(len(cg.ids_gs)) :
            cg.ids.append(cg.ids_gs[i])
            
    if checklist[4]==1 :
        for i in range(len(cg.ids_cov)) :
            if cg.ids_cov[i] in cg.ids :
                continue
            else :
                cg.ids.append(cg.ids_cov[i])
                
    if checklist[16] == 1 :
        for i in range(len(cg.ids_inter_rot)) :
            if cg.ids_inter_rot[i] in cg.ids :
                continue
            else :
                cg.ids.append(cg.ids_inter_rot[i])    


    cg.ids.sort()
    
    #initialise deltas    
    
    print("ALL IDS:")
    print(cg.ids)
    
    #generate gs(mu) and covariance. Sampled is initialised to 0, target is read from cgna+ 
    if checklist[3] == 1:        
        
        for i in range(cg.Nseq) :
            cg.dimension.append((cg.fin_j[i]-cg.in_j[i]+1)*(len(cg.ids)))
        
        cg.mu_sampled = np.zeros(cg.dimension[cg.seq_id], dtype = float)
        cg.cov_sampled = np.zeros((cg.dimension[cg.seq_id],cg.dimension[cg.seq_id]), dtype = float)
        cg.cov0_sampled = np.zeros((cg.dimension[cg.seq_id],cg.dimension[cg.seq_id]), dtype = float)
        
        for i in range(cg.Nseq) :
        
            print("DIMENSION Seq "+str(i)+": " + str(cg.dimension[i]))       

            
            if cg.ave == True and cg.diag == True:
                
                tm, tcov = get_cgdna_pars.get_target_mean_and_covariance_diag_ave((cg.fin_j[i]-cg.in_j[i]+1), cg.ids)
                
                cg.target_mu.append(tm)
                cg.target_cov.append(tcov)
                
            elif cg.ave == False and cg.diag == False:
                
                tm, tcov = get_cgdna_pars.get_target_mean_and_covariance(cg.seq[i], cg.ids, cg.inj, cg.jfe)
                
                cg.target_mu.append(tm)
                cg.target_cov.append(tcov)
                
            elif cg.ave == False and cg.diag == True:
                
                tm, tcov = get_cgdna_pars.get_target_mean_and_covariance(cg.seq[i], cg.ids, cg.inj, cg.jfe)
                
                cg.target_mu.append(tm)
                cg.target_cov.append(tcov)
                
                for k in range(len(cg.target_cov[i])) :
                    for l in range(len(cg.target_cov[i])) :
                        if k != l :
                            cg.target_cov[i][k,l] = 0.
                        
                        
        print("TARGET GS: see file target_gs.txt")
        
        for l in range(cg.Nseq) :
            ofile = open("target_gs_Seq"+str(l)+".txt", 'w')
            for i in range(len(cg.target_mu[l])) :
                print(str(cg.target_mu[l][i]), file = ofile)
            ofile.close()
            
        print("TARGET COV: : see file target_cov.txt")
        for l in range(cg.Nseq) :
            ofile = open("target_cov_Seq"+str(l)+".txt",'w')
            for i in range(len(cg.target_cov[l])) :
                string = ""
                for j in range(len(cg.target_cov[l])) :            
                    string += str(cg.target_cov[l][i,j]) + " "
                print(string,file=ofile)    
            ofile.close()
            
            
    if checklist[6] == 1:
        print("INITIAL SNAPSHOT: "+str(cg.in_snap))
        print("Ignoring all sampled snapshopts < "+ str(cg.in_snap) +" in the optimisation.")
    else :
        print("OPTION. No in snapshot specified. Using all snapshots (is the trajectory equilibrated?).")
        print("Usage:")
        print("IN_SNAP in_snap")
        
    if checklist[7] == 1:
        print("NUMBER OF REPETITIONS: "+str(cg.Nreps))
    else :
        print("OPTION. No Nreps specified. Running only one simulation replica.")
        print("Usage:")
        print("IN_SNAP in_snap")
        
        
    if checklist[8] == 1:
        
        cg.par_dimension = len(cg.par_codename)
        for i in range(cg.par_dimension) :
            cg.used.append(False)
        
        print("PARAMETERS used in the optimisation: ")
        print("OPTIMISE - "+str(cg.par_dimension)+":")
        for par in cg.par_codename :
            print(par)      
        print("AUXILIARY (for continuity) - "+str(len(cg.continuity_par_codenames))+":")
        if len(cg.continuity_par_codenames) == 0:
            print("None")
        else:
            for i in range(len(cg.continuity_par_codenames)) :
                print(cg.continuity_par_codenames[i]+" = "+str(cg.continuity_par_values[i]))      
       
        
    else :
        print("MANDATORY. No parameters for optimisation specified.")
        print("Usage. For optimisation parameters:")
        print("cname par_name1 OPTIMISE")
        print("cname par_name2 OPTIMISE")
        print("...")
        print("Usage. For auxiliary parameters (for imposing continuity constraints):")
        print("cname par_name1 CONTINUITY value")
        print("cname par_name2 CONTINUITY value")
        print("...")
        
        return False
    
    
    if checklist[9] == 1:
       if cg.ave == True:
           print("MODE: AVERAGE")
       else :
           print("MODE: SD")
    else:
        print("OPTION. Mode not specified.")
        print("Using default mode average.")
        print("Usange: ")
        print("MODE ave/sd")
        
        
    if checklist[10] == 1:
        print("Optimisation algorithm selected: "+str(cg.algo))
    else :
        print("OPTION. No optimisation algorithm specified. Using default: "+ str(cg.algo))
        print("Usage:")
        print("ALGO algorithm")
        print("Available algorithms:")
        print("See scipy.optimise.minimise for documentation:")
        print("---1---")
        print("algorithm = nelder-mead")
        print("Options:")
        print("NEVA neva")
        print("---2---")
        print("algorithm = L-BFGS-B")
        print("Options:")
        print("NEVA neva")
        print("LBFGSB_EPS eps (default 0.01)")
        print("LBFGSB_IPRINT iprint (default 1)")        
        
        
    if checklist[11] == 1:
        print("MAX number of ITERATIONS of the optimisation algorithm: "+str(cg.miter))
    else :
        print("OPTION. No max iter specified. Using default value miter = "+ str(cg.miter))
        print("Note: it might be that MAX number of function evaluations is used, instead")
        print("Usage:")
        print("MAXITER miter")
        
    if checklist[15] == 1:
        print("Using weight for GS part of the relative entropy: "+str(cg.weight_gs))
    else :
        print("OPTION. No weight for GS part of the relative entropy specified. Using default weight_gs = "+ str(cg.weight_gs))
        print("Usage:")
        print("WEIGHT_GS weight_gs")
        
    if checklist[12] == 1:
        print("MAX number of function evaluations of the optimisation algorithm: "+str(cg.neva))
    else :
        print("OPTION. No function evaluations specified. Using default value miter = "+ str(cg.neva))
        print("Note: it might be that MAX number of iterations is used, instead")
        print("Usage:")
        print("NEVA neva")

    if cg.algo == "L-BFGS-B" :
        if checklist[13] == 1:
            print("Eps for L-BFGS-B algorithm: "+str(cg.LBFGSB_eps))
        else :
            print("OPTION. No eps specified for L-BFGS-B algorithm. Using default value eps = "+ str(cg.LBFGSB_eps))
            print("Note: it might be that MAX number of iterations is used, instead")
            print("Usage:")
            print("LBFGSB_EPS eps")
            
        if checklist[14] == 1:
            print("iprint for L-BFGS-B algorithm: "+str(cg.LBFGSB_iprint))
        else :
            print("OPTION. No iprint specified for L-BFGS-B algorithm. Using default value iprint = "+ str(cg.LBFGSB_iprint))
            print("Note: it might be that MAX number of iterations is used, instead")
            print("Usage:")
            print("LBFGSB_iprint iprint")
            
    if checklist[17] == 1:
        print("Breaking stacking symmetry (AA/TT only)")
    else :
        print("OPTION. Using symmetric stacking.")
        print("If you want to break the AA/TT symmetry,")
        print("Usage:")
        print("SYMM_STCK 1")
            
    

    
    return True


#given a junction trajectory (read_oxdna_trajectory_standard_order), store specific internal coordinates in
#global variable internal_coords
def store_internal_coord(traj,ids,in_j,fin_j,in_snap,overwrite=True) :
       
    coords = []
    
    Nsnaps = len(traj)
    Njuns = len(traj[0])
    
    alpha = 5.*math.pi/180 #angles in cgna are in radiants/5
    
    for i in range(in_snap,Nsnaps) :

        coord = []
        for j in range(Njuns) :
            #print("snap:"+str(i)+" junn:"+str(j))
            if j < in_j or j > fin_j :
                continue
            for z in range(len(ids)) :
                if ids[z] == 0 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[0]*alpha)
                elif ids[z] == 1 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[1]*alpha)
                elif ids[z] == 2 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[2]*alpha)
                elif ids[z] == 3 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[0])
                elif ids[z] == 4 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[1])
                elif ids[z] == 5 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[2])
                    
                elif ids[z] == 6 :
                    coord.append(traj[i][j].inter_coord.rot[0]*alpha)
                elif ids[z] == 7 :
                    coord.append(traj[i][j].inter_coord.rot[1]*alpha)
                elif ids[z] == 8 :
                    coord.append(traj[i][j].inter_coord.rot[2]*alpha)
                elif ids[z] == 9 :
                    coord.append(traj[i][j].inter_coord.tran[0])
                elif ids[z] == 10 :
                    coord.append(traj[i][j].inter_coord.tran[1])
                elif ids[z] == 11 :
                    coord.append(traj[i][j].inter_coord.tran[2])
                    
        if overwrite == False :
            cg.internal_coords.append(coord)
        else :
            coords.append(coord)
           
    if overwrite == True :
        cg.internal_coords = coords
                    
    return

def ave_and_cov_stored() :
      

    Nsnaps = len(cg.internal_coords)
    Ncoords = len(cg.internal_coords[0])
    
    Nave = Nsnaps
    
    for i in range(Nsnaps) :
        if 998.9 < cg.energy_sampled[i] and cg.energy_sampled[i] > 999.1 :
            Nave -= 1
    
    for i in range(Ncoords) :
        cg.mu_sampled[i] = 0.
        for j in range(Ncoords) :
            cg.cov_sampled[i][j] = 0.
            cg.cov0_sampled[i][j] = 0.
       
    for i in range(Nsnaps) :
        if 998.9 < cg.energy_sampled[i] and cg.energy_sampled[i] > 999.1 :
            continue
        for j in range(Ncoords) :
            cg.mu_sampled[j] += cg.internal_coords[i][j]/Nave
    
    for i in range(Nsnaps) :
        if 998.9 < cg.energy_sampled[i] and cg.energy_sampled[i] > 999.1 :
            continue
        for j in range(Ncoords) :
            for z in range(j,Ncoords) :
                cg.cov_sampled[j][z] += (cg.internal_coords[i][j] - cg.mu_sampled[j])*(cg.internal_coords[i][z] - cg.mu_sampled[z])/Nave
                cg.cov0_sampled[j][z] += cg.internal_coords[i][j]*cg.internal_coords[i][z]/Nave
    
    for j in range(Ncoords) :
        for z in range(j+1,Ncoords) :
            cg.cov_sampled[z][j] = cg.cov_sampled[j][z]
            cg.cov0_sampled[z][j] = cg.cov0_sampled[j][z] 
    
    return


def compute_deltas(First) :
    #seq = ""
    cg.Deltas = np.zeros((len(cg.ids),4,4),dtype=float)
    if First:
        cg.bsteps_counts = np.zeros((4,4), dtype=int)
    for i in range(len(cg.Deltas)):
        for j in range(4) :
            for k in range(4):
                cg.Deltas[i][j][k] = 0.
    
    for i in range(cg.in_j[cg.seq_id], cg.fin_j[cg.seq_id]) :
        #seq+=cg.seq[cg.seq_id]
        #print("Seq "+str(cg.seq_id)+ " index " +str(i) + " " + cg.seq[cg.seq_id][i])
        if First:
            cg.bsteps_counts[cg.base_to_id(cg.seq[cg.seq_id][i])][cg.base_to_id(cg.seq[cg.seq_id][i+1])]+=1
            
        for j in range(len(cg.ids)) :
            delta = cg.mu_curr[(i-cg.in_j[cg.seq_id])*len(cg.ids)+j]-cg.target_mu[cg.seq_id][(i-cg.in_j[cg.seq_id])*len(cg.ids)+j]
            cg.Deltas[j][cg.base_to_id(cg.seq[cg.seq_id][i])][cg.base_to_id(cg.seq[cg.seq_id][i+1])]+=delta
            #print("Delta"+str(j)+ " " +str(cg.seq[cg.seq_id][i])+str(cg.seq[cg.seq_id][i+1])+ " ("+str(cg.base_to_id(cg.seq[cg.seq_id][i]))+str(cg.base_to_id(cg.seq[cg.seq_id][i+1]))+") " +str(delta))
            
    


#given a parameter, check if continuity constraints should be imposed
def impose_continuity(par_cname,p_id,pars) :
    vals = par_cname.split('_')
    auxiliars = []
    aux_values = []
    output = []
    
    f1 = False
    f2 = False
    f4 = False
    
    #print("Continuity pars")
    
    #this is the scaling factor of the HYDRO and STCK (e.g. HYDR_A_T and STCK_G_A)
    if len(vals) == 3 and (vals[0] == "HYDR" or vals[0] == "STCK") and (vals[1] in cg.bases) and (vals[2] in cg.bases) :
        output.append('No')
        return output
    
    
    #check if the parameter is in f1 modulation and impose continuity of the modulation and its first derivative
    #i.e. from r0, A and Rc derive rl,rh,bl,bh,rcl,rch (check continuity_constraints.py)
    
    r0 = 0.
    rc = 0.
    a = 0.  
    
    #print("searching")

    if vals[1] == 'R0' and vals[0] != 'FENE' and vals[0] != 'CRST': #FENE and CRST have different modulations or the radial part
        f1 = True
        r0 = pars[p_id]
        auxiliars.append('A')
        auxiliars.append('RC')
    elif vals[1] == 'RC' :
        f1 = True
        rc = pars[p_id]
        auxiliars.append('R0')
        auxiliars.append('A')
    elif vals[1] == 'A':   
        f1 = True
        a = pars[p_id]
        auxiliars.append('R0')
        auxiliars.append('RC')
    if f1 :
        #print("f1")
        #check if we are also optimising one of the auxiliary parameters
        for i in range(len(auxiliars)) : 
            found = False
            aux_pname = vals[0]+'_'+auxiliars[i]
            for j in range(2,len(vals)) :
                aux_pname = aux_pname + '_' + vals[j]
            print(aux_pname)
            for j in range(len(cg.continuity_par_codenames)) :
                if aux_pname == cg.continuity_par_codenames[j] :
                    aux_values.append(cg.continuity_par_values[j])
                    found = True
            if found == False :
                for j in range(len(cg.par_codename)) :
                    if aux_pname == cg.par_codename[j] :
                        aux_values.append(pars[j])
                        cg.used[j] = True
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
    
    
    #check if the parameter is in f2 modulation
    if len(vals) >= 2:
        if vals[1] == 'R0' and vals[0] == "CRST" :
            print("f2")
            f2 = True
            r0 = pars[p_id]

    if f2 :

        rc,rl,rh,bl,bh,rcl,rch,continuous = continuity_constraints.continuity_f2_fixed_width(r0) #fix rc so that the amplitude (and depth) is not affected 
        
        output.append('f2')
        output.append(rc)
        output.append(rl)
        output.append(rh)        
        output.append(bl)
        output.append(bh)
        output.append(rcl)
        output.append(rch)
        
        return output
       

    #check if the parameter is in f4 modulation
    a = 0.

    if len(vals) > 2 :
        if vals[2] == 'A' and (vals[1] != 'HYDR' or vals[1] != 'STCK') and  vals[0] != 'FENE':
            f4 = True
            a = pars[p_id]

    if f4 :

        dts,dtc,b= continuity_constraints.continuity_f4(a)
        
        output.append('f4')
        output.append(dts)
        output.append(dtc)
        output.append(b)
        
        return output
    
    if f1 == False and f2 == False and f4 == False:
        output.append('No')
        return output
         
   
def build_initial_simplex_for_nm(x0,up_bnds,low_bnds) :
    in_simplex = []
    
    x = copy.deepcopy(x0)
    for i in range(len(x)) :
        x[i]*=1.0
    in_simplex.append(x)
    
    
    for i in range(len(x0)) :
        x = copy.deepcopy(x0)
                
        vals1 = cg.par_codename[i].split('_')
        
        if vals1[0] == "FENE" and vals1[1] == "R0" :    #stricter for FENE_R0 (+-3%), to avoid problems with the FENE potential

            found = False

            for k in range(len(cg.ids)):
                if cg.ids[k] == 8:
                    found = True
                    if cg.Deltas[k][cg.base_to_id(vals1[2])][cg.base_to_id(vals1[3])]< 0:
        
                        if x[i]*1.02 < up_bnds[i] : 
                            x[i]*=1.02
                        else :
                            x[i] = up_bnds[i] - 0.001
                    else :
                        if x[i]*0.98 > low_bnds[i] : 
                            x[i]*=0.98
                        else :
                            x[i] = low_bnds[i] + 0.001     
            if found == False :
                if x[i]*1.02 < up_bnds[i] : 
                    x[i]*=1.02
                else :
                    x[i] = up_bnds[i] - 0.001
                
        
        elif vals1[0] == "FENE" and vals1[1] == "DELTA" :    #stricter for FENE_DELTA (+-5%), to avoid problems with the FENE potential
            
            if x[i]*1.02 < up_bnds[i] : 
                x[i]*=1.02
            else :
                x[i] = up_bnds[i] - 0.001
                
        
            
        elif vals1[0] == "STCK" and vals1[1] == "R0" :
            
            found = False
                
            for k in range(len(cg.ids)):
                if cg.ids[k] == 11:
                    found = True
                    if cg.Deltas[k][cg.base_to_id(vals1[2])][cg.base_to_id(vals1[3])] < 0:
        
                        if x[i]*1.1 < up_bnds[i] : 
                            x[i]*=1.1
                        else :
                            x[i] = up_bnds[i] - 0.001
                    else :
                        if x[i]*0.9 > low_bnds[i] : 
                            x[i]*=0.9
                        else :
                            x[i] = low_bnds[i] + 0.001 
                            
            if found == False:
                if x[i]*1.1 < up_bnds[i] : 
                    x[i]*=1.1
                else :
                    x[i] = up_bnds[i] - 0.001
            
        elif vals1[0] == "STCK" and vals1[2] == "T0" and abs(x[i]) < 0.01:
            
            if x[i]+0.2 < up_bnds[i] : 
                x[i]+=0.2
            else :
                x[i] = up_bnds[i] - 0.001
            
        elif vals1[0] == "CRST" and vals1[1] == "THETA4" and vals1[2] == "T0" :   #stricter for FENE_R0 (+-3%), to avoid problems with the FENE potential
            

            found = False
            for k in range(len(cg.ids)):
                if cg.ids[k] == 7:
                    found = True
                    if vals1[3] == '55':
                        if cg.Deltas[k][cg.base_to_id(vals1[4])][3-cg.base_to_id(vals1[5])] < 0:
        
                            if x[i]+0.2 < up_bnds[i] : 
                                x[i]+=0.2
                            else :
                                x[i] = up_bnds[i] - 0.001
                        else :
                            if x[i]-0.2 > low_bnds[i] : 
                                x[i]-=0.2
                            else :
                                x[i] = low_bnds[i] + 0.001
                    elif vals1[3] == '33':
                        if cg.Deltas[k][3-cg.base_to_id(vals1[5])][cg.base_to_id(vals1[4])] < 0:
        
                            if x[i]+0.2 < up_bnds[i] : 
                                x[i]+=0.2
                            else :
                                x[i] = up_bnds[i] - 0.001
                        else :
                            if x[i]-0.2 > low_bnds[i] : 
                                x[i]-=0.2
                            else :
                                x[i] = low_bnds[i] + 0.001
            if found == False:
                if x[i]+0.2 < up_bnds[i] : 
                    x[i]+=0.2
                else :
                    x[i] = up_bnds[i] - 0.001                
                        
            
                
        elif vals1[0] == "CRST" and vals1[1] == "THETA4" and vals1[2] == "A" :   #stricter for FENE_R0 (+-3%), to avoid problems with the FENE potential
            
           if x[i]*3.0 < up_bnds[i] : 
               x[i]*=3.0
           else :
               x[i] = up_bnds[i] - 0.001
               
               
        elif vals1[0] == "STCK" and (vals1[1] == "THETA4" or vals1[1] == "THETA5") and vals1[2] == "A":
            
            
            found = False
            for k in range(len(cg.ids)):
                if cg.ids[k] == 1:
                    found = True
                    if cg.Deltas[k][cg.base_to_id(vals1[2])][cg.base_to_id(vals1[3])] > 0:
        
                        if x[i]*1.5 < up_bnds[i] : 
                            x[i]*=1.5
                        else :
                            x[i] = up_bnds[i] - 0.001
                    else :
                        if x[i]*0.75 > low_bnds[i] : 
                            x[i]*=0.75
                        else :
                            x[i] = low_bnds[i] + 0.001
            if found == False:
                if x[i]*1.5 < up_bnds[i] : 
                    x[i]*=1.5
                else :
                    x[i] = up_bnds[i] - 0.001
                            
        elif vals1[0] == "HYDR" and vals1[1] == "THETA4" and vals1[2] == "A":
            
            found = False
            for k in range(len(cg.ids)):
                if cg.ids[k] == 1:
                    found = True
                    delta = 0
                    for m in range(4):
                        for n in range(4):
                            delta += cg.Deltas[k][m][n]/16.
                    
                    if delta < 0:
        
                        if x[i]*1.5 < up_bnds[i] : 
                            x[i]*=1.5
                        else :
                            x[i] = up_bnds[i] - 0.001
                    else :
                        if x[i]*0.75 > low_bnds[i] : 
                            x[i]*=0.75
                        else :
                            x[i] = low_bnds[i] + 0.001 
                            
            if found == False :
                if x[i]*1.5 < up_bnds[i] : 
                    x[i]*=1.5
                else :
                    x[i] = up_bnds[i] - 0.001
                
            
        else :                    
            if x[i]*1.5 < up_bnds[i] : 
                x[i]*=1.5
            else :
                x[i] = up_bnds[i] - 0.001
                
        in_simplex.append(x)
        
        
    ofile = open("Deltas_and_in_simplexes.txt", 'a')
    
    print("ITER "+ str(cg.Niter), file=ofile)
    print("DELTAS", file=ofile)
    print("",file=ofile)
    for i in range(len(cg.Deltas)):
        if cg.ids[i] == 0:
            print("buckle")
        if cg.ids[i] == 1:
            print("propeller")
        if cg.ids[i] == 2:
            print("opening")
        if cg.ids[i] == 3:
            print("shear")
        if cg.ids[i] == 4:
            print("stretch")
        if cg.ids[i] == 5:
            print("stagger")
        if cg.ids[i] == 6:
            print("tilt")
        if cg.ids[i] == 7:
            print("roll")
        if cg.ids[i] == 8:
            print("twist")
        if cg.ids[i] == 9:
            print("shift")
        if cg.ids[i] == 10:
            print("slide")
        if cg.ids[i] == 11:
            print("rise")
            
        for j in range(4):
            string = ""
            for k in range(4):
                string+= " "+str(cg.Deltas[i][j][k])
            print(string, file=ofile)
        print("",file=ofile)
        
    print("INITIAL SIMPLEX", file=ofile)
    
    i = 0
    string = ""
    for j in range(len(in_simplex[i])):
        string += " " + " {:.3f}".format(in_simplex[i][j])
    print(string, file=ofile)
    
    string = ""
    for i in range(1, len(in_simplex)):
        string += " " + " {:.3f}".format(in_simplex[i][i-1])
    print(string, file=ofile)
            
    print("",file=ofile)
    
    ofile.close()        
    
    
    return in_simplex
             
        
    
#update sequence dependent file with new values (par) of the optimisation parameters + continuity constraints and symmetries.
def update_rew_seq_dep_file_ave(par) :
    
    ofile = open('oxDNA_sequence_dependent_parameters_tmp.txt','w')
    ifile = open('oxDNA_sequence_dependent_parameters.txt','r')    #file with fixed sequence dependent parameters
    
    #continuity for f1
    
    if len(par) != len(cg.par_codename) :
        print("Something is not right with the parameters! Check codename file.")
    
    else :
        for line in ifile.readlines() :
            print(line.strip('\n'), file=ofile)
            
        print('\n', file=ofile)
        
        for i in range(len(cg.par_codename)) :
            
            if cg.par_codename[i] == "FENE_DELTA":
                print(cg.par_codename[i]+" = "+str(par[i]),file=ofile)
                print(cg.par_codename[i]+"2 = "+str(par[i]*par[i]),file=ofile)
                print('\n', file=ofile)
                continue
            
            for l in range(4):
                for m in range(4) :
                    print(cg.par_codename[i]+"_"+cg.bases[l]+"_"+cg.bases[m]+" = "+str(par[i]),file=ofile)
            
            #impose continuity!
            if cg.used[i] == False :
                output = impose_continuity(cg.par_codename[i],i,par)
                vals = cg.par_codename[i].split('_')
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
                    
                elif output[0] == 'f2':
                    string = vals[0]+"_RC"
                    if len(vals) >= 3 :
                        if vals[2] == '33' or vals[2] == '55' :
                            string = string + '_' + vals[2]
                    names.append(string)
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
                            print(names[k]+"_"+cg.bases[l]+"_"+cg.bases[m]+" = "+str(output[k+1]),file=ofile)
                print('\n', file=ofile)
                    
    ofile.close()
    ifile.close()
    
    return



#update sequence dependent file with new values (par) of the optimisation parameters + continuity constraints and symmetries. SD version
def update_rew_seq_dep_file(par) :
    
    print("WARNING: F2 SD continuity is not implemented yet!!!!")
    
    ofile = open('oxDNA_sequence_dependent_parameters_tmp.txt','w')
    ifile = open('oxDNA_sequence_dependent_parameters.txt','r')    #file with fixed sequence dependent parameters
    
    #continuity for f1
    
    if len(par) != len(cg.par_codename) :
        print("Something is not right with the parameters! Check codename file.")
    
    else :
        for line in ifile.readlines() :
            print(line.strip('\n'), file=ofile)
            
        print('\n', file=ofile)
        
        for i in range(len(cg.par_codename)) :
            print(cg.par_codename[i]+" = "+str(par[i]),file=ofile)
            vals = cg.par_codename[i].split('_')
            name = vals[0]
            for k in range(1,len(vals)-2) :
                name = name + "_"+vals[k]
                
            if vals[1] == 'DELTA' :
                print(name+"2"+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1]+" = "+str(par[i]*par[i]),file=ofile)
            
            #symmetries
            if vals[0] == 'STCK':
                
                """
                id1 = cg.base_to_id(vals[len(vals)-3])
                id2 = cg.base_to_id(vals[len(vals)-2])
                id3 = cg.base_to_id(vals[len(vals)-1])
                
                cb1 = cg.bases[3-id3]
                cb2 = cg.bases[3-id2]
                cb3 = cg.bases[3-id1]
                
                
                if vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'G' :
                    print(name+"_C_C"+" = "+str(par[i]),file=ofile)                      
                elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'C' :
                    print(name+"_G_G"+" = "+str(par[i]),file=ofile)
                    
                elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'A' :
                    print(name+"_T_C"+" = "+str(par[i]),file=ofile)
                elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'C' :
                    print(name+"_G_A"+" = "+str(par[i]),file=ofile)
                    
                elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'G' :
                    print(name+"_C_T"+" = "+str(par[i]),file=ofile)
                elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'T' :
                    print(name+"_A_G"+" = "+str(par[i]),file=ofile)
                    
                elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'G' :
                    print(name+"_C_A"+" = "+str(par[i]),file=ofile)
                elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'A' :
                    print(name+"_T_G"+" = "+str(par[i]),file=ofile)
                    
                elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'T' :
                    print(name+"_A_C"+" = "+str(par[i]),file=ofile)
                elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'C' :
                    print(name+"_G_T"+" = "+str(par[i]),file=ofile)
                
                if cg.symm_stck:
                    if vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'A' :
                        print(name+"_T_T"+" = "+str(par[i]),file=ofile)
                    elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'T' :
                        print(name+"_A_A"+" = "+str(par[i]),file=ofile)
                
                if cb1 != vals[len(vals)-3] and cb2 != vals[len(vals)-2] and cb3 != vals[len(vals)-1]:
                    print(name+"_"+cb1+"_"+cb2+"_"+cb3+" = "+str(par[i]),file=ofile)
                """
            #symmetries
            elif vals[0] == 'FENE':
                
                if vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'G' :
                    print(name+"_C_C"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_C_C"+" = "+str(par[i]*par[i]),file=ofile)                      
                elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'C' :
                    print(name+"_G_G"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_G_G"+" = "+str(par[i]*par[i]),file=ofile)      
                    
                elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'A' :
                    print(name+"_T_C"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_T_C"+" = "+str(par[i]*par[i]),file=ofile)     
                elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'C' :
                    print(name+"_G_A"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_G_A"+" = "+str(par[i]*par[i]),file=ofile)      
                    
                elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'G' :
                    print(name+"_C_T"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_C_T"+" = "+str(par[i]*par[i]),file=ofile)      
                elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'T' :
                    print(name+"_A_G"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_A_G"+" = "+str(par[i]*par[i]),file=ofile)    
                    
                elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'G' :
                    print(name+"_C_A"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_C_A"+" = "+str(par[i]*par[i]),file=ofile)  
                elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'A' :
                    print(name+"_T_G"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_T_G"+" = "+str(par[i]*par[i]),file=ofile)     
                    
                elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'T' :
                    print(name+"_A_C"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_A_C"+" = "+str(par[i]*par[i]),file=ofile)    
                elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'C' :
                    print(name+"_G_T"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_G_T"+" = "+str(par[i]*par[i]),file=ofile)  
                    
                elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'A' :
                    print(name+"_T_T"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_T_T"+" = "+str(par[i]*par[i]),file=ofile)      
                elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'T' :
                    print(name+"_A_A"+" = "+str(par[i]),file=ofile)
                    if vals[1] == 'DELTA' :       
                        print(name+"2"+"_A_A"+" = "+str(par[i]*par[i]),file=ofile)      
                    
            elif vals[0] == 'CRST' or vals[0] == 'HYDR':
                if vals[len(vals)-2] !=  vals[len(vals)-1]  :
                    print(name+"_"+vals[len(vals)-1]+"_"+vals[len(vals)-2]+" = "+str(par[i]),file=ofile)
            #impose continuity!
            if cg.used[i] == False :
                output = impose_continuity(cg.par_codename[i],i,par)
                #vals = par_codename[i].split('_')
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
                    
                #update seq dep file + impose symmetries
                for k in range(len(names)) :
                    if vals[0] == 'STCK' :
                        print(names[k]+"_"+vals[len(vals)-3]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1]+" = "+str(output[k+1]),file=ofile)
                        """
                        if vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'G' :
                            print(names[k]+"_C_C"+" = "+str(output[k+1]),file=ofile)
                        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'C' :
                            print(names[k]+"_G_G"+" = "+str(output[k+1]),file=ofile)
                            
                        elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'A' :
                            print(names[k]+"_T_C"+" = "+str(output[k+1]),file=ofile)
                        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'C' :
                            print(names[k]+"_G_A"+" = "+str(output[k+1]),file=ofile)
                            
                        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'G' :
                            print(names[k]+"_C_T"+" = "+str(output[k+1]),file=ofile)
                        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'T' :
                            print(names[k]+"_A_G"+" = "+str(output[k+1]),file=ofile)
                            
                        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'G' :
                            print(names[k]+"_C_A"+" = "+str(output[k+1]),file=ofile)
                        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'A' :
                            print(names[k]+"_T_G"+" = "+str(output[k+1]),file=ofile)
                            
                        elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'T' :
                            print(names[k]+"_A_C"+" = "+str(output[k+1]),file=ofile)
                        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'C' :
                            print(names[k]+"_G_T"+" = "+str(output[k+1]),file=ofile)
                            
                        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'A' :
                            print(names[k]+"_T_T"+" = "+str(output[k+1]),file=ofile)
                        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'T' :
                            print(names[k]+"_A_A"+" = "+str(output[k+1]),file=ofile)
                        """
                            
                    elif vals[0] == 'CRST' or vals[0] == 'HYDR':
                        print(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1]+" = "+str(output[k+1]),file=ofile)
                        if vals[len(vals)-2] !=  vals[len(vals)-1]  :
                            print(names[k]+"_"+vals[len(vals)-1]+"_"+vals[len(vals)-2]+" = "+str(output[k+1]),file=ofile)
                    else :
                        print(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1]+" = "+str(output[k+1]),file=ofile)
                    
                    
    ofile.close()
    ifile.close()
    
    return



# Computes stiffness matrix(m) (see Eq. 10 and 12 of Enrico's 2017 paper) from covariance
def Stiff(cov,m) :
    
    M = np.zeros((3,3),dtype = float)

    a=0.34
    
    n = len(cg.ids_cov)
    
    #print("n: " +str(n))
    
    Nbp = int(len(cov)/n)
    
    #print("Nbp: " +str(Nbp))
    
    for i in range(Nbp-m) :
        M_i = np.zeros((3,3),dtype = float)
        for j in range(i,i+m) :
            for k in range(i,i+m) :
                #print("i: " +str(i)+" j: " +str(j))
                row = -1
                M_jk = np.zeros((3,3),dtype = float)
                for l in range(n) :
                    if cg.ids_cov[l] in cg.ids_inter_rot :
                        row += 1
                        col = -1
                        for z in range(n) :
                            if cg.ids_cov[z] in cg.ids_inter_rot :
                                col+=1                        
                                #print(str(row)+" " +str(col))
                                M_jk[row,col] = cov[j*n+l,k*n+z]/(25*a*a)
                #if i == Nbp-m-1 :
                #    print("j: "+str(j))
                #    print(M_jk*(25*a*a))
                M_i += M_jk
            
        #M_i[0,1] = 0.
        #M_i[1,0] = 0.
        #M_i[0,2] = 0.
        #M_i[2,0] = 0.
        #print("i: "+str(i))
        #print(M_i*(25*a*a)/(m))
        M += np.linalg.inv(M_i)*(m/a)/(Nbp-m)
    
    return M

#Compute persistence lengths from stiffness matrix
def lps(M) :
    
    A1 = M[0,0]
    A2 = M[1,1]
    C = M[2,2]
    G = M[1,2]
    
    lb = 2*A1*(A2-G*G/C)/(A1+A2-G*G/C)
    lt = 2*C*(1-G*G/A2/C)
    
    return lb,lt
                    

def print_matrix(M):
    
    for i in range(len(M)):
        string = ""
        for j in range(len(M[i])) :
            string += str(M[i][j]) + " "
        print(string)
        
    return


def callbackF(par) :
    
    cg.Niter += 1
    
    ofile = open("parameters_v_iter.txt", 'a')
    
    string = str(cg.Niter)
    
    for i in range(len(par)):
        
        string = string + " " + str(par[i])
    
    print(string, file=ofile)
    
    ofile.close()
    
    ofile = open("S_v_iter.txt", 'a')
    
    if cg.Niter == 1:
        print('#niter = number of iterations, neva = number of S evaluations, S = cost function')
        print('#niter neva S', file=ofile)        
    
    string = str(cg.Niter) + " " + str(cg.curr_feva) + " "+ str(cg.S_curr)
    
    print(string, file=ofile)
    
    ofile.close()
    
    
    return



#Compute Relative Entropy.
def Relative_entropy_wRew(par,stop,par0):
    
    stop[0]=cg.comm.bcast(stop[0], root=0)  #this is used to stop all processes (the while loop in main cycle)
                                            #at the end of optimisation stop[0] is set to 1 and communicated to all processes                                          
    S = 0.
    
    if stop[0] == 0 :
                    
        if cg.rank == 0:
            
            frac = []
            for i in range(len(par)):
                frac.append(par[i]/par0[i])            
            
            """
            for k in range(len(par)):
                par[k] *= par0[k]    
            """
            print("parameters")
            #print(par)
            print(["{0:0.3f}".format(i) for i in par])
            #print(par0)
            #print(["{0:0.3f}".format(i) for i in par0])
            print("fraction (par/par0):")
            #print(frac)
            print(["{0:0.3f}".format(i) for i in frac])
            
            #update parameters file (only once, at rank 0)
            if cg.ave:
                update_rew_seq_dep_file_ave(par)
            else:
                update_rew_seq_dep_file(par)
                
        #print("We are here 0 rank " +str(cg.rank))
        #bcast par from rank 0 (where optimisation is performed) to other cpus
        par=cg.comm.bcast(par, root=0)
        #print("We are there 0 rank " +str(cg.rank))
        #print("communicated par")
        #print(par)
        
        
        
        ######################
        #### REWEIGHT GS AND COV
        ######################
        
        l = cg.seq_id
        rep = cg.rep_id     
        

        #compute new energy (with par)         
        energy1 = []
                 
        read = False
        
        cg.stop_flag = 0

        with oxpy.Context():      
            
             
             #print("We are here -1 rank " +str(cg.rank))
             #read input script specifying sequence dependent file
             inp = oxpy.InputFile()
             inp.init_from_filename("./Seq"+str(l)+"/Rep"+str(rep)+"/input2.an")
             
             #backend = oxpy.analysis.AnalysisBackend(inp)
             
             try :
                 #create analysys backend
                 backend = oxpy.analysis.AnalysisBackend(inp)
             except :
                 print("Could not start oxpy. Throwing a stop flag.")
                 cg.stop_flag = 1

             if cg.stop_flag == 0 :
                 obs = backend.config_info().observables
                 
                 #print("We are there -1 rank " +str(cg.rank))
                 
                 counts = -1
                 
                 #print("We are here rank " +str(cg.rank))
                 
                 while 1==1 : 
                     try:
                         read =  backend.read_next_configuration()
                     except:
                         counts+=1
                         energy1.append(999)
                         print("Warning: exception in oxpy energy computation; reweighting. Seq "+str(l)+", Rep "+str(rep)+", conf" +str(counts))
                         continue
                     if read == False :
                         break
                     counts+=1
                     
                     if(counts < cg.in_snap) :
                         continue
                     a = float(obs[0].get_output_string(backend.conf_step).split()[0])
                     #print(counts)
                     if math.isnan( a ) or (abs((cg.Njuns[l]+1)*20*a-cg.energy_sampled[counts-cg.in_snap])>70):    #avoid nans and overflows
                         energy1.append(999)
                         #print("We are here 0 rank " +str(cg.rank))
                     else :
                         energy1.append((cg.Njuns[l]+1)*20*a)
                         
        counts_disc = 0
                         
        if cg.stop_flag == 0 :              
           for k in range(len(energy1)) :
               if energy1[k] < 999.01 and energy1[k] > 998.99 :
                   counts_disc += 1
           
        #sanity check: if more than half of the configurations are discarded, then the parameters are too extreme. 
        #Return S = 1000000, works with Nelder-Mead.
           
           if counts_disc >= len(energy1)*0.5 :
                print("Warning: too many discarded configurations")
                cg.stop_flag = 1
                         
                         
        #if something goes horribly wrong, return S = 10^6
        cg.stop_flag = cg.comm.allreduce(cg.stop_flag,op=MPI.SUM)
        
        if cg.stop_flag > 0:
            if cg.rank == 0:
                print("Warning!: Got a stop flag in the reweighting. Returning S = 10^6")
            return 1000000
                     
        #print("We are there rank " +str(cg.rank))
        counts_disc = 0                     

        #reweight for rep rep and seq l
        
        cov = np.zeros([cg.dimension[l],cg.dimension[l]], dtype=float)    
        mu = np.zeros(cg.dimension[l], dtype=float)    
         
        #<e^-DH>
        av_e_to_deltaH = np.zeros(cg.dimension[l], dtype=float) 
         
        #reweight mean for seq l rep rep
        for i in range(len(cg.internal_coords)) :
     
             if (energy1[i] > 999.01 or energy1[i] < 998.99) and (cg.energy_sampled[i] > 999.01 or cg.energy_sampled[i] < 998.99):
                 deltaH = (energy1[i] - cg.energy_sampled[i])
                 if math.isnan( deltaH ) :
                     print("rank "+ str(cg.rank) + " " + str(i) + " " + str(deltaH))                     
                 
                 for j in range(len(cg.internal_coords[i])) :
     
                         mu[j] += cg.internal_coords[i][j]*math.exp(-deltaH)
                         
                         if math.isnan(math.exp(-deltaH)) :
                             print("Exp is nan: delta = "+str(deltaH))
                         
                         av_e_to_deltaH[j] += math.exp(-deltaH)
             
        #reduce gs to seq leader and compute rew gs for seq l (i.e. sum over reps)
        
        #print("We are here 1 rank " +str(cg.rank))
        
        mu = cg.comm_seq.reduce(mu,op=MPI.SUM, root=0)
        av_e_to_deltaH = cg.comm_seq.allreduce(av_e_to_deltaH,op=MPI.SUM)
        
        #print("We are there 1 rank " +str(cg.rank))
        
        if cg.rank_seq == 0: #same as cg.rank in cg.leaders
        

            for i in range(len(mu)) :
                 mu[i] = mu[i]*(1./av_e_to_deltaH[i])
                 
        
        mu = cg.comm_seq.bcast(mu, root=0)   #mu, for each sequence, is now the average over all reps
        #av_e_to_deltaH = cg.comm_seq.bcast(av_e_to_deltaH, root=0)   #same as mu: sum over all reps
        
        #reweight covariance
        for i in range(len(cg.internal_coords)) :
             
             if (energy1[i] > 999.01 or energy1[i] < 998.99) and (cg.energy_sampled[i] > 999.01 or cg.energy_sampled[i] < 998.99):
                 
                 for j in range(len(cg.internal_coords[i])) :
                         for z in range(j,len(cg.internal_coords[i])) :
                             
                             deltaH =  energy1[i] - cg.energy_sampled[i]                    
                             cov[j,z]+=(cg.internal_coords[i][j]-mu[j])*(cg.internal_coords[i][z]-mu[z])*math.exp(-deltaH)
        
        #reduce cov to seq leader
        
        #print("We are here 2 rank " +str(cg.rank))
        
        cov = cg.comm_seq.reduce(cov,op=MPI.SUM, root=0)
        
        #print("We are there 2 rank " +str(cg.rank))
          
        if cg.rank_seq == 0:         
            for i in range(len(mu)):
     
                 for j in range(i,len(mu)):   
                     cov[i,j] = cov[i,j]/(1.0*av_e_to_deltaH[i])
                     if i != j:
                         cov[j,i] = cov[i,j]

        ####################
        ### COMPUTE RELAIVE ENTROPY
        ###################

        #we use seq leaders for this (i.e. we are done with summing over reps of same sequence)        
        
        if cg.rank_seq == 0:
            
            cg.curr_mu = mu
            
        
            """
            print("LENGTHS")
            print(len(mu[0]))
            print(len(cg.target_mu[0]))
            print(len(cov[0][0]))
            print(len(cg.target_cov[0][0]))
            """
            
            
            #initialise average observables
            ave_C = 0
            ave_Ar = 0
            ave_G = 0
            ave_lb = 0
            ave_lt = 0
            
            opti_lp = cg.opti_lp # change to true for lp 
            
            for z in range(len(cg.ids_inter_rot)) :
                if cg.ids_inter_rot[z] in cg.ids_cov :
                    continue
                else :
                    print("Warning: Not all inter rotations are used for covariance optimisation. Cannot tune stiffness matrix (persistence length).")
                    opti_lp = False
                    
            #compute relative entropy
            
            S = 0. 
            
        
            
            print("Seq: "+str(l))    
        
            #initialise    
        
            ave_target_mu = np.zeros(len(cg.ids),dtype=float)
            ave_target_cov = np.zeros((len(cg.ids),len(cg.ids)),dtype=float)
            ave_cov = np.zeros((len(cg.ids),len(cg.ids)),dtype=float)
        
            ave_mu = []
        
            for i in range(len(cg.ids)) :
                 ave_mu.append(0.)
                 
            for i in range(len(mu)) :
                 ave_mu[i%len(cg.ids)] += mu[i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
            
            for i in range(len(cov[0])) :
                ave_target_mu[i%len(cg.ids)] += cg.target_mu[l][i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
                if cg.diag == True:
                    ave_cov[i%len(cg.ids),i%len(cg.ids)] += cov[i,i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
                    ave_target_cov[i%len(cg.ids),i%len(cg.ids)] += cg.target_cov[l][i,i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
            
            #GS
            
            delta_mu = np.zeros(len(mu), dtype = float)
            ave_delta_mu = np.zeros(len(cg.ids), dtype = float)
            
            if cg.ave == True:
                
                ave_delta_mu = ave_mu - ave_target_mu
                
                for i in range (len(ave_delta_mu)) :
                    if cg.ids[i] in cg.ids_gs:
                        continue
                    else:
                        ave_delta_mu[i] = 0       
            
                S += 0.5*(np.dot(np.dot(ave_delta_mu.transpose(),np.linalg.inv(ave_target_cov)),ave_delta_mu))*cg.weight_gs
            
            
            else :
                
                delta_mu = mu - cg.target_mu[l]
                
                for i in range (len(delta_mu)) :
                    if cg.ids[i%len(cg.ids)] in cg.ids_gs:
                        continue
                    else:
                        delta_mu[i] = 0       
            
                S += 0.5*(np.dot(np.dot(delta_mu.transpose(),np.linalg.inv(cg.target_cov[l])),delta_mu))*cg.weight_gs
                
            print("SEQUEUCE "+str(l))
            print("###############")
            print("Complete rew mu: ")
            #print(mu)
            string = ""
            counts_coo = -1
            for coo in mu:
                counts_coo+=1
                if counts_coo % len(cg.ids) == 0 :
                    string = ""
                string += " {:.3f}".format(coo)
                if counts_coo % len(cg.ids) == len(cg.ids)-1 :
                    print(string) 
            
            print("Delta (target - rew mu): ")
            #print(cg.target_mu[l]-mu)
            
            string = ""
            counts_coo = -1
            for coo in mu:
                counts_coo+=1
                if counts_coo % len(cg.ids) == 0 :
                    string = ""
                string += " " + " {:.3f}".format((cg.target_mu[l][counts_coo]-coo))
                if counts_coo % len(cg.ids) == len(cg.ids)-1 :
                    print(string) 
                    
                        
            print("Delta Ratio Delta/Delta0 (Delta0 = sampled): ")
            #print(cg.target_mu[l]-mu)
            
            string = ""
            counts_coo = -1
            for coo in mu:
                counts_coo+=1
                if counts_coo % len(cg.ids) == 0 :
                    string = ""
                string += " " + " {:.3f}".format((cg.target_mu[l][counts_coo]-coo)/(cg.target_mu[l][counts_coo]-cg.mu_sampled[counts_coo]))
                if counts_coo % len(cg.ids) == len(cg.ids)-1 :
                    print(string) 

            
            
            print("Rew ave_mu: ")
            print(ave_mu)
            
            print("Complete rew cov: ")
            print(cov)    
            
            print("Ave rew cov: ")
            print(ave_cov)    
            
            #COVARIANCE + STIFFNESS MATRIX (LONG RANGE)
            
            if len(cg.ids_cov) > 0 :        
                
    
                if cg.ave == False:            
                    #add symmetric term of the gs part of S. If we are not adding the covariance part, the reweighted covariance does not appear in S, and we are optimising esclusively the gs
                    invc = np.linalg.inv(cov)
                    for i in range(len(invc)):
                        for j in range(len(invc)):
                            if i!=j :
                                invc[i,j] = 0
                    S += 0.5*(np.dot(np.dot(delta_mu.transpose(),invc),delta_mu))*cg.weight_gs      
                    
                else:
                    S += 0.5*(np.dot(np.dot(ave_delta_mu.transpose(),np.linalg.inv(ave_cov)),ave_delta_mu))*cg.weight_gs
    
                    
                #reduce covariance for stiff part: remove from the covariance all coordinates we want to exclude from cov optimisation
                
                red_n = (cg.fin_j[l]-cg.in_j[l]+1)*(len(cg.ids)-len(cg.ids_cov))
                
                reduced_cov = np.zeros((cg.dimension[l]-red_n,cg.dimension[l]-red_n),dtype=float)
                reduced_target_cov = np.zeros((cg.dimension[l]-red_n,cg.dimension[l]-red_n),dtype=float)
                
                ave_reduced_cov = np.zeros((len(cg.ids_cov),len(cg.ids_cov)),dtype=float)
                ave_reduced_target_cov = np.zeros((len(cg.ids_cov),len(cg.ids_cov)),dtype=float)
                
                nrow=-1
                
                for i in range(cg.dimension[l]):
                    if cg.ids[i%len(cg.ids)]  in cg.ids_cov:
            
                        nrow+=1
                        ncol=-1
                        for j in range(cg.dimension[l]):
                            if cg.ids[j%len(cg.ids)]  in cg.ids_cov:
                                ncol+=1
                                reduced_cov[nrow,ncol] = cov[i,j]
                                reduced_target_cov[nrow,ncol] = cg.target_cov[l][i,j]
                nrow=-1          
                for i in range(len(cg.ids)) :
                    if cg.ids[i] in cg.ids_cov :
                        nrow+=1
                        ncol=-1
                        for j in range(len(cg.ids)):
                            if cg.ids[j]  in cg.ids_cov:
                                ncol+=1
                                ave_reduced_cov[nrow,ncol] = ave_cov[i,j]
                                ave_reduced_target_cov[nrow,ncol] = ave_target_cov[i,j]
                                
                                
                                
                #COMPUTE PERSISTENCE LENGTHS AND STIFFNESS MODULI        
                ###WORKS AS INTENDED ONLY IF IDS_COV = 6,7,8 (i.e. inter rotations)
                ###TO EXTEND THIS WE HAVE TO REDUCE THE COVARIANCE SO THAT WE HAVE ONLY THE INTER ROTATIONS    
                
                if opti_lp :                    
                    
                    if(cg.Njuns[l] != cg.Njuns[0]) :
                        print("Warning: not all sequences are of the same length.") 
                        print("m in the computation of asymptotic stffness is not the same for all sequences (Njuns - 4)")
                    
                    
                    M = Stiff(reduced_cov,cg.Njuns[l]-cg.inj-cg.jfe -3)                                
                    
                    lb,lt = lps(M)
                    
                    C = M[2,2]
                    Ar = M[1,1]
                    G = M[1,2]
                    
                    ave_C = C
                    ave_Ar = Ar
                    ave_G = G
                    ave_lb = lb
                    ave_lt = lt
                    
                    
                    print("Long range (m="+str(cg.Njuns[l]-cg.inj-cg.jfe -3)+"):")
                    print("lb: "+str(lb))     
                    print("lt: "+str(lt))
                    
                    print("C: "+str(C))     
                    print("Ar: "+str(Ar))
                        
                
                #diagonal cov
                
                for i in range(len(ave_reduced_cov)):
                    for j in range(len(ave_reduced_cov)):
                        if i != j:
                            ave_reduced_cov[i,j] = 0.
                            
                for i in range(len(reduced_cov)):
                    for j in range(len(reduced_cov)):
                        if i != j:
                            reduced_cov[i,j] = 0.
                
                if cg.ave == True:
                    
                    #ave_reduced_target_cov[1,1]*=1.1 #rescale Ar for tuning bending persistence length
                    
                    S += 0.5*(np.dot(np.linalg.inv(ave_reduced_cov),ave_reduced_target_cov).trace()+np.dot(np.linalg.inv(ave_reduced_target_cov),ave_reduced_cov).trace()-2*len(ave_reduced_cov))
        
                else:
                    S += 0.5*(np.dot(np.linalg.inv(reduced_cov),reduced_target_cov).trace()+np.dot(np.linalg.inv(reduced_target_cov),reduced_cov).trace()-2*len(reduced_cov))
                   
                print("COV-reduced rew cov: ")
                #print_matrix(reduced_cov)
                print(reduced_cov)  
                print("COV-reduced rew ave_cov: ")
                print(ave_reduced_cov) 
                
            if opti_lp and len(cg.ids_cov) > 0 :
                #this term is a weighted sum of squared distances (l/lt+lt/l-2 = (l-lt)**2/(lt*l)). The weigth is as in the stiff term of the likelihood (see below).
                #The complete cost function is a hybrid: weighted sum of likelihoods for each sequence + squared distance for the average persistence lengths.
                #S += 0.5*(lb/cg.target_lb + cg.target_lb/lb + lt/cg.target_lt + cg.target_lt/lt - 4)
                
                
                #reduce moduli to rank 0 and compute average 
                
                ave_C = cg.comm_leaders.reduce(ave_C,op=MPI.SUM, root=0)
                ave_Ar = cg.comm_leaders.reduce(ave_Ar,op=MPI.SUM, root=0)
                ave_G = cg.comm_leaders.reduce(ave_G,op=MPI.SUM, root=0)
                ave_lb = cg.comm_leaders.reduce(ave_lb,op=MPI.SUM, root=0)
                ave_lt = cg.comm_leaders.reduce(ave_lt,op=MPI.SUM, root=0)
                
                ave_C /= cg.Nseq
                ave_Ar /= cg.Nseq
                ave_G /= cg.Nseq                
                ave_lb /= cg.Nseq
                ave_lt /= cg.Nseq
                
                
                print("Average long range stiffness:")
                print("lb: "+str(ave_lb))     
                print("lt: "+str(ave_lt))
                
                print("C: "+str(ave_C))     
                print("Ar: "+str(ave_Ar))
                
                scaling_factor = cg.Nseq
                
                """
                if cg.ave == False:
                    for l in range(cg.Nseq) :
                        scaling_factor = scaling_factor* (cg.Njuns[l]-cg.inj-cg.jfe -3)
                """
                
                S += 0.5*(ave_C/cg.target_C + cg.target_C/ave_C + ave_Ar/cg.target_Ar + cg.target_Ar/ave_Ar - 4)*scaling_factor    
        
    
            #finally, reduce S to rank 0, which runs the optimisation        
            
            
            print("seq, rep: " + str(cg.seq_id) + ", " + str(cg.rep_id) + ". S: " + str(S))
            
            if S > 100000 or S < 0 :
                print("S overflow. Setting it to 10^6")
                S = 1000000
            
            S = cg.comm_leaders.reduce(S,op=MPI.SUM, root=0)
            
            if cg.rank == 0 :                
                print("tot S: "+str(S))
                cg.S_curr = S
                cg.curr_feva += 1
                
    
    return S