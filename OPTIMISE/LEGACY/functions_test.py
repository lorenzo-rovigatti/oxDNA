#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:00:12 2023

@author: yqb22156
"""

import numpy as np
import math
import oxpy
import continuity_constraints
import config as cg
import get_cgdna_pars

#Parse config file (read parameters)
#Returns False if mandatory parameters are missing 
def read_config(cfile_name) :
    
    cfile = open(cfile_name,'r')
    
    checklist = np.zeros(11, dtype = int) #check if some mandatory parameters are missing
    
    for line in cfile.readlines() :
        vals = line.split()
        if len(vals) == 0:
            continue
        if vals[0][0] == '#':
            continue        
        #print(vals)
        
        #read sequence
        if(vals[0] == 'SEQ'):
            cg.seq = vals[1]
            cg.Njuns = len(vals[1])-1
            checklist[0] = 1
        
        #read initial and final junctions ids
        if(vals[0] == "IN_J") :
            cg.in_j = int(vals[1])
            checklist[1] = 1
            
        if(vals[0] == "FIN_J") :
            cg.fin_j = int(vals[1])
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
                checklist[8] = 1
            elif(vals[1]=="sd") :
                cg.ave = False
                checklist[9] = 1
                
        if(vals[0] == "MAXITER") :
            cg.miter = int(vals[1])
            checklist[10] = 1




    #CHECK AND PRINT
    
    if checklist[0] == 1:
        print("SEQUENCE: "+cg.seq)
        print("Njunctions: "+str(cg.Njuns))
    else:
        print("MANDATORY. Sequence missing from input file.")
        print("Usage:")
        print("SEQ seq")
        return False
    
    if checklist[1] == 1:
        print("IN JUNCTION: "+str(cg.in_j))
        print("Ignoring all junctions < "+ str(cg.in_j) +" in the optimisation.")
    else :
        print("OPTION. No in junction specified. Using default 0")
        print("Usage:")
        print("IN_J in_junction")
        
    if checklist[2] == 1:
        print("FIN JUNCTION: "+str(cg.fin_j))
        print("Ignoring all junctions > "+ str(cg.fin_j) +" in the optimisation.")
    else :
        print("OPTION. No fin junction specified. Using default value Njun = "+ str(cg.Njuns))
        print("Usage:")
        print("FIN_J fin_junction")
        
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
    else:
        print("OPTION. Ids covariance missing from config file.")
        print("Optimising ground state only")
        print("Usage:")
        print("ids_cov id1 id2 id3 ...")
        
    
    if checklist[3] == 1 and (checklist[4] == 1 or checklist[5] == 1):    #adjust this
        for i in range(len(cg.ids_gs)) :
            cg.ids.append(cg.ids_gs[i])
        for i in range(len(cg.ids_cov)) :
            if cg.ids_cov[i] in cg.ids :
                continue
            else :
                cg.ids.append(cg.ids_cov[i])
        cg.ids.sort()
        
        print("ALL IDS:")
        print(cg.ids)
        
    if checklist[3] == 1 and checklist[4] == 0:
        cg.ids = cg.ids_gs
        
    #generate gs(mu) and covariance. Sampled is initialised to 0, target is read from cgna+ 
    if checklist[1] == 1 and checklist[2] == 1 and checklist[3] == 1:
        cg.dimension = (cg.fin_j-cg.in_j+1)*(len(cg.ids))
        
        print("DIMENSION: " + str(cg.dimension))
        
        cg.mu_sampled = np.zeros(cg.dimension, dtype = float)
        cg.cov_sampled = np.zeros((cg.dimension,cg.dimension), dtype = float)
        
        if cg.ave == True and cg.diag == True:
            
            cg.target_mu, cg.target_cov = get_cgdna_pars.get_target_mean_and_covariance_diag_ave((cg.fin_j-cg.in_j+1), cg.ids)
            
        elif cg.ave == False and cg.diag == False:
            
            cg.target_mu, cg.target_cov = get_cgdna_pars.get_target_mean_and_covariance(cg.seq, cg.ids, cg.in_j, cg.Njuns-cg.fin_j)
            
        elif cg.ave == False and cg.diag == True:
            
            cg.target_mu, cg.target_cov = get_cgdna_pars.get_target_mean_and_covariance(cg.seq, cg.ids, cg.in_j, cg.Njuns-cg.fin_j)
            
            for i in(len(cg.target_cov)) :
                for j in(len(cg.target_cov)) :
                    if i != j :
                        cg.taget_cov[i,j] = 0.
                        
                        
        print("TARGET GS: see file target_gs.txt")
        
        ofile = open("target_gs.txt", 'w')
        for i in range(len(cg.target_cov)) :
            print(str(cg.target_mu[i]), file = ofile)
        ofile.close()
            
        print("TARGET COV: : see file target_cov.txt")
        ofile = open("target_cov.txt", 'w')
        for i in range(len(cg.target_cov)) :
            string = ""
            for j in range(len(cg.target_mu)) :            
                string += str(cg.target_cov[i,j]) + " "
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
        print("MAX number of ITERATIONS of the optimisation algorithm (Nelder-Mead): "+str(cg.miter))
    else :
        print("OPTION. No max iter specified. Using default value miter = "+ str(cg.miter))
        print("Usage:")
        print("MAXITER miter")
    
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
    
    for i in range(Ncoords) :
        cg.mu_sampled[i] = 0.
        for j in range(Ncoords) :
            cg.cov_sampled[i][j] = 0.
    
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            cg.mu_sampled[j] += cg.internal_coords[i][j]/Nsnaps
    
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            for z in range(j,Ncoords) :
                cg.cov_sampled[j][z] += (cg.internal_coords[i][j] - cg.mu_sampled[j])*(cg.internal_coords[i][z] - cg.mu_sampled[z])/Nsnaps
    
    for j in range(Ncoords) :
        for z in range(j+1,Ncoords) :
            cg.cov_sampled[z][j] = cg.cov_sampled[j][z]    
    
    return


#given a parameter, check if continuity constraints should be imposed
def impose_continuity(par_cname,p_id,pars) :
    vals = par_cname.split('_')
    auxiliars = []
    aux_values = []
    output = []
    
    f1 = False
    f2 = False
    f4 = False
    
    #f1   
    r0 = 0.
    rc = 0.
    a = 0.

    if vals[1] == 'R0' and vals[0] != 'FENE' and vals[0] != 'CRST':
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
        print(vals)
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
    
    
    #f2
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
       

    #f4
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
            #symmetries
            if vals[0] == 'STCK' or vals[0] == 'FENE':
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
                    
                elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'A' :
                    print(name+"_T_T"+" = "+str(par[i]),file=ofile)
                elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'T' :
                    print(name+"_A_A"+" = "+str(par[i]),file=ofile)
                    
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
                        print(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1]+" = "+str(output[k+1]),file=ofile)
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
                            
                    elif vals[0] == 'CRST' or vals[0] == 'HYDR':
                        print(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1]+" = "+str(output[k+1]),file=ofile)
                        if vals[len(vals)-2] !=  vals[len(vals)-1]  :
                            print(names[k]+"_"+vals[len(vals)-1]+"_"+vals[len(vals)-2]+" = "+str(output[k+1]),file=ofile)
                    else :
                        print(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1]+" = "+str(output[k+1]),file=ofile)
                    
                    
    ofile.close()
    ifile.close()
    
    return



#Compute mean and covariance with parameters par, by reweighting mean and covariance with parameters par0.
#The data used to compute ensamble avarages is stored inside the global list data.
#Data must be sampled with parameters par0 before reweighting.
#The reweighting procedure is for multivariate normal distributions.
def reweight_cov_and_mu(par,par0) :
    
    print(par)
    print(par0)
    
    cov = np.zeros([cg.dimension,cg.dimension], dtype=float)    
    mu = np.zeros(cg.dimension, dtype=float)    
    
    #<e^-DH>
    av_e_to_deltaH = np.zeros(cg.dimension, dtype=float) 

    #data are sampled with par0 before calling reweighting, and stored globally.
    #This is to avoid sampling multiple times when unnecessary.
    #Fits well with what we want to do with oxDNA
    
    energy1 = []
    
    #update_rew_seq_dep_file(par)
    if cg.ave:
        update_rew_seq_dep_file_ave(par)
    else:
        update_rew_seq_dep_file(par)
    
    inp = []
    backend = []
    obs = []
    
    read = False
    
    for i in range(cg.Nreps) :
        with oxpy.Context():
            
            #read input script specifying sequence dependent file
            inp.append(oxpy.InputFile())
            
            inp[i].init_from_filename("./Rep"+str(i)+"/input2.an")
            #create analysys backend
            backend.append(oxpy.analysis.AnalysisBackend(inp[i]))
        
            obs.append(backend[i].config_info().observables)
            
            counts = -1
            
            while 1==1 : 
                try:
                    read =  backend[i].read_next_configuration()
                except:
                    counts+=1
                    energy1.append(999)
                    print("Warning: exception in oxpy energy computation; reweighting. Rep "+str(i)+",conf" +str(counts))
                    continue
                if read == False :
                    break
                counts+=1
                
                if(counts < cg.in_snap) :
                    continue
                a = float(obs[i][0].get_output_string(backend[i].conf_step).split()[0])
                energy1.append((cg.Njuns+1)*20*a)
                
   # print(energy1)    
        
    ave_mu = []
    for i in range(len(cg.ids)) :
        ave_mu.append(0.)
    
    #reweight mean
    for i in range(len(cg.internal_coords)) :

        if energy1[i] != 999 and cg.energy_sampled[i] != 999:
            deltaH = (energy1[i] - cg.energy_sampled[i])
            
            #print(deltaH)
            
            for j in range(len(cg.internal_coords[i])) :

                    mu[j] += cg.internal_coords[i][j]*math.exp(-deltaH)
                    av_e_to_deltaH[j] += math.exp(-deltaH)
        
    
    for i in range(len(mu)) :
        mu[i] = mu[i]*(1./av_e_to_deltaH[i])
        ave_mu[i%len(cg.ids)] += mu[i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
    
    
    #reweight covariance
    for i in range(len(cg.internal_coords)) :
        if energy1[i] != 999 and cg.energy_sampled[i] != 999:
            for j in range(len(cg.internal_coords[i])) :
                    #if mode = diagonal (i.e. looking only at the diagonal of the full covariance matrix)
                    if cg.diag:
                        deltaH =  energy1[i] - cg.energy_sampled[i]                    
                        cov[j,j]+=(cg.internal_coords[i][j]-mu[j])*(cg.internal_coords[i][j]-mu[j])*math.exp(-deltaH)
                        #XXXXXX DO THIS FOR LONG RANGE STIFFNESS
                        if(j+1 < len(cg.internal_coords[i])) :
                            cov[j,j+1]+=(cg.internal_coords[i][j]-mu[j])*(cg.internal_coords[i][j+1]-mu[j+1])*math.exp(-deltaH)
                    else :
                        for z in range(j,len(cg.internal_coords[i])) :
                            deltaH =  energy1[i] - cg.energy_sampled[i]                    
                            cov[j,z]+=(cg.internal_coords[i][j]-mu[j])*(cg.internal_coords[i][z]-mu[z])*math.exp(-deltaH)
                    
                    
    for i in range(len(mu)):
        if cg.diag:
            cov[i,i] = cov[i,i]/(1.0*av_e_to_deltaH[i])
            if(i+1 < len(mu)):
               cov[i,i+1] = cov[i,i+1]/(1.0*av_e_to_deltaH[i])
            for j in range(i,len(mu)):   
                if i != j:
                    cov[j,i] = cov[i,j]
        else:
            for j in range(i,len(mu)):   
                cov[i,j] = cov[i,j]/(1.0*av_e_to_deltaH[i])
                if i != j:
                    cov[j,i] = cov[i,j]
    """
    print("mu: ")
    print(mu)
    
    print("ave_mu: ")
    print(ave_mu)
    
    print("cov: ")
    print(cov)
    """
    
    return cov,mu,ave_mu



# Computes Delta(q) for a given snapshot and q, where Delta is the displacement of the chosen set of internal coordinates 
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
        print("i: "+str(i))
        print(M_i*(25*a*a)/(m))
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
                    
    

#Compute Relative Entropy.
def Relative_entropy_wRew(par,par0):
    
    cov,mu,ave_mu = reweight_cov_and_mu(par,par0)
    
    #compute relative entropy
    
    S = 0. 
    """
    #compute averages
    ave_cov = []
    ave_target_cov = []
    ave_target_mu = []
    
    for i in range(len(cg.ids)) :
        ave_target_mu.append(0.)
        ave_cov.append(0.)
        ave_target_cov.append(0.)
    """
    
    ave_target_mu = np.zeros(len(cg.ids),dtype=float)
    ave_target_cov = np.zeros((len(cg.ids),len(cg.ids)),dtype=float)
    ave_cov = np.zeros((len(cg.ids),len(cg.ids)),dtype=float)

    
    for i in range(len(cov[0])) :
        ave_target_mu[i%len(cg.ids)] += cg.target_mu[i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
        if cg.diag == True:
            ave_cov[i%len(cg.ids),i%len(cg.ids)] += cov[i,i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
            ave_target_cov[i%len(cg.ids),i%len(cg.ids)] += cg.target_cov[i,i]/(1.*len(cg.internal_coords[0])/(1.*len(cg.ids)))
    
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
    
        S += 0.5*(np.dot(np.dot(ave_delta_mu.transpose(),np.linalg.inv(ave_target_cov)),ave_delta_mu))
    
    
    else :
        
        delta_mu = mu - cg.target_mu
        
        for i in range (len(delta_mu)) :
            if i%len(cg.ids) in cg.ids_gs:
                continue
            else:
                delta_mu[i] = 0       
    
        S += 0.5*(np.dot(np.dot(delta_mu.transpose(),np.linalg.inv(cg.target_cov)),delta_mu))
        
    print("Complete rew mu: ")
    print(mu) 
    
    print("Target mu: ")
    print(cg.target_mu)
    
    print("Rew ave_mu: ")
    print(ave_mu)
    
    print("Complete rew cov: ")
    print(cov)    
    
    print("Ave rew cov: ")
    print(ave_cov)    
    
    #STIFF
    
    if len(cg.ids_cov) > 0 :        
        
        if cg.ave == False:            
            #add symmetrised term of the gs part of S. If we are not adding the covariance part, the reweighted covariance does not appear in S, and we are optimising esclusively the gs
            invc = np.linalg.inv(cov)
            for i in range(len(invc)):
                for j in range(len(invc)):
                    if i!=j :
                        invc[i,j] = 0
            S += 0.5*(np.dot(np.dot(delta_mu.transpose(),invc),delta_mu))            
            
        else:
            S += 0.5*(np.dot(np.dot(ave_delta_mu.transpose(),np.linalg.inv(ave_cov)),ave_delta_mu))
            
        #reduce covariance for stiff part: remove from the covariance all coordinates we want to exclude from cov optimisation
        
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
            
            M = Stiff(reduced_cov,cg.Njuns-10)
            lb,lt = lps(M)
            
            
            print("Long range (m="+str(cg.Njuns-10)+"):")
            print("lb: "+str(lb))     
            print("lt: "+str(lt))
            
            #this term is a weighted sum of squared distances (l/lt+lt/l-2 = (l-lt)**2/(lt*l)). The weigth is as in the stiff term of the likelihood (see below).
            #The complete cost function is a hybrid: weighted sum of likelihoods for each sequence + squared distance for the average persistence lengths.
            S += 0.5*(lb/cg.target_lb + cg.target_lb/lb + lt/cg.target_lt + cg.target_lt/lt - 4)
                
        
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
            S += 0.5*(np.dot(np.linalg.inv(ave_reduced_cov),ave_reduced_target_cov).trace()+np.dot(np.linalg.inv(ave_reduced_target_cov),ave_reduced_cov).trace()-2*len(ave_reduced_cov))

        else:
            S += 0.5*(np.dot(np.linalg.inv(reduced_cov),reduced_target_cov).trace()+np.dot(np.linalg.inv(reduced_target_cov),reduced_cov).trace()-2*len(reduced_cov))
           
        print("COV-reduced rew cov: ")
        print(reduced_cov)   
        print("COV-reduced rew ave_cov: ")
        print(ave_reduced_cov) 
    
    print("S: "+str(S))
        
    return S