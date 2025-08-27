#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:00:12 2023

@author: yqb22156
"""


import numpy as np
import copy
import math
import oxpy
import config_multi as cg
import SantaLucia
import functions_multi
from estimate_Tm import get_Tm_and_width
import pysed
from mpi4py import MPI

#Parse config file (read parameters)
#Returns False if mandatory parameters are missing 

def read_config(cfile_name) :
    
    cfile = open(cfile_name,'r')
    
    if cg.first_step_flag == False:
        tfile = open("in_Ts.txt",'r')
        Ts_lines = tfile.readlines()
    Ts_index = -1
    
    checklist = np.zeros(12, dtype = int) #check if some mandatory parameters are missing
    
    Seq_counter = -1
    
    for line in cfile.readlines() :
        vals = line.split()
        if len(vals) == 0:
            continue
        if vals[0][0] == '#':
            continue        
        #print(vals)
        
        #read sequence
        if(vals[0] == 'SEQ'):
            
            Seq_counter += 1            
            Ts_index += 1
            
            if len(vals) != 9 :
                print("Invalid SEQ sintax.")
                print("Usage:")
                print("SEQ seq box_size Cs simT0 Nts DT opfile wfile")
                print("box_size = box_size")
                print("Cs = ssalt concentration")
                print("simT0 = simulation temperature first step")
                print("NTs = number of rew temperatures")
                print("DT = DT between rew temperatures")
                print("opfile = order param file")
                print("wfile = weight_file")
                return False            
            
            checklist[0] = 1
            
            cg.seq.append(vals[1])
            cg.Njuns.append(len(vals[1])-1)
            #cg.energy_sampled.append([])
            #cg.hbs_sampled.append([])
            
            Ct = 2./pow(float(vals[2]),3.)*2.6868
            Cs = float(vals[3])
            
            cg.Ct.append(Ct)            
            cg.Cs.append(Cs)


            SLt = SantaLucia.melting_temperature(vals[1],Ct,Cs)
            SimT = 0
            if cg.first_step_flag == True:
                SimT = float(vals[4])     
            else:
                SimT = float(Ts_lines[Ts_index])

            cg.target_mTs.append(SLt)
            
            NT = int(vals[5])
            DT = float(vals[6])
            
            cg.simTs.append(SimT)            
            cg.n_Ts.append(NT)


            """
            #Does the range of selected temperatures cover both simT and target T?
            #if not, adjust DT

            if Delta > NT*DT+16 :
                
                DT_new = (NT*Delta+16)/NT
                
                print("Warning: DT is too small to cover simulated T and target mT.")
                print("SimT - Target T: " + str(SimT-SLt))
                print("Enlarging DT to: " + str(DT_new))
            """

            #DT = DT_new
            cg.DTs.append(DT)


            #create T range around T0

            #T0 = (SimT-SLt)/2.
            T0 = SimT
            
            rew_Ts = []
                        
            lr = 0
            ur = 0
            
            if (NT)%2 == 0 :
                lr = -int((NT)/2-1)
                ur = int((NT)/2+1)
            else :
                lr = -int((NT-1)/2)
                ur = int((NT-1)/2+1)
                
            
            offset = 0.
            if DT*lr+T0<0 :  #oxpy does not like T < 0C and T > 100C
                offset = -(DT*lr+T0-0.5)
            if DT*ur+T0>100 :
                offset = - (DT*ur+T0-100) - 0.5
            
            for l in range(lr,ur) :
                rew_Ts.append( DT*l+T0  + offset)
                
                
            cg.rew_Ts.append(rew_Ts)
                
            
            #read weights

            wfile = open("Seq"+str(Seq_counter)+"/Rep0/"+vals[8],'r')
            cg.fin_wfiles.append("Seq"+str(Seq_counter)+"/"+vals[8])
          
            weights = []
            
            for line in wfile.readlines() :
                vals = line.split()
                if len(vals) == 0:
                    continue
                if vals[0][0] == '#':
                    continue  
                
                weights.append(float(vals[1]))
                
            cg.weights.append(weights)
            wfile.close()
            
            
            #initialise op_hists
            
            s_op_hist_un = []
            s_op_hist_b = []
            s_op_hist_rew_Ts = []
            
            for l in range(len(weights)) :
                
                s_op_hist_un.append(0)
                s_op_hist_b.append(0)
                
                s_op_hist_rew_T = []
                
                for m in range(len(rew_Ts)):
                
                    s_op_hist_rew_T.append(0)
                    
                s_op_hist_rew_Ts.append(s_op_hist_rew_T)
                
                        
            cg.sampled_op_histo_b.append(s_op_hist_b)
            cg.sampled_op_histo_un.append(s_op_hist_un)
            cg.sampled_op_hist_all_Ts.append(s_op_hist_rew_Ts)      
            
        
        
        #read snapshots to discard (equilibration)
        if(vals[0] == "IN_SNAP") :
            cg.in_snap = int(vals[1]) 
            checklist[1] = 1
            
        #read how many replicas
        if(vals[0] == "REPS") :
            cg.Nreps = int(vals[1])
            checklist[2] = 1

        #read parameters to use for optimisation
        if(vals[0] == "cname") :
            #print(vals)
            if vals[2] == 'OPTIMISE' :
                cg.par_codename.append(vals[1])
                checklist[3] = 1
            else :
                cg.continuity_par_codenames.append(vals[1])
                cg.continuity_par_values.append(float(vals[3]))
                
        if(vals[0] == "ALGO") :
            cg.algo = vals[1]
            checklist[4] = 1
            
        if(vals[0] == "MAXITER") :
            cg.miter = int(vals[1])
            checklist[5] = 1
            
        if(vals[0] == "NEVA") :
            cg.neva = int(vals[1])
            checklist[6] = 1
            
        if(vals[0] == "LBFGSB_EPS") :
            cg.LBFGSB_eps = float(vals[1])
            checklist[7] = 1            
            
        if(vals[0] == "LBFGSB_IPRINT") :
            cg.LBFGSB_iprint = int(vals[1])
            checklist[8] = 1
            
        #symm stacking default = True
        if(vals[0] == "SYMM_STCK") :
            if int(vals[1]) == 0 :
                cg.symm_stck = True
            else :
                cg.symm_stck = False
                checklist[9] = 1
    
    cfile.close()
    if cg.first_step_flag == False:
        tfile.close()

    #CHECK AND PRINT
    
    if checklist[0] == 1:
        
        cg.Nseq = len(cg.seq)
        
        for i in range(cg.Nseq):
                print("SEQUENCE " + str(i) + ": "+cg.seq[i])
                print("Njunctions Seq" + str(i) + ": "+str(cg.Njuns[i]))
                print("Single strand concentration: " + str(cg.Ct[i]) + " M")
                print("Salt concentration: " + str(cg.Cs[i]) + " M")
                print("Sim temperature: " + str(cg.simTs[i]) + " M")
                
                string = ""
                
                for l in range(len(cg.rew_Ts[i])) :
                                       
                    string = string + str(cg.rew_Ts[i][l]) + " "
                    
                print("Reweighted temperatures: " + string)
                
                string = ""
                
                for l in range(len(cg.weights[i])) :
                    string = string + str(cg.weights[i][l]) + " "
                    
                print("Weights: " +  string)
                    
                
    else:
        print("MANDATORY. Sequences missing from input file.")
        print("Usage:")
        print("SEQ seq")
        return False
    
            
    if checklist[1] == 1:
        print("INITIAL SNAPSHOT: "+str(cg.in_snap))
        print("Ignoring all sampled snapshopts < "+ str(cg.in_snap) +" in the optimisation.")
    else :
        print("OPTION. No in snapshot specified. Using all snapshots (is the trajectory equilibrated?).")
        print("Usage:")
        print("IN_SNAP in_snap")
        
    if checklist[2] == 1:
        print("NUMBER OF REPETITIONS: "+str(cg.Nreps))
    else :
        print("OPTION. No Nreps specified. Running only one simulation replica.")
        print("Usage:")
        print("IN_SNAP in_snap")
        
        
    if checklist[3] == 1:
        
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
        
        
    if checklist[4] == 1:
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
        
        
    if checklist[5] == 1:
        print("MAX number of ITERATIONS of the optimisation algorithm: "+str(cg.miter))
    else :
        print("OPTION. No max iter specified. Using default value miter = "+ str(cg.miter))
        print("Note: it might be that MAX number of function evaluations is used, instead")
        print("Usage:")
        print("MAXITER miter")

        
    if checklist[6] == 1:
        print("MAX number of function evaluations of the optimisation algorithm: "+str(cg.neva))
    else :
        print("OPTION. No function evaluations specified. Using default value miter = "+ str(cg.neva))
        print("Note: it might be that MAX number of iterations is used, instead")
        print("Usage:")
        print("NEVA neva")

    if cg.algo == "L-BFGS-B" :
        if checklist[7] == 1:
            print("Eps for L-BFGS-B algorithm: "+str(cg.LBFGSB_eps))
        else :
            print("OPTION. No eps specified for L-BFGS-B algorithm. Using default value eps = "+ str(cg.LBFGSB_eps))
            print("Note: it might be that MAX number of iterations is used, instead")
            print("Usage:")
            print("LBFGSB_EPS eps")
            
        if checklist[8] == 1:
            print("iprint for L-BFGS-B algorithm: "+str(cg.LBFGSB_iprint))
        else :
            print("OPTION. No iprint specified for L-BFGS-B algorithm. Using default value iprint = "+ str(cg.LBFGSB_iprint))
            print("Note: it might be that MAX number of iterations is used, instead")
            print("Usage:")
            print("LBFGSB_iprint iprint")
            
    if checklist[9] == 1:
        print("Breaking stacking symmetry (AA/TT only)")
    else :
        print("OPTION. Using symmetric stacking.")
        print("If you want to break the AA/TT symmetry,")
        print("Usage:")
        print("SYMM_STCK 1")            
            
    
    return True



   #change temperature to reweighting target in input2.an

def update_T_input_file(T,file_name) :
    
    strold1 = r"(T) = .*(K)$"  #could be K or 
    strold2 = r"(T) = .*(C)$"  #could be C
    strnew = "T = "+str(T)+"C"    
    
    pysed.replace(strold1,strnew,file_name)  
    pysed.replace(strold2,strnew,file_name)   
    
    return

def print_final_melting_temperatures(mTs) :
    
    ofile = open("final_rew_mTs.txt",'w')
    for i in range(len(mTs)) :
        if i%cg.Nreps == 0:
            print(mTs[i], file=ofile)
    ofile.close()
    
    return

def estimate_new_wfile(hist,wfile_name) :
    

    max_entry = max(hist)
    w = []
    for i in hist :
        if i == 0:
            print("WARNING! There is a zero in the w file. Can't produce final w file")
            w = cg.weights[cg.seq_id]
            break
        else:
            w.append(max_entry/i)
        
    ofile = open(wfile_name,'w')
    
    for i in range(len(w)):
        print(str(i) + " " + str(w[i]),file=ofile)
        
    ofile.close()
    
    return True



def check_update() :
    if cg.Diter_Trange % cg.nevery_uprewTs == 0:
        cg.Diter_Trange = 0
        return True
    else:
        return False


#update rewighting temperature range every cg.nevery_uprewTs iteration steps
def update_rew_Ts(l,update_flag):
    
    if update_flag:       
        
        print("Updating rew Ts")
        
        T0 = cg.current_mT
        
        NT = cg.n_Ts[l]
        lr = 0
        ur = 0
        
        if (NT)%2 == 0 :
            lr = -int((NT)/2-1)
            ur = int((NT)/2+1)
        else :
            lr = -int((NT-1)/2)
            ur = int((NT-1)/2+1)
            
        counts = 0
        for k in range(lr,ur) :
            
            offset = 0.
            if cg.DTs[l]*lr+T0<0 :  #oxpy does not like T < 0C and T > 100C
                offset = -(cg.DTs[l]*lr+T0-0.5)
            if cg.DTs[l]*ur+T0>100 :
                offset = - (cg.DTs[l]*ur+T0-100) - 0.5
            
            cg.rew_Ts[l][counts] = cg.DTs[l]*k+T0 + offset
            counts += 1
        
    return


def callbackF(par) :
    
    cg.Niter += 1
    cg.Diter_Trange += 1
    
    cg.update_rews = True
    
    ofile = open("parameters_v_iter.txt", 'a')
    
    string = str(cg.Niter)
    
    for i in range(len(par)):
        
        string = string + " " + str(par[i])
    
    print(string, file=ofile)
    
    ofile.close()
    
    return
   
    
def Cost_function_mT(par,stop,par0):
    
    stop[0]=cg.comm.bcast(stop[0], root=0)  #this is used to stop all processes (the while loop in main cycle)
         
    par_fr = copy.deepcopy(par)
    for k in range(len(par)):
        par[k] *= par0[k]

    #at the end of optimisation stop[0] is set to 1 and communicated to all processes   
    C = 0
    
    if stop[0] == 0 :
        
        update_flag = False
        
        if cg.rank == 0:    #check if we have to update rew Ts. This depends on the number of iteration of the opti algorithm, so rank 0
            if cg.Diter_Trange > 0 :
                update_flag = check_update()
            
        
        l = cg.seq_id
        n = cg.rep_id        
        
        update_flag=cg.comm.bcast(update_flag, root=0) #bradcast update flag
                
        update_rew_Ts(l,update_flag)
    
        #update seq dep files. Has to be done just once
        if cg.rank == 0:
            print("par0")
            print(par0)
            print("par")
            print(par)
            print("fraction")            
            print(par_fr)
                
            if cg.ave:
                functions_multi.update_rew_seq_dep_file_ave(par)
            else:
                functions_multi.update_rew_seq_dep_file(par)

                
        #bcast or Barrier make sure everything is synchronised
        par=cg.comm.bcast(par, root=0)
        #cg.comm.Barrier()
            
        hist_Ts = np.zeros([cg.n_Ts[l],len(cg.weights[l])], dtype=float)    
                   
       
        for m in range(cg.n_Ts[l]) :
            
            energy_ratio = 300.0/(cg.rew_Ts[l][m]+273.15) #300K/rewT in K 
            #print("Energy_ratio: " + str(energy_ratio))
        
            energy1 = []
            
            inp = []
            backend = []
            obs = []            
            read = False
            
                
            file_name ="./Seq"+str(l)+"/Rep"+str(n)+"/input2_melting.an"
            
            #print(file_name)
            
            update_T_input_file(cg.rew_Ts[l][m],file_name)    #change temperature to reweighting target
            print("updating T to "+ str(cg.rew_Ts[l][m]))
            
            with oxpy.Context():              
                
                #read input script specifying sequence dependent file
                inp = oxpy.InputFile()
                
                inp.init_from_filename("./Seq"+str(l)+"/Rep"+str(n)+"/input2_melting.an")
                #create analysys backend
                backend = oxpy.analysis.AnalysisBackend(inp)
            
                obs = backend.config_info().observables
                
                counts = -1
                
                while 1==1 : 
                    try:
                        read =  backend.read_next_configuration()
                    except:
                        counts+=1
                        energy1.append(999)
                        print("Warning: exception in oxpy energy computation; reweighting. Seq "+str(l)+", Rep "+str(n)+", conf" +str(counts))
                        continue
                    if read == False :
                        break
                    counts+=1
                    
                    if(counts < cg.in_snap) :
                        continue
                    a = float(obs[0].get_output_string(backend.conf_step).split()[0])
                    #stk = float(obs[nrep][1].get_output_string(backend[nrep].conf_step).split()[2])   #total stacking energy per nucleotide
                    energy1.append((cg.Njuns[l]+1)*20*a*energy_ratio)
                    #energy1_stk.append((cg.Njuns[l]+1)*20*stk)
                
            
            #reweight histogram to new parameters and temperature
            for i in range(len(energy1)) :
        
                if energy1[i] != 999 and cg.energy_sampled[i] != 999:
                    
                    deltaH = (energy1[i] - cg.energy_sampled[i])
                    
                    #if m == 5:
                    #print(cg.rank,energy1[i],cg.energy_sampled[i],deltaH)
       
                    hist_Ts[m][cg.hbs_sampled[i]] += math.exp(-deltaH)/cg.weights[l][cg.hbs_sampled[i]]   #XXXtocheck do we have to compute hbs again?
                    #hist_Ts[m][cg.hbs_sampled[l][i]] += 1./cg.weights[l][cg.hbs_sampled[l][i]]   #XXXtocheck do we have to compute hbs again?
            
        #print(hist_Ts)
    
        
        hist_Ts = cg.comm_seq.reduce(hist_Ts,op=MPI.SUM, root=0)

        
        #if cg.rank in cg.leaders: #these are the leaders
        if cg.rank_seq == 0: #same as cg.rank in cg.leaders
        
            #print(hist_Ts)
        
            #sum histograms of different repetitions of a given sequence
            print("summing on rank " + str(cg.rank))
        
            mT, mT_w = get_Tm_and_width(cg.rew_Ts[l],hist_Ts)
        
            cg.current_hist_Ts = hist_Ts
            
            cg.current_mT = mT
            
            cg.current_mT=cg.comm_seq.bcast(cg.current_mT, root=0) 
            
            if mT <= cg.rew_Ts[l][cg.n_Ts[l]-1]+0.001 and mT >= cg.rew_Ts[l][cg.n_Ts[l]-1]-0.001 :
            
                loc_C = 100000
                
            elif  mT <= cg.rew_Ts[l][0]+0.001 and mT >= cg.rew_Ts[l][0]-0.001 :
            
                loc_C = 100000
                
            else :
        
                loc_C = (mT-cg.target_mTs[l])*(mT-cg.target_mTs[l])
            
            print("seq, rep:" + str(cg.seq_id) + ", " + str(cg.rep_id) + ". RewT: " + str(mT) + ". width: " + str(mT_w) + ". target: " + str(cg.target_mTs[l])+ ". loc C: " + str(loc_C))
        
            print("rank_leaders: " + str(cg.rank_leaders))

            C = cg.comm_leaders.reduce(loc_C,op=MPI.SUM, root=0)    #summ cost function terms (only leaders)
            
        else :
            cg.current_mT=cg.comm_seq.bcast(cg.current_mT, root=0)            
            
        if cg.rank == 0 :
            
            print("tot C = " +str(C))
    
    return C

