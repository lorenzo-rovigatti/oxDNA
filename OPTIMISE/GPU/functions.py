import numpy as np
import math
import parameters_list as parl
import get_cgdna_pars
import config as cg
import cost_function as cfun

PARS_LIST = parl.PARS_LIST
par_index = parl.par_index

bases = ['A','C','G','T']
#map bases to id
def base_to_id(b) :
    if b == 'A':
        return 0
    elif  b == 'C':
        return 1
    elif  b == 'G':
        return 2
    elif  b == 'T':
        return 3
    else:
        return -1


#MODULATIONS

def Vmod(theta, a, theta0) :
    f = 1-a*(theta-theta0)**2
    return f

def Vsmooth(x,b,xc):
    f = b*(xc-x)**2
    return f

def Morse(r,epsilon,r0,a) :
    f = epsilon*(1-math.exp(-(r-r0)*a))**2
    return f

def Vharmonic(x,x0):
    f = 0.5*(x0-x)**2
    return f



 ###################################################################################################
 ############## READ OPTIMISATION PARAMETERS########################################################
 ###################################################################################################


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
            cfun.internal_coords.append([])
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
            cfun.add_opti_par(vals[1])
            checklist[8] = 1
            
        #options
        #NOTE: MODE ave NOT IMPLEMENTED
        if(vals[0] == "MODE") :
            if(vals[1]=="ave") :
                cg.ave = True
                checklist[9] = 1
            elif(vals[1]=="sd") :
                cg.ave = False
                checklist[9] = 1
        #NOTE: ALGO NOT IMPLEMENTED YET
        if(vals[0] == "ALGO") :
            cg.algo = vals[1]
            checklist[10] = 1
            
        if(vals[0] == "MAXITER") :
            cg.miter = int(vals[1])
            checklist[11] = 1
            
        if(vals[0] == "NEVA") :
            cg.neva = int(vals[1])
            checklist[12] = 1
        #NOTE: THESE OPTIONS WERE FOR THE OLDER VERSION OF THE CODE (WITHOUT TORCH)
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
                
        #optimise persistence lengths 0==true
        if(vals[0] == "T") :
            cg.T = float(vals[1])               
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
    
    print("ALL IDS:")
    print(cg.ids)
        
               
    #generate gs(mu) and covariance. Target is read from cgna+ 
    if checklist[3] == 1:
        for i in range(cg.Nseq) :
            cg.dimension.append((cg.fin_j[i]-cg.in_j[i]+1)*(len(cg.ids)))
        
            print("DIMENSION Seq "+str(i)+": " + str(cg.dimension[i]))
    
            
            if cg.ave == True and cg.diag == True:
                
                tm, tcov = get_cgdna_pars.get_target_mean_and_covariance_diag_ave((cg.fin_j[i]-cg.in_j[i]+1), cg.ids)
                
                cfun.target_mu.append(tm)
                cfun.target_cov.append(tcov)
                
            elif cg.ave == False and cg.diag == False:
                
                tm, tcov = get_cgdna_pars.get_target_mean_and_covariance(cg.seq[i], cg.ids, cg.inj, cg.jfe)
                
                cfun.target_mu.append(tm)
                cfun.target_cov.append(tcov)
                
            elif cg.ave == False and cg.diag == True:
                
                tm, tcov = get_cgdna_pars.get_target_mean_and_covariance(cg.seq[i], cg.ids, cg.inj, cg.jfe)
                
                cfun.target_mu.append(tm)
                cfun.target_cov.append(tcov)
                
                for k in range(len(cfun.target_cov[i])) :
                    for l in range(len(cfun.target_cov[i])) :
                        if k != l :
                            cfun.target_cov[i][k,l] = 0.
                        
                        
        print("TARGET GS: see file target_gs.txt")
        
        for l in range(cg.Nseq) :
            ofile = open("target_gs_Seq"+str(l)+".txt", 'w')
            for i in range(len(cg.target_mu[l])) :
                print(str(cfun.target_mu[l][i]), file = ofile)
            ofile.close()
            
        print("TARGET COV: : see file target_cov.txt")
        for l in range(cg.Nseq) :
            ofile = open("target_cov_Seq"+str(l)+".txt",'w')
            for i in range(len(cfun.target_cov[l])) :
                string = ""
                for j in range(len(cfun.target_cov[l])) :            
                    string += str(cfun.target_cov[l][i,j]) + " "
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
        
        print("PARAMETERS used in the optimisation: ")
        print("OPTIMISE - "+str(cg.par_dimension)+":")
       
        
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
        print("Temperature = "+str(cg.T))
    else:
        print("No temperature specified. Using default T = 0.1 (300K). Be sure you run the simulations at 300K.")
    
    return True


 ###################################################################################################
 ############## READ PARAMETERS FROM model.h E SD PARAMETERS FILE ##################################
 ###################################################################################################

#find initial values of the parameters from model.h file (if it gets updated, changes are read without modifying the code)
def read_vanilla_parameters(mfile) :
        
    pars_from_modelh = []
    vals_from_modelh = []

    for line in mfile.readlines() :
        vals = line.strip().split()
        for i in range(len(PARS_LIST)):
            if len(vals) >= 3:
                if PARS_LIST[i] == vals[1] or (PARS_LIST[i]+"_OXDNA2" == vals[1] and PARS_LIST[i] != "HYDR_EPS" ) : #second condition is for FENE_R0
                    pars_from_modelh.append(PARS_LIST[i])
                    if vals[2]=="PI":
                        vals_from_modelh.append(math.pi)
                    elif vals[2]=="(PI*0.5f)":
                        vals_from_modelh.append(math.pi*0.5)
                    elif len(vals) == 5 :
                        #print(vals[2]+" " +vals[3]+" "+vals[4])
                        if vals[2]+" " +vals[3]+" "+vals[4]=="(PI - 2.35f)":
                            vals_from_modelh.append(math.pi-2.35)    
                        if vals[2]+" " +vals[3]+" "+vals[4]=="(PI - 0.875f)":
                            vals_from_modelh.append(math.pi-0.875)  
                    else:
                        vals_from_modelh.append(float(vals[2][:-1]))
                    break
                
    return pars_from_modelh, vals_from_modelh
                

#read parameters from SD file
def read_pars_from_SD_file(SDfile) :
    
    #over_indices: 
    #0 - parameter index
    #1 - tetramer type (convert base 4 to base 10)
    #over_vals: corresponding parameter value
    
    over_indices = []
    over_vals = []
    stck_fact_eps = 0.18
    stck_fact_eps_read = False

    for line in SDfile.readlines() :
        vals = line.strip().split()
        if len(vals) == 0:
            continue
        if vals[0] == "STCK_FACT_EPS":
            stck_fact_eps = float(vals[2])
            stck_fact_eps_read = True
            continue
        vals_cn = vals[0].split("_")
        if vals_cn[0] == "STCK" and len(vals_cn) == 3:
            for i in range(4):
                for j in range(4):
                    ty=i+base_to_id(vals_cn[1])*4+base_to_id(vals_cn[2])*4*4+j*4*4*4    #from base 4 to base 10
                    oi = [PARS_LIST.index("STCK_EPS"),ty]
                    over_indices.append(oi)
                    over_vals.append(float(vals[2]))
        elif vals_cn[0] == "HYDR" and len(vals_cn) == 3:
            for i in range(4):
                for j in range(4):
                    ty=i+base_to_id(vals_cn[1])*4+base_to_id(vals_cn[2])*4*4+j*4*4*4    #from base 4 to base 10
                    oi = [PARS_LIST.index("HYDR_EPS"),ty]
                    over_indices.append(oi)
                    over_vals.append(float(vals[2]))
        elif vals_cn[0] == "HYDR" or vals_cn[0] == "CRST" :
            par_name = vals_cn[0]
            for i in range(1,len(vals_cn)-2):
                par_name += "_"+vals_cn[i]
            for i in range(4):
                for j in range(4):
                    ty=i+base_to_id(vals_cn[len(vals_cn)-2])*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4+j*4*4*4    #from base 4 to base 10
                    oi = [PARS_LIST.index(par_name),ty]
                    over_indices.append(oi)
                    over_vals.append(float(vals[2]))
            if vals_cn[1] == "THETA2" or vals_cn[1] == "THETA5" or vals_cn[1] == "THETA7":
                par_name = vals_cn[0]
                if vals_cn[1] == "THETA2" : 
                     par_name += "_THETA3"
                elif vals_cn[1] == "THETA5" : 
                     par_name += "_THETA6"
                elif vals_cn[1] == "THETA7" : 
                     par_name += "_THETA8"
                for i in range(2,len(vals_cn)-2):
                     par_name += "_"+vals_cn[i]
                for i in range(4):
                    for j in range(4):
                        ty=i+base_to_id(vals_cn[len(vals_cn)-2])*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4+j*4*4*4    #from base 4 to base 10
                        oi = [PARS_LIST.index(par_name),ty]
                        over_indices.append(oi)
                        over_vals.append(float(vals[2]))

        elif vals_cn[0] == "STCK" or vals_cn[0] == "FENE" :
            par_name = vals_cn[0]
            for i in range(1,len(vals_cn)-4):
                par_name += "_"+vals_cn[i]
            ty=base_to_id(vals_cn[len(vals_cn)-4])+base_to_id(vals_cn[len(vals_cn)-3])*4+base_to_id(vals_cn[len(vals_cn)-2])*4*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4*4    #from base 4 to base 10
            oi = [PARS_LIST.index(par_name),ty]
            over_indices.append(oi)
            over_vals.append(float(vals[2]))
            if vals_cn[1] == "THETA5" :
                par_name = vals_cn[0]+"_THETA6"
                for i in range(2,len(vals_cn)-4):
                     par_name += "_"+vals_cn[i]
                     oi = [PARS_LIST.index(par_name),ty]
                     over_indices.append(oi)
                     over_vals.append(float(vals[2]))
                     
                     
    return over_indices, over_vals, stck_fact_eps_read, stck_fact_eps


#over_indices and over_vals are indices and values of parameters to overwrite
#note: overwrite pars must include STCK_x_y
#T = temperature in oxdna units. This is needed to correctly set STCK_EPS
def init_oxpars(pars_mh, vals_mh, over_indices, over_vals,T,stck_fact_eps) :
    
    OXPS_zero = np.zeros((len(PARS_LIST),256),dtype=float)
    shifts = np.zeros((2,256),dtype=float) #0 = hydr, 1 = stck
    
    for i in range(len(PARS_LIST)) :
        
        if PARS_LIST[i] == "STCK_EPS":
            for j in range(256) :
                OXPS_zero[i][j] = 1.
            continue
        if PARS_LIST[i] == "HYDR_EPS":
            for j in range(256) :
                OXPS_zero[i][j] = 0.
            continue

        index = pars_mh.index(PARS_LIST[i])
        val = vals_mh[index]
        
        for j in range(256) :
            OXPS_zero[i][j] = val
    
    #here we use the initial custom parameters, must include stck_x_y and hydr_x_y
    if [4,48] not in over_indices :
        print("No HYDR_x_y in SD file. Terminating")
        exit(1)
    if [44,0] not in over_indices :
        print("No STCK_x_y in SD file. Terminating")
        exit(1)
    for i in range(len(over_indices)) :
        OXPS_zero[over_indices[i][0]][over_indices[i][1]] = over_vals[i]
        
    #set eps and shifts
    for j in range(256) :
        #hydr
        shifts[0][j] = Morse(OXPS_zero[par_index[7]][j],OXPS_zero[par_index[4]][j],OXPS_zero[par_index[5]][j],OXPS_zero[par_index[6]][j])
        #stacking
        OXPS_zero[par_index[44]][j] =  OXPS_zero[par_index[44]][j]* (1.0 - stck_fact_eps + (T * 9.0 * stck_fact_eps))
        shifts[1][j] = Morse(OXPS_zero[par_index[47]][j],OXPS_zero[par_index[44]][j],OXPS_zero[par_index[45]][j],OXPS_zero[par_index[46]][j])
    
    return OXPS_zero, shifts

###################################################################################################
############## READ TRAJECTORY AND COMPUTE OXDNA COORDINATES (i.e angles and distances) ###########
###################################################################################################

#compute rcut
#hydr distance cutoff. If rhydr is larger than this, then both hydr and crst are zero.
def find_cuts_for_lists(OXPS_zero) :

    rcut_high = 0.
    rcut_low = 1000.
    
    
    for j in range(len(OXPS_zero[par_index[13]])):
        if OXPS_zero[par_index[13]][j] > rcut_high:
            rcut_high = OXPS_zero[par_index[13]][j]
        if OXPS_zero[par_index[85]][j] > rcut_high:
            rcut_high = OXPS_zero[par_index[85]][j]
        if OXPS_zero[par_index[124]][j] > rcut_high:
            rcut_high = OXPS_zero[par_index[124]][j]
        if OXPS_zero[par_index[12]][j] < rcut_low:
            rcut_low = OXPS_zero[par_index[12]][j]
        if OXPS_zero[par_index[84]][j] < rcut_low:
            rcut_low = OXPS_zero[par_index[84]][j]
        if OXPS_zero[par_index[123]][j] < rcut_low:
            rcut_low = OXPS_zero[par_index[123]][j]
    
    rcut_high = rcut_high + 0.000005
    rcut_low = rcut_low - 0.000005
    
    return rcut_low , rcut_high


class topo :
    def __init__(self, nid, sid, bty, do, up) :
        self.id = nid
        self.strand_id = sid
        self.base_type = bty
        self.down_id = do
        self.up_id = up


pos_bb1 = [-0.34,-0.34,-0.34,-0.34]
pos_bb2 = [0.3408,0.3408,0.3408,0.3408]
#note: oxdna2
pos_stck = [0.34,0.34,0.34,0.34]
pos_hydr = [0.4,0.4,0.4,0.4]
#note: oxdna3
#pos_stck = [0.37,0.37,0.37,0.37]
#pos_hydr = [0.43,0.37,0.43,0.37]
class nucleotide :
    def __init__(self, C, BV, N, ty) :
        t = base_to_id(ty)
        self.type = t
        self.c = C
        self.bv = BV
        self.n = N
        self.norv = np.cross(N,BV)
        self.bb = C + pos_bb1[t]*BV + pos_bb2[t]*self.norv
        self.stck = C + pos_stck[t]*BV
        self.hydr = C + pos_hydr[t]*BV
    

def read_oxdna_trajectory_dist_and_angles(rcut_low, rcut_high, tr_file, topo_file):
    
    rcut_sq_high = rcut_high*rcut_high
    rcut_sq_low = rcut_low*rcut_low
    
    
    #define tensors

    #bonded
    fene_r = []
    stck_r = []
    th4_bn = []
    th5 = []
    th6 = []
    cosphi1 = []
    cosphi2 = []
    types_bn = []

    #unbonded

    hydr_r = []
    th1 = []
    th2 = []
    th3 = []
    th4_unbn = []
    th7 = []
    th8 = []
    types_unbn = []
    
     
    Nb = 0
    Ns = 0
    nid = 0
    
    topology = []
    counts = 0
    for line in topo_file.readlines() :
        counts+=1
        vals = line.split()
        if len(vals) == 2 :
            Nb = int(vals[0])
            Ns = int(vals[1])
            if Ns != 2 :
                print("Number of strands in not 2.")
                exit(1)
        else :
           to = topo(nid, int(vals[0]), vals[1], int(vals[2]), int(vals[3]))
           topology.append(to)
           nid += 1
           
    counts = 0
    for line in tr_file.readlines():
        a = line.strip()[0]
        if a == 't':
            nid = 0
            config = []
            counts+=1
            #print(counts)
            #if counts == 501:
            #   break
            counts_b = -1
        else:
            if a != 'b' and a != 'E':
                counts_b += 1
                vals = line.split()
                c = np.array([float(vals[0]),float(vals[1]), float(vals[2])])
                bv = np.array([float(vals[3]),float(vals[4]), float(vals[5])])
                n = np.array([float(vals[6]),float(vals[7]), float(vals[8])])
                config.append(nucleotide(c,bv,n, topology[counts_b].base_type))
                nid += 1
                #configuration read. Computing coordinates
                if len(config) == Nb :
                    
                    #compute coordinates
                    fene_r_1conf = []
                    stck_r_1conf = []
                    th4_bn_1conf = []
                    th5_1conf = []
                    th6_1conf = []
                    cosphi1_1conf = []
                    cosphi2_1conf = []
                    types_1conf = []

                    #bonded basepairs
                    #find pairs and tetramer type
                    for i in range(len(topology)) :
                        ty0 = 0
                        ty1 = 0
                        ty2 = 0
                        ty3 = 0
                        if topology[i].up_id == -1: continue
                        n1 = config[i]
                        ty1 = base_to_id(topology[i].base_type)
                        n2 = n1

                        if topology[i].down_id != -1:
                         for j in range(len(topology)) :
                             if topology[j].id == topology[i].down_id:
                                ty0 = base_to_id(topology[j].base_type)
                                break
                        for j in range(len(topology)) :
                            if topology[j].id == topology[i].up_id:
                               ty2 = base_to_id(topology[j].base_type)
                               n2 = config[j]
                               if topology[j].up_id == -1:
                                   ty3 = ty0
                               else:
                                   for z in range(len(topology)) :
                                       if topology[z].id == topology[j].up_id:
                                          ty3 = base_to_id(topology[z].base_type)
                                          break
                                      
                               break
                        if topology[i].down_id == -1:
                            ty0 = ty3
                            
                        ty = ty0+ty1*4+ty2*4*4+ty3*4*4*4 #tetramer type in base 10
                        
                        types_1conf.append(ty)
                        
                        #compute bnd pair coordinates
                        
                        fene_r_1conf.append( np.linalg.norm(n1.bb-n2.bb) )
                        stck_r_1conf.append( np.linalg.norm(n1.stck-n2.stck) )
                        th4_bn_1conf.append( np.arccos(np.dot(n1.n,n2.n)) )
                        

            			#note rbb is the distance between backbone sites  in oxdna1 (see standalone)!
                        bp1 = n1.c - 0.4*n1.bv
                        bp2 = n2.c - 0.4*n2.bv
                        rbb = (bp1 - bp2)/np.linalg.norm((bp1 - bp2))
                        rstck = (n1.stck - n2.stck)/np.linalg.norm((n1.stck - n2.stck))
                        
                        th5_1conf.append( np.arccos(-np.dot(n2.n,rstck)))
                        th6_1conf.append( np.arccos(-np.dot(n1.n,rstck)))
                        cosphi1_1conf.append( np.dot(n2.norv,rbb))
                        cosphi2_1conf.append( np.dot(n1.norv,rbb))


                    types_bn.append(types_1conf)
                    fene_r.append(fene_r_1conf)
                    stck_r.append(stck_r_1conf)
                    th4_bn.append(th4_bn_1conf)
                    th5.append(th5_1conf)
                    th6.append(th6_1conf)
                    cosphi1.append(cosphi1_1conf)
                    cosphi2.append(cosphi2_1conf)
            
                    hydr_r_1conf = []
                    th1_1conf = []
                    th2_1conf = []
                    th3_1conf = []
                    th4_unbn_1conf = []
                    th7_1conf = []
                    th8_1conf = []
                    types_unbn_1conf = []
            
                    #TODO UNDBONDED
                    for i in range(len(topology)) :
                        ty0 = 0
                        ty1 = 0
                        ty2 = 0
                        ty3 = 0
                        n1 = config[i]
                        ty1 = base_to_id(topology[i].base_type)
                        for z in range(len(topology)) :
                           if topology[z].id == topology[i].down_id:
                              ty0 = base_to_id(topology[z].base_type)
                              break
                        ty0 = base_to_id(topology[i].base_type)
                        for j in range(len(topology)) :
                            n2 = config[j]
                            if np.dot(n1.hydr - n2.hydr, n1.hydr - n2.hydr) < rcut_sq_low: continue #verlet cutoff
                            if np.dot(n1.hydr - n2.hydr, n1.hydr - n2.hydr) > rcut_sq_high: continue #verlet cutoff
                            if topology[j].id <= topology[i].id: continue #ordered ids (avoid pairs repetition)
                            if topology[j].id == topology[i].down_id or topology[j].id == topology[i].up_id: continue #no bonded pairs
                            
                            ty2 = base_to_id(topology[j].base_type)
                            
                            if topology[j].up_id == -1:
                               ty3 = ty0
                            else:
                                for z in range(len(topology)) :
                                   if topology[z].id == topology[j].up_id:
                                      ty3 = base_to_id(topology[z].base_type)
                                      break
                        
                            if topology[i].down_id == -1:
                               ty0 = ty3
                               
                            ty = ty0+ty1*4+ty2*4*4+ty3*4*4*4 #tetramer type in base 10
                        
                            types_unbn_1conf.append(ty)
                        
                            #compute unbnd pair coordinates
                            hydr_r_1conf.append(np.linalg.norm((n1.hydr - n2.hydr)))
                            rhydr = (n1.hydr - n2.hydr)/np.linalg.norm((n1.hydr - n2.hydr))

                            th1_1conf.append(np.acos(-np.dot(n1.bv,n2.bv)))
                            th3_1conf.append(np.acos(-np.dot(n1.bv,rhydr)))
                            th2_1conf.append(np.acos(np.dot(n2.bv,rhydr)))
                        
                            th4_unbn_1conf.append(np.acos(np.dot(n1.n,n2.n)))
                            th8_1conf.append(np.acos(-np.dot(n1.n,rhydr)))
                            th7_1conf.append(np.acos(np.dot(n2.n,rhydr)))
                        

                    types_unbn.append(types_unbn_1conf)
                    hydr_r.append(hydr_r_1conf)

                    th1.append(th1_1conf)
                    th2.append(th2_1conf)
                    th3.append(th3_1conf)
                    
                    th4_unbn.append(th4_unbn_1conf)
                    th7.append(th7_1conf)
                    th8.append(th8_1conf)
                    
                    
    #make unbnd tensor square. Extra unbnd pairs have zero interaction energy.

    max_ints = 0
    for j in range(len(types_unbn)):
        if len(types_unbn[j]) > max_ints:
           max_ints = len(types_unbn[j])
    print("max unbn pairs: "+str(max_ints))
    for j in range(len(types_unbn)):
        for z in range(len(types_unbn[j]), max_ints):
            types_unbn[j].append(0)
            hydr_r[j].append(0.)

            th1[j].append(0.)
            th2[j].append(0.)
            th3[j].append(0.)

            th4_unbn[j].append(0.)
            th7[j].append(0.)
            th8[j].append(0.)
                             
                    
                    
        #bonded
        fene_r = []
        stck_r = []
        th4_bn = []
        th5 = []
        th6 = []
        cosphi1 = []
        cosphi2 = []
        types_bn = []

        #unbonded

        hydr_r = []
        th1 = []
        th2 = []
        th3 = []
        th4_unbn = []
        th7 = []
        th8 = []
        types_unbn = []                
                    
    return fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3, th4_unbn, th7, th8, types_unbn


###################################################################################################
############## STORE INTERNAL COORDINATES (i.e angles and distances) ##############################
###################################################################################################

#given a junction trajectory (read_oxdna_trajectory_standard_order), store specific internal coordinates in
#global variable internal_coords
def store_internal_coord(traj,seq,ids,in_j,fin_j,in_snap,overwrite=True) :
       
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
            cfun.internal_coords[seq].append(coord)
        else :
            coords.append(coord)
           
    if overwrite == True :
        cfun.internal_coords[seq] = coords
                    
    return

def ave_and_cov_sampled() :
      
    mu_sampled = []
    cov_sampled = []
    
    for l in range(cg.Nseq) :
        Nsnaps = len(cfun.internal_coords[l])
        Ncoords = len(cfun.internal_coords[l][0])
        co = []
        ma = []
        for m in range(Ncoords) :
            co.append(0.)
        for m in range(Ncoords) :
            ma.append(co)
            
        mu_sampled.append(co)
        cov_sampled.append(ma)
    
    for l in range(cg.Nseq) :
        Nsnaps = len(cfun.internal_coords[l])
        Ncoords = len(cfun.internal_coords[l][0])
        
        for i in range(Ncoords) :
            mu_sampled[l][i] = 0.
            for j in range(Ncoords) :
                cov_sampled[l][i][j] = 0.
        
        for i in range(Nsnaps) :
            for j in range(Ncoords) :
                mu_sampled[l][j] += cfun.internal_coords[l][i][j]/Nsnaps
        
        for i in range(Nsnaps) :
            for j in range(Ncoords) :
                for z in range(j,Ncoords) :
                    cov_sampled[l][j][z] += (cfun.internal_coords[l][i][j] - mu_sampled[l][j])*(cfun.internal_coords[l][i][z] - mu_sampled[l][z])/Nsnaps
        
        for j in range(Ncoords) :
            for z in range(j+1,Ncoords) :
                cov_sampled[l][z][j] = cov_sampled[l][j][z]    
    
    return mu_sampled, cov_sampled
