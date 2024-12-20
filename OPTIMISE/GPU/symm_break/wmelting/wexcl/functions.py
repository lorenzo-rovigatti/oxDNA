import numpy as np
import math
import parameters_list as parl
import get_cgdna_pars
import config as cg
import cost_function as cfun
import matplotlib.pyplot as plt
import torch


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
############## READ OPTIMISATION CONFIG    ########################################################
###################################################################################################


#Parse config file (read parameters)
#Returns False if mandatory parameters are missing
def read_config(cfile_name) :

    cfile = open(cfile_name,'r')

    checklist = np.zeros(19, dtype = int) #check if some mandatory parameters are missing

    for line in cfile.readlines() :
        vals = line.strip().split()
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

        #Switch on and off lp optimisation.
        #Default = False
        if(vals[0] == 'LP') :
            if vals[1] == "0" or vals[1] == "True" or vals[1] == "true" or vals[1] == "TRUE" :
                cg.opti_lp = True
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
        if(vals[0] == "TEMPERATURE") :
            if vals[1][len(vals[1])-1] == "K":
                cg.T = float(vals[1][:-1])/300.*0.1
            elif vals[1][len(vals[1])-1] == "C":
                cg.T = (float(vals[1][:-1])+273.15)/300.*0.1
            else:
                print("Can't recognise temperature unit. K=kelvin, C=celsius")
                print("e.g. T 300K or T 25C")
                print("Using T = 300K")

            checklist[17] = 1

        #optimise persistence lengths 0==true
        if(vals[0] == "MODELH") :
            cg.modelh = vals[1]

            checklist[18] = 1

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


    if cg.opti_lp == True :
        print("Optimising persistence lengths (i.e. Ar and C)")
    else:
        print("Neglecting persistence lengths (i.e. Ar and C)")
        print("Usage:")
        print("LPS i")
        print("i = 0 or i = true for optimising lps.")


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

    if cg.opti_lp == True :
        for i in range(len(cg.ids_lp)) :
            if cg.ids_lp[i] in cg.ids :
                continue
            else :
                cg.ids.append(cg.ids_lp[i])

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
                            cfun.target_cov[i][k][l] = 0.


        print("TARGET GS: see file target_gs.txt")

        for l in range(cg.Nseq) :
            ofile = open("target_gs_Seq"+str(l)+".txt", 'w')
            for i in range(len(cfun.target_mu[l])) :
                print(str(cfun.target_mu[l][i]), file = ofile)
            ofile.close()

        print("TARGET COV: : see file target_cov.txt")
        for l in range(cg.Nseq) :
            ofile = open("target_cov_Seq"+str(l)+".txt",'w')
            for i in range(len(cfun.target_cov[l])) :
                string = ""
                for j in range(len(cfun.target_cov[l])) :
                    string += str(cfun.target_cov[l][i][j]) + " "
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

        cg.par_dimension = len(cfun.OPT_PAR_LIST)

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

    if checklist[18] == 1:
        print("Using model.h file: "+cg.modelh)
    else:
        print("No model.h file specified")
        print("Using default: "+cg.modelh)

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
        if vals[0][0] == '#':
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
            if vals_cn[0] == "HYDR" or vals_cn[1] == "K":
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
            else :
                for i in range(1,len(vals_cn)-4):
                    par_name += "_"+vals_cn[i]
                ty=base_to_id(vals_cn[len(vals_cn)-4])+base_to_id(vals_cn[len(vals_cn)-3])*4+base_to_id(vals_cn[len(vals_cn)-2])*4*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4*4
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
                    for i in range(2,len(vals_cn)-4):
                         par_name += "_"+vals_cn[i]
                         ty=base_to_id(vals_cn[len(vals_cn)-4])+base_to_id(vals_cn[len(vals_cn)-3])*4+base_to_id(vals_cn[len(vals_cn)-2])*4*4+base_to_id(vals_cn[len(vals_cn)-1])*4*4*4
                         oi = [PARS_LIST.index(par_name),ty]
                         over_indices.append(oi)
                         over_vals.append(float(vals[2]))

        elif vals_cn[0] == "STCK" or vals_cn[0] == "FENE" or vals_cn[0] = "EXCL":
            if vals_cn[0] == "FENE" and vals_cn[1] == "EPS" : par_name = "FENE_EPS"
            else:
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
    #tmp_ind = [0,4,44,48]
    for i in range(len(over_indices)) :
        #if over_indices[i][0] not in tmp_ind : continue
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

    #+-0.2 because we are changing crst
    rcut_high = rcut_high + 0.2
    rcut_low = rcut_low - 0.2

    return rcut_low , rcut_high

#down id = bonded nucleotide, n3 direction
#up id = bonded nucleotide, n5 direction
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

def read_melting_seqs_and_mTs(ifile) :
    mTs = []
    seqs = []
    en_offsets = []
    en0 = []

    for line in ifile.readlines():
        vals = line.strip().split()
        seqs.append(vals[0])
        mTs.append(float(vals[1]))
        en_offsets.append(float(vals[2]))
        en0.append(float(vals[3]))

    return seqs, mTs, en_offsets, en0

def read_topology_from_file(topo_file) :

    topology = []
    counts = 0
    nid = 0
    Nb = 0
    Ns = 0

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

    return topology, Nb, Ns

def read_topology_from_list(topo_data) :

    Nb = len(topo_data)*2
    Ns = 2
    topology = []
    nid = 0
    for l in range(len(topo_data)) :
        if l < len(topo_data)-1 : to = topo(nid, 1, topo_data[l], nid-1,nid+1)
        else : to = topo(nid, 1, topo_data[l], nid-1,-1)
        nid+=1
        topology.append(to)
    for l in range(len(topo_data)) :
        if l == 0: to = topo(nid, 2, cg.bases[3-cg.bases.index(topo_data[len(topo_data)-l-1])], -1,nid+1)
        elif l < len(topo_data)-1 : to = topo(nid, 2, cg.bases[3-cg.bases.index(topo_data[len(topo_data)-l-1])], nid-1,nid+1)
        else : to = topo(nid, 2, cg.bases[3-cg.bases.index(topo_data[len(topo_data)-l-1])], nid-1,-1)
        nid+=1
        topology.append(to)

    return topology, Nb, Ns

def read_oxdna_trajectory_dist_and_angles(rcut_low, rcut_high, tr_file, topo_data, isfile=True):

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

    ba_ba_r_bn = []
    ba_bb_r_bn = []
    bb_ba_r_bn = []


    #unbonded

    hydr_r = []
    th1 = []
    th2 = []
    th3 = []
    th4_unbn = []
    th7 = []
    th8 = []
    types_unbn_33 = []
    types_unbn_55 = []

    bb_bb_r_unbn = []
    ba_bb_r_unbn = []
    bb_ba_r_unbn = []


    Nb = 0
    Ns = 0
    nid = 0
    topology = None

    if isfile : topology, Nb, Ns = read_topology_from_file(topo_data)
    else: topology, Nb, Ns = read_topology_from_list(topo_data)

    counts = 0
    counts_b = -1
    config = []
    nid = 0
    for line in tr_file.readlines():
        a = line.strip()[0]
        if counts <= cg.in_snap and isfile == True:
            if a == 't': counts += 1
            continue
        if a == 't':
            nid = 0
            config = []
            counts+=1
            counts_b = -1
        else:
            if a != 'b' and a != 'E':
                #if isfile == False: topology = topologies[counts-1]
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
                    ba_ba_bn_1conf = []
                    ba_bb_bn_1conf = []
                    bb_ba_bn_1conf = []

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
                        prod = np.dot(n1.n,n2.n)
                        if prod > 1: prod -= 1e-12  #this is for the 0th configurations
                        if prod < -1: prod += 1e-12

                        th4_bn_1conf.append( np.arccos(prod) )


            	        #note rbb is the distance between backbone sites in oxdna1, not in oxdna2/3 (see standalone)!
                        bp1 = n1.c - 0.4*n1.bv
                        bp2 = n2.c - 0.4*n2.bv
                        rbb = (bp1 - bp2)/np.linalg.norm((bp1 - bp2))
                        rstck = (n1.stck - n2.stck)/np.linalg.norm((n1.stck - n2.stck))


                        th5_1conf.append( np.arccos(-np.dot(n2.n,rstck)))
                        th6_1conf.append( np.arccos(-np.dot(n1.n,rstck)))
                        cosphi1_1conf.append( np.dot(n2.norv,rbb))
                        cosphi2_1conf.append( np.dot(n1.norv,rbb))

                        ba_ba_r_bn_1conf.append(np.linalg.norm(n1.hydr-n2.hydr))
                        ba_bb_r_bn_1conf.append(np.linalg.norm(n1.hydr-n2.bb))
                        bb_ba_r_bn_1conf.append(np.linalg.norm(n1.bb-n2.hydr))


                    types_bn.append(types_1conf)
                    fene_r.append(fene_r_1conf)
                    stck_r.append(stck_r_1conf)
                    th4_bn.append(th4_bn_1conf)
                    th5.append(th5_1conf)
                    th6.append(th6_1conf)
                    cosphi1.append(cosphi1_1conf)
                    cosphi2.append(cosphi2_1conf)
                    ba_ba_r_bn.append(ba_ba_r_bn_1conf)
                    ba_bb_r_bn.append(ba_bb_r_bn_1conf)
                    bb_ba_r_bn.append(bb_ba_r_bn_1conf)

                    hydr_r_1conf = []
                    th1_1conf = []
                    th2_1conf = []
                    th3_1conf = []
                    th4_unbn_1conf = []
                    th7_1conf = []
                    th8_1conf = []
                    types_unbn_1conf_33 = []
                    types_unbn_1conf_55 = []

                    #UNDBONDED
                    for i in range(len(topology)) :
                        ty0_33 = 0
                        ty0_55 = 0
                        ty1 = 0
                        ty2 = 0
                        ty3_33 = 0
                        ty3_55 = 0
                        n1 = config[i]
                        ty1 = base_to_id(topology[i].base_type)
                        for z in range(len(topology)) :
                           if topology[z].id == topology[i].down_id:
                              ty0_55 = base_to_id(topology[z].base_type)
                           if topology[z].id == topology[i].up_id:
                              ty0_33 = base_to_id(topology[z].base_type)

                        for j in range(len(topology)) :
                            n2 = config[j]
                            if np.dot(n1.hydr - n2.hydr, n1.hydr - n2.hydr) < rcut_sq_low: continue #verlet cutoff
                            if np.dot(n1.hydr - n2.hydr, n1.hydr - n2.hydr) > rcut_sq_high: continue #verlet cutoff
                            if topology[j].id <= topology[i].id: continue #ordered ids (avoid pairs repetition)
                            if topology[j].id == topology[i].down_id or topology[j].id == topology[i].up_id: continue #no bonded pairs

                            ty2 = base_to_id(topology[j].base_type)

                            if topology[j].up_id == -1:
                               ty3_33 = 0
                            if topology[j].down_id == -1:
                               ty3_55 = 0
                            else:
                                for z in range(len(topology)) :
                                   if topology[z].id == topology[j].down_id:
                                      ty3_55 = base_to_id(topology[z].base_type)
                                   if topology[z].id == topology[j].up_id:
                                      ty3_33 = base_to_id(topology[z].base_type)

                            if topology[i].down_id == -1:
                               ty0_55 = 0
                            if topology[i].up_id == -1:
                               ty0_33 = 0

                            ty_33 = ty0_33+ty1*4+ty2*4*4+ty3_33*4*4*4 #tetramer type in base 10
                            ty_55 = ty0_55+ty1*4+ty2*4*4+ty3_55*4*4*4 #tetramer type in base 10

                            types_unbn_1conf_33.append(ty_33)
                            types_unbn_1conf_55.append(ty_55)

                            #compute unbnd pair coordinates
                            hydr_r_1conf.append(np.linalg.norm((n1.hydr - n2.hydr)))
                            rhydr = (n1.hydr - n2.hydr)/np.linalg.norm((n1.hydr - n2.hydr))
                            prod = -np.dot(n1.bv,n2.bv)
                            if prod > 1: prod -= 1e-12  #this is for the 0th configurations
                            if prod < -1: prod += 1e-12
                            th1_1conf.append(np.arccos(prod))
                            prod = -np.dot(n1.bv,rhydr)
                            if prod > 1: prod -= 1e-12  #this is for the 0th configurations
                            if prod < -1: prod += 1e-12
                            th3_1conf.append(np.arccos(prod))
                            prod = np.dot(n2.bv,rhydr)
                            if prod > 1: prod -= 1e-12  #this is for the 0th configurations
                            if prod < -1: prod += 1e-12
                            th2_1conf.append(np.arccos(prod))

                            prod = np.dot(n1.n,n2.n)
                            if prod > 1: prod -= 1e-12  #this is for the 0th configurations
                            if prod < -1: prod += 1e-12
                            th4_unbn_1conf.append(np.arccos(prod))
                            th8_1conf.append(np.arccos(-np.dot(n1.n,rhydr)))
                            th7_1conf.append(np.arccos(np.dot(n2.n,rhydr)))

                            bb_bb_r_unbn_1conf.append(np.linalg.norm(n1.bb-n2.bb))
                            ba_bb_r_unbn_1conf.append(np.linalg.norm(n1.hydr-n2.bb))
                            bb_ba_r_unbn_1conf.append(np.linalg.norm(n1.bb-n2.hydr))


                    types_unbn_33.append(types_unbn_1conf_33)
                    types_unbn_55.append(types_unbn_1conf_55)

                    hydr_r.append(hydr_r_1conf)

                    th1.append(th1_1conf)
                    th2.append(th2_1conf)
                    th3.append(th3_1conf)

                    th4_unbn.append(th4_unbn_1conf)
                    th7.append(th7_1conf)
                    th8.append(th8_1conf)

                    bb_bb_r_unbn.append(bb_bb_r_unbn_1conf)
                    ba_bb_r_unbn.append(ba_bb_r_unbn_1conf)
                    bb_ba_r_unbn.append(bb_ba_r_unbn_1conf)


    return fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, ba_ba_r_bn, ba_bb_r_bn, bb_ba_r_bn, types_bn, \
           bb_bb_r_bn, ba_bb_r_unbn, bb_ba_r_unbn, hydr_r, th1, th2, th3, th4_unbn, th7, th8, types_unbn_33, types_unbn_55


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
                    coord.append(float(traj[i][j].base_pair1.intra_coord.rot[0]*alpha))
                elif ids[z] == 1 :
                    coord.append(float(traj[i][j].base_pair1.intra_coord.rot[1]*alpha))
                elif ids[z] == 2 :
                    coord.append(float(traj[i][j].base_pair1.intra_coord.rot[2]*alpha))
                elif ids[z] == 3 :
                    coord.append(float(traj[i][j].base_pair1.intra_coord.tran[0]))
                elif ids[z] == 4 :
                    coord.append(float(traj[i][j].base_pair1.intra_coord.tran[1]))
                elif ids[z] == 5 :
                    coord.append(float(traj[i][j].base_pair1.intra_coord.tran[2]))

                elif ids[z] == 6 :
                    coord.append(float(traj[i][j].inter_coord.rot[0]*alpha))
                elif ids[z] == 7 :
                    coord.append(float(traj[i][j].inter_coord.rot[1]*alpha))
                elif ids[z] == 8 :
                    coord.append(float(traj[i][j].inter_coord.rot[2]*alpha))
                elif ids[z] == 9 :
                    coord.append(float(traj[i][j].inter_coord.tran[0]))
                elif ids[z] == 10 :
                    coord.append(float(traj[i][j].inter_coord.tran[1]))
                elif ids[z] == 11 :
                    coord.append(float(traj[i][j].inter_coord.tran[2]))

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

        mu_sampled_one = np.zeros(Ncoords,dtype=float)
        cov_sampled_one = np.zeros((Ncoords,Ncoords),dtype=float)

        mu_sampled.append(mu_sampled_one)
        cov_sampled.append(cov_sampled_one)

        for i in range(Ncoords) :
            mu_sampled[l][i] = 0.
            for j in range(Ncoords) :
                cov_sampled[l][i][j] = 0.

        for i in range(Nsnaps) :
            for j in range(Ncoords) :
                mu_sampled[l][j] += cfun.internal_coords[l][i][j]/Nsnaps

        for j in range(Ncoords) :
            for z in range(Ncoords) :
                for i in range(Nsnaps) :
                    cov_sampled[l][j][z] += (cfun.internal_coords[l][i][j] - mu_sampled[l][j])*(cfun.internal_coords[l][i][z] - mu_sampled[l][z])/Nsnaps


    return mu_sampled, cov_sampled

#compute sequence averaged ground state
def Sequence_ave_GS(mu_sampled) :

    mu_sdave_sampled = np.zeros((cg.Nseq,len(cg.ids)),dtype=float)

    for l in range(cg.Nseq) :
        for i in range(len(mu_sampled[l])) :
            mu_sdave_sampled[l][i%len(cg.ids)] += mu_sampled[l][i]/len(mu_sampled[l])*len(cg.ids)

    return mu_sdave_sampled

#FROM INTERNAL COORDINATES, COMPUTE ANGLES FOR PERSISTENCE LENGTH, COSTHETAB and COSOMEGAT

e3s = []

#store normals of bps
def store_normals(traj,seq,in_j,fin_j,in_snap,overwrite=True) :

    global e3s

    print(seq)

    e3s_loc = []

    Nsnaps = len(traj)
    Njuns = len(traj[0])

    for i in range(in_snap,Nsnaps) :
        e3 = []
        for j in range(Njuns) :
            if j < in_j or j > fin_j :
                continue
            #print(traj[i][j].base_pair1.frame.orientation[:,2])
            e3.append(traj[i][j].base_pair1.frame.orientation[:,2].real)

        if overwrite == False :
            e3s[seq].append(e3)
        else :
            e3s_loc.append(e3)

    if overwrite == True :
        e3s[seq] = e3s_loc

    return

#compute cos theta and cos omega3 for persistence lengths
def plengths_angles() :

    m = cg.lp_m


    print(len(e3s))
    print(len(e3s[0]))
    print(len(e3s[0][0]))
    print(len(e3s[0][0][0]))

    costb = []
    for l in range(cg.Nseq): costb.append([])

    for l in range(cg.Nseq):    #loop sequences
        N = int(len(cfun.internal_coords[l][0])/len(cg.ids))
        for i in range(len(e3s[l])):	#loop sampled confs
            ctb = 0.
            counts = 0
            for j in range(N-m):
    	        #print(i,j)
    	        #print(np.dot(e3s[i][j],e3s[i][j+m]))
    	        ctb += np.dot(e3s[l][i][j],e3s[l][i][j+m])
    	        counts+=1

            costb[l].append(ctb/counts)

    #print(costb)

    cosot = []
    for l in range(cg.Nseq): cosot.append([])

    twist_idx = cg.ids.index(8)

    for l in range(cg.Nseq):    #loop sequences
        N = int(len(cfun.internal_coords[l][0])/len(cg.ids))
        ave_twist = cfun.save_mu[l][twist_idx]
        for i in range(len(e3s[l])):       #loop sampled confs
            cot = 0.
            counts = 0
            for j in range(N-m):
                om3 = 0.
                for z in range(j,j+m):
                    om3 += (cfun.internal_coords[l][i][z*len(cg.ids)+twist_idx]-ave_twist)

                cot += math.cos(om3/5.) #/5. is from cgna to radiant
                counts+=1

            cosot[l].append(cot/counts)

    #print(cosot)

    return costb, cosot

#PLOT SAMPLED GS AND STD


def unscrumble_gs(gs) :

    unsc_gs_all = []


    for j in range(len(cg.ids)) :

        unsc_gs = []

        for z in range(len(gs)) :
            if z%len(cg.ids) == j :
                unsc_gs.append(gs[z])

        unsc_gs_all.append(unsc_gs)

    return unsc_gs_all


def plot_gs_sampled(mu_sampled,seqid) :

    unscr_gs_sampled = unscrumble_gs(mu_sampled)
    unscr_gs_target = unscrumble_gs(cfun.target_mu[seqid])

    for j in range(len(unscr_gs_sampled)) : #Nids

        coord_name = ""

        if cg.ids[j] == 0:
            coord_name = "buckle"
        if cg.ids[j] == 1:
            coord_name = "propeller"
        if cg.ids[j] == 2:
            coord_name = "opening"
        if cg.ids[j] == 3:
            coord_name = "shear"
        if cg.ids[j] == 4:
            coord_name = "stretch"
        if cg.ids[j] == 5:
            coord_name = "stagger"
        if cg.ids[j] == 6:
            coord_name = "tilt"
        if cg.ids[j] == 7:
            coord_name = "roll"
        if cg.ids[j] == 8:
            coord_name = "twist"
        if cg.ids[j] == 9:
            coord_name = "shift"
        if cg.ids[j] == 10:
            coord_name = "slide"
        if cg.ids[j] == 11:
            coord_name = "rise"

        x = []
        for z in range(len(unscr_gs_sampled[j])) :
            x.append(z)

        ys = unscr_gs_sampled[j]
        yt = unscr_gs_target[j]


        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        # plt.title(r"All: $R=7$ $\phi=0.364$")
        ax.title.set_fontsize(20)
        ax.set_xlabel(r"Sequence",fontsize=20)
        ax.set_ylabel(coord_name,fontsize=20)
        #ax.set_ylim(0,160)
        #ax.set_xlim(-1.2,1.2)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        ax.plot(x, yt, 'b', label="cgna+")
        ax.plot(x, ys, 'r', label="oxdna")
        ax.legend(fontsize = 20)
        plt.savefig("Seq"+str(seqid)+"_"+coord_name+".pdf",bbox_inches='tight',pad_inches=0.05)

        plt.close()

    return


def unscrumble_cov_diag(cov) :

    unsc_cov_all = []

    for j in range(len(cg.ids)) :

        unsc_cov = []

        for z in range(len(cov)) :
            if z%len(cg.ids) == j :
                unsc_cov.append(cov[z][z])

        unsc_cov_all.append(unsc_cov)

    return unsc_cov_all


def plot_std_sampled(cov_sampled, seqid) :

    unscr_cov_sampled = unscrumble_cov_diag(cov_sampled)
    unscr_cov_target = unscrumble_cov_diag(cfun.target_cov[seqid])

    for j in range(len(unscr_cov_sampled)) : #Nids

        coord_name = ""

        if cg.ids[j] == 0:
            coord_name = "buckle"
        if cg.ids[j] == 1:
            coord_name = "propeller"
        if cg.ids[j] == 2:
            coord_name = "opening"
        if cg.ids[j] == 3:
            coord_name = "shear"
        if cg.ids[j] == 4:
            coord_name = "stretch"
        if cg.ids[j] == 5:
            coord_name = "stagger"
        if cg.ids[j] == 6:
            coord_name = "tilt"
        if cg.ids[j] == 7:
            coord_name = "roll"
        if cg.ids[j] == 8:
            coord_name = "twist"
        if cg.ids[j] == 9:
            coord_name = "shift"
        if cg.ids[j] == 10:
            coord_name = "slide"
        if cg.ids[j] == 11:
            coord_name = "rise"

        x = []
        for z in range(len(unscr_cov_sampled[j])) :
            x.append(z)

        ys = unscr_cov_sampled[j]
        yt = unscr_cov_target[j]

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        # plt.title(r"All: $R=7$ $\phi=0.364$")
        ax.title.set_fontsize(20)
        ax.set_xlabel(r"Sequence",fontsize=20)
        ax.set_ylabel("Cov, diag "+ coord_name,fontsize=20)
        #ax.set_ylim(0,160)
        #ax.set_xlim(-1.2,1.2)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        ax.plot(x, yt, 'b', label="cgna+")
        ax.plot(x, ys, 'r', label="oxdna")
        ax.legend(fontsize = 20)
        plt.savefig("Seq"+str(seqid)+"_STD_"+coord_name+".pdf",bbox_inches='tight',pad_inches=0.05)

        plt.close()

    return


###################################################################################################
############## PRINT FINAL PARAMETERS FILE ########################################################
###################################################################################################

#takes a list with the optimised parameters and the initial SD file, and produces the final SD dep file.
def print_final_pfile(FOPARS,infile) :

    ofile = open("oxDNA_sequence_dependent_parameters_fin.txt",'w')

    ids = np.array(cfun.OPT_PAR_LIST)[:,0]
    if (45 in ids) or (1 in ids) :
        ids = np.append(ids,78)
        ids = np.append(ids,117)

    #update melting as well
    for id_m in cfun.OPT_PAR_LIST_m :
        ids = np.append(ids,id_m[0])


    #add enslaved excluded volume

    for i in range(155:182) :
        ids = np.append(ids,i)

    #Collect all par ids to update ( optimised + dependencies (e.g. continuity) )
    #NOTE: only f1 and f4 continuity are implemented!!

    #f1
    f1_used = [[False,False],[False,False]]

    f1_r0_id = [5,45]
    f1_a_id = [6,46]
    f1_rc_id = [7,47]
    f1_bl_id = [8,48]
    f1_bh_id = [9,49]
    f1_rl_id = [10,50]
    f1_rh_id = [11,51]
    f1_rcl_id = [12,52]
    f1_rch_id = [13,53]

    #f2
    f2_r0_id = [78,117]
    f2_rc_id = [79,118]
    f2_bl_id = [80,119]
    f2_bh_id = [81,120]
    f2_rl_id = [82,121]
    f2_rh_id = [83,122]
    f2_rcl_id = [84,123]
    f2_rch_id = [85,124]

    #f4
    f4_a_id = [15,20,25,30,35,40,55,60,65,87,92,97,102,107,112,126,131,136,141,146,151]
    f4_b_id = [16,21,26,31,36,41,56,61,66,88,93,98,103,108,113,127,132,137,142,147,152]
    f4_ts_id = [17,22,27,32,37,42,57,62,67,89,94,99,104,109,114,128,133,138,143,148,153]
    f4_tc_id = [18,23,28,33,38,43,58,63,68,90,95,100,105,110,115,129,134,139,144,149,154]

    ids_to_update = []

    for ID in ids:

        f1 = False
        f2 = False
        f4 = False

        if ID in ids_to_update:
            continue
        else:
            ids_to_update.append(ID)
            #continuity

            #delta
            if ID == 2:
                ids_to_update.append(3)

            #f1
            idx = None
            if ID in f1_r0_id:
               idx = f1_r0_id.index(ID)
               if f1_used[idx][1]:
                   continue
               else:
                   f1 = True
                   f1_used[idx][0] = True

            if ID in f1_a_id:
               idx = f1_a_id.index(ID)
               if f1_used[idx][0]:
                   continue
               else:
                   f1 = True
                   f1_used[idx][1] = True

            #f2
            if ID in f2_r0_id:
               f2 = True
               idx = f2_r0_id.index(ID)
            #f4
            if ID in f4_a_id:
               f4 = True
               idx = f4_a_id.index(ID)

            if f1:
                ids_to_update.append(f1_rc_id[idx])
                ids_to_update.append(f1_bl_id[idx])
                ids_to_update.append(f1_bh_id[idx])
                ids_to_update.append(f1_rl_id[idx])
                ids_to_update.append(f1_rh_id[idx])
                ids_to_update.append(f1_rcl_id[idx])
                ids_to_update.append(f1_rch_id[idx])

            if f2:
                ids_to_update.append(f2_rc_id[idx])
                ids_to_update.append(f2_bl_id[idx])
                ids_to_update.append(f2_bh_id[idx])
                ids_to_update.append(f2_rl_id[idx])
                ids_to_update.append(f2_rh_id[idx])
                ids_to_update.append(f2_rcl_id[idx])
                ids_to_update.append(f2_rch_id[idx])

            if f4:
                ids_to_update.append(f4_b_id[idx])
                ids_to_update.append(f4_ts_id[idx])
                ids_to_update.append(f4_tc_id[idx])

    #CREATE TENSOR WITH FINAL VALUES OF ALL PARAMETERS
    #we do that on the cpu and copy it to the cpu

    CURR_PARS = torch.tensor(cfun.CURR_PARS_m,device=cfun.device)
    PARS_OPTI = torch.tensor(FOPARS,device=cfun.device)

    CURR_PARS.put_(cfun.UPDATE_MAP, PARS_OPTI)

    #impose symmetries
    VALS = torch.gather( torch.reshape(PARS_OPTI,(-1,)),0,cfun.SYMM_LIST )
    CURR_PARS.put_(cfun.SYMM_LIST_SYMM,VALS)

    #constraints - no continuity

    #delta
    aveD = torch.mean(CURR_PARS[2])
    #print("Old delta:", aveD)
    CURR_PARS[2] = CURR_PARS[2] - aveD + cfun.AVE_DELTA
    #print("New delta:", torch.mean(CURR_PARS[2]))

    #Make STCK_THETA5A average

    aveTH5A = torch.mean(CURR_PARS[60])
    CURR_PARS[60] = aveTH5A

    CURR_PARS[3] = torch.square( CURR_PARS[2] )

    #crst r0
    #crst_r0 = sqrt( stck_r0^2+hydr_r0^2/2*(1+cos(2asin(sqrt(fene_ro^2-stck_r0^2)))) )

    #Constraints - no continuity
    fene_r02_crst = torch.square(torch.clone(CURR_PARS[1][cfun.CRST_TETRA_TYPE_33])+torch.clone(CURR_PARS[1][cfun.CRST_TETRA_TYPE_33_SYMM]))*0.25
    stck_r02_crst = torch.square(torch.clone(CURR_PARS[45][cfun.CRST_TETRA_TYPE_33])+torch.clone(CURR_PARS[45][cfun.CRST_TETRA_TYPE_33_SYMM]))*0.25

    #0.08 <- if hydr_r0 = 0.4. otherwise this has to be changed
    CURR_PARS[78][cfun.CRST_TETRA_TYPE_33] = torch.sqrt( stck_r02_crst+0.08*(1+torch.cos(2*torch.arcsin(0.5*torch.sqrt(fene_r02_crst-stck_r02_crst)))) )
    CURR_PARS[78][cfun.CRST_TETRA_TYPE_33_SYMM] = torch.sqrt( stck_r02_crst+0.08*(1+torch.cos(2*torch.arcsin(0.5*torch.sqrt(fene_r02_crst-stck_r02_crst)))) )

    fene_r02_crst = torch.square(torch.clone(CURR_PARS[1][cfun.CRST_TETRA_TYPE_55])+torch.clone(CURR_PARS[1][cfun.CRST_TETRA_TYPE_55_SYMM]))*0.25
    stck_r02_crst = torch.square(torch.clone(CURR_PARS[45][cfun.CRST_TETRA_TYPE_55])+torch.clone(CURR_PARS[45][cfun.CRST_TETRA_TYPE_55_SYMM]))*0.25

    CURR_PARS[117][cfun.CRST_TETRA_TYPE_55] = torch.sqrt( stck_r02_crst+0.08*(1+torch.cos(2*torch.arcsin(0.5*torch.sqrt(fene_r02_crst-stck_r02_crst)))) )
    CURR_PARS[117][cfun.CRST_TETRA_TYPE_55_SYMM] = torch.sqrt( stck_r02_crst+0.08*(1+torch.cos(2*torch.arcsin(0.5*torch.sqrt(fene_r02_crst-stck_r02_crst)))) )


    #Excluded volume
    #fixing base-base bonded interaction according to the stacking resting distance
    CURR_PARS[171] = CURR_PARS[45] - 0.07 #171 = bonded baba excl volume, 45 = stck r0
    CURR_PARS[cfun.f3_R_ID] = CURR_PARS[cfun.f3_S_ID]*0.97

    # constraints - continuity
    #f1

    CURR_PARS[cfun.f1_RL_ID] = cfun.OFFSET_f1_RL[cfun.OFFSET_f1_ID] + CURR_PARS[cfun.f1_R0_ID]
    CURR_PARS[cfun.f1_RH_ID] = cfun.OFFSET_f1_RH[cfun.OFFSET_f1_ID] + CURR_PARS[cfun.f1_R0_ID]
    CURR_PARS[cfun.f1_RC_ID] = cfun.OFFSET_f1_RC[cfun.OFFSET_f1_ID] + CURR_PARS[cfun.f1_R0_ID]

    EXP1 = torch.exp( -CURR_PARS[cfun.f1_A_ID]*(CURR_PARS[cfun.f1_RL_ID]-CURR_PARS[cfun.f1_R0_ID]) )
    EXP2 = torch.exp( -CURR_PARS[cfun.f1_A_ID]*(CURR_PARS[cfun.f1_RC_ID]-CURR_PARS[cfun.f1_R0_ID]) )
    EXP3 = torch.exp( -CURR_PARS[cfun.f1_A_ID]*(CURR_PARS[cfun.f1_RH_ID]-CURR_PARS[cfun.f1_R0_ID]) )

    CURR_PARS[cfun.f1_BL_ID] = torch.square( CURR_PARS[cfun.f1_A_ID] )*torch.square( EXP1*(1-EXP1) )/( torch.square(1-EXP1) - torch.square(1-EXP2) )
    CURR_PARS[cfun.f1_BH_ID] = torch.square( CURR_PARS[cfun.f1_A_ID] )*torch.square( EXP3*(1-EXP3) )/( torch.square(1-EXP3) - torch.square(1-EXP2) )

    CURR_PARS[cfun.f1_RCL_ID] = CURR_PARS[cfun.f1_RL_ID] - CURR_PARS[cfun.f1_A_ID]/CURR_PARS[cfun.f1_BL_ID]*( EXP1*(1-EXP1) )
    CURR_PARS[cfun.f1_RCH_ID] = CURR_PARS[cfun.f1_RH_ID] - CURR_PARS[cfun.f1_A_ID]/CURR_PARS[cfun.f1_BH_ID]*( EXP3*(1-EXP3) )

    #f2
    CURR_PARS[cfun.f2_RL_ID] = cfun.OFFSET_f2_RL[cfun.OFFSET_f2_ID] + CURR_PARS[cfun.f2_R0_ID]
    CURR_PARS[cfun.f2_RH_ID] = cfun.OFFSET_f2_RH[cfun.OFFSET_f2_ID] + CURR_PARS[cfun.f2_R0_ID]
    CURR_PARS[cfun.f2_RC_ID] = cfun.OFFSET_f2_RC[cfun.OFFSET_f2_ID] + CURR_PARS[cfun.f2_R0_ID]

    TERM1 = CURR_PARS[cfun.f2_RL_ID]-CURR_PARS[cfun.f2_R0_ID]
    TERM2 = CURR_PARS[cfun.f2_RH_ID]-CURR_PARS[cfun.f2_R0_ID]
    TERM3 = CURR_PARS[cfun.f2_RC_ID]-CURR_PARS[cfun.f2_R0_ID]

    CURR_PARS[cfun.f2_RCL_ID] = CURR_PARS[cfun.f2_RL_ID] - (torch.square(TERM1)-torch.square(TERM3))/TERM1
    CURR_PARS[cfun.f2_RCH_ID] = CURR_PARS[cfun.f2_RH_ID] - (torch.square(TERM2)-torch.square(TERM3))/TERM2

    CURR_PARS[cfun.f2_BL_ID] = -0.5*TERM1/( CURR_PARS[cfun.f2_RCL_ID]-CURR_PARS[cfun.f2_RL_ID] )
    CURR_PARS[cfun.f2_BH_ID] = -0.5*TERM2/( CURR_PARS[cfun.f2_RCH_ID]-CURR_PARS[cfun.f2_RH_ID] )


    #f3

    lj_x = torch.clone( torch.pow(CURR_PARS[cfun.f3_S_ID]/CURR_PARS[cfun.f3_R_ID], 6) )

    g1 = 4*( torch.square(lj_x) - lj_x )
    g2 = 12/CURR_PARS[cfun.f3_S_ID]*( 2*torch.square(lj_x)-lj_x )

    CURR_PARS[cfun.f3_RC_ID] = g1/g2
    CURR_PARS[cfun.f3_B_ID] = torch.square(g2)/g1

    #f4
    CURR_PARS[cfun.f4_TS_ID] = torch.sqrt(0.81225/CURR_PARS[cfun.f4_A_ID])
    CURR_PARS[cfun.f4_TC_ID] = 1./CURR_PARS[cfun.f4_A_ID]/CURR_PARS[cfun.f4_TS_ID]
    CURR_PARS[cfun.f4_B_ID] = CURR_PARS[cfun.f4_A_ID]*CURR_PARS[cfun.f4_TS_ID]/(CURR_PARS[cfun.f4_TC_ID]-CURR_PARS[cfun.f4_TS_ID])


    FIN_PARS = torch.tensor(CURR_PARS,device='cpu')

    #PARSE SD IN FILE, COPY WHAT HASN'T CHANGED, UPDATE WHAT HAS CHANGED

    for line in infile.readlines() :
        vals = line.strip().split()
        if len(vals) == 0:
            print(line.strip(),file=ofile)
            continue
        if vals[0][0] == '#':
            print(line.strip(),file=ofile)
            continue
        if vals[0] == "STCK_FACT_EPS":
            print(line.strip(),file=ofile)
            continue

        vals_cn = vals[0].split("_")

        #4D parameters
        if (vals_cn[0] == "STCK" and len(vals_cn) > 3) or vals_cn[0] == "FENE" or (vals_cn[0] == "CRST" and len(vals_cn) > 6):
            if vals_cn[0] == "FENE" and vals_cn[1] == "EPS" : par_name = "FENE_EPS"
            else:
                par_name = vals_cn[0]
                for i in range(1,len(vals_cn)-4):
                    par_name+="_"+vals_cn[i]
            index = PARS_LIST.index(par_name)

            if index in ids_to_update:
                ty3 = base_to_id(vals_cn[len(vals_cn)-1])
                ty2 = base_to_id(vals_cn[len(vals_cn)-2])
                ty1 = base_to_id(vals_cn[len(vals_cn)-3])
                ty0 = base_to_id(vals_cn[len(vals_cn)-4])

                ty = ty0+ty1*4+ty2*16+ty3*64

                print(vals[0] + " = " + str(float(FIN_PARS[index,ty])),file=ofile)
            else:
                print(line.strip(),file=ofile)
        #2D parameters
        else:
            par_name = vals_cn[0]
            for i in range(1,len(vals_cn)-2):
                par_name+="_"+vals_cn[i]

            if par_name == "STCK" or par_name == "HYDR":
                par_name += "_EPS"

            index = PARS_LIST.index(par_name)
            if index in ids_to_update:
                ty2 = base_to_id(vals_cn[len(vals_cn)-1])
                ty1 = base_to_id(vals_cn[len(vals_cn)-2])

                ty = ty1*4+ty2*16

                print(vals[0] + " = " + str(float(FIN_PARS[index,ty])),file=ofile)
            else:
                print(line.strip(),file=ofile)

    ofile.close()


###################################################################################################
############## OPTIMISE ONLY AVERAGE (i.e. no SD) FOR GIVEN COORDINATE ############################
###################################################################################################

#makes gs of specified coordinates (e.g. propeller) flat (i.e. average). Must be called before initialising GPU tensors
def make_it_flat_GS(ids) :

    for l in range(cg.Nseq) :
        ave_sampled = []
        Ncoords = int(len(cfun.internal_coords[l][0])/len(cg.ids))
        #print("check counts: ", Ncoords)

        #make sampled specified coordinates flat (average)
        for i in range(len(ids)):
            ave_sampled.append(0.)
            ave_target.append(0.)

        for n in range(len(cfun.internal_coords[l])):
           for i in range(len(ids)):
               ave_sampled[i] = 0.

           for i in range(len(cfun.internal_coords[l][n])) :
               coord_type = cg.ids[i%len(cg.ids)]
               if coord_type in ids:
                   idx = ids.index(coord_type)
                   ave_sampled[idx] += cfun.internal_coords[l][n][i]/Ncoords

           for i in range(len(cfun.internal_coords[l][n])) :
               coord_type = cg.ids[i%len(cg.ids)]
               if coord_type in ids:
                   idx = ids.index(coord_type)
                   cfun.internal_coords[l][n][i] = ave_sampled[idx]

    #make target coordinates flat
    ave_target = []

    for i in range(len(ids)):
        ave_target.append(0.)

    for l in range(cg.Nseq):
        Ncoords += int(len(cfun.target_mu[l])/len(cg.ids))

    for l in range(cg.Nseq):
        for i in range(len(cfun.target_mu[l])) :
           coord_type = cg.ids[i%len(cg.ids)]
           if coord_type in ids:
               idx = ids.index(coord_type)
               ave_target[idx] += cfun.internal_coords[l][n][i]/Ncoords

    for l in range(cg.Nseq):
        for i in range(len(cfun.target_mu[l])) :
           coord_type = cg.ids[i%len(cg.ids)]
           if coord_type in ids:
               idx = ids.index(coord_type)
               cfun.target_mu[l][i] = ave_target[idx]

    return
