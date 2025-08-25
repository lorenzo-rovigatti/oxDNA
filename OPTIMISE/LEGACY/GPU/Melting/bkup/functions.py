import numpy as np
import math
import parameters_list as parl
import get_cgdna_pars
import config as cg
import cost_function as cfun
import matplotlib.pyplot as plt
import torch
import SantaLucia


PARS_LIST = parl.PARS_LIST
par_index = parl.par_index

bases = ['A','C','G','T','E']
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
    elif  b == 'E':
        return 4
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

    #first, let't check if we used parallel tempering or not!
    #depending on that, the structure of the config file is different

    for line in cfile.readlines() :
        vals = line.split()
        if len(vals) == 0:
            continue
        if vals[0][0] == '#':
            continue

        if vals[0] == "PARALLEL_TEMPERING" :
            if int(vals[1]) == 0 :
                cg.parallel_tempering = True
            else:
                cg.parallel_tempering = False

    if cg.parallel_tempering == True :
        print("Configurations sampled with parallel tempering!")
        print("Trying to read weights and temperatures for each replica of each sequence.")
    else:
        print("No parallel tempering detected. Expected one weight file and one temperature for each sequence.")

    cfile.close()

    cfile = open(cfile_name,'r')

    if cg.first_step_flag == False:
        tfile = open("in_Ts.txt",'r')
        Ts_lines = tfile.readlines()

    checklist = np.zeros(18, dtype = int) #check if some mandatory parameters are missing

    Seq_counter_n5 = -1
    Ts_index_n5 = -1

    Seq_counter_n8 = -1
    Ts_index_n8 = -1

    Seq_counter_n15 = -1
    Ts_index_n15 = -1


    #flags for delta times
    stf = [False, False, False]
    dtf = [False, False, False]
    dpf = [False, False, False]
    dsf = [False, False, False]


    for line in cfile.readlines() :
        vals = line.split()
        if len(vals) == 0:
            continue
        if vals[0][0] == '#':
            continue
        #print(vals)

        ws_n5_read = False
        ws_n8_read = False
        ws_n15_read = False

        #read sequence
        if(vals[0] == 'SEQ'):

            #Seq_counter += 1
            #Ts_index += 1

            if cg.parallel_tempering == False :

                if len(vals) != 9 :
                    print("Invalid SEQ sintax.")
                    print("Usage:")
                    print("SEQ seq box_size Cs simT0 Nts DT opfile wfile")
                    print("box_size = box_size")
                    print("Cs = salt concentration")
                    print("simT0 = simulation temperature first step")
                    print("NTs = number of rew temperatures")
                    print("DT = DT between rew temperatures")
                    print("opfile = order param file")
                    print("wfile = weight_file")
                    return False

            else:

                if len(vals) < 9 :
                    print("Invalid SEQ sintax.")
                    print("Usage:")
                    print("SEQ seq box_size Cs Nreps simTs DT opfile wfiles")
                    print("box_size = box_size")
                    print("Cs = salt concentration")
                    print("Nreps = number of replicas (for PT)")
                    print("simTs = simulation temperatures, one for each replica")
                    print("DT = DT between rew temperatures")
                    print("opfile = order param file")
                    print("wfiles = weight_files, one for each replica")
                    return False


            checklist[0] = 1

            nbp = len(vals[1])

            if nbp == 5:

                cg.good_n5.append(True)

                Seq_counter_n5 += 1
                Ts_index_n5 += 1

                cg.seq_n5.append(vals[1])
                cg.Njuns_n5.append(len(vals[1])-1)

                box = np.zeros(3)
                box[0] = float(vals[2])
                box[1] = float(vals[2])
                box[2] = float(vals[2])

                cg.boxes_n5.append(box)

                Ct = 2./pow(float(vals[2]),3.)*2.6868
                Cs = float(vals[3])

                cg.Ct_n5.append(Ct)
                cg.Cs_n5.append(Cs)

                SLt = SantaLucia.melting_temperature(vals[1],Ct,Cs)
                cfun.target_Tms_n5.append(SLt)

                if cg.parallel_tempering == False:
                    SimT = 0
                    if cg.first_step_flag == True:
                        SimT = float(vals[4])
                    else:
                        SimT = float(Ts_lines[Ts_index])

                    NT = int(vals[5])
                    DT = float(vals[6])

                    cfun.sim_Ts_n5.append(SimT)
                    cfun.n_Ts_n5.append(NT)

                    #DT = DT_new
                    cfun.D_Ts_n5.append(DT)

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
                    if DT*lr+T0<0 :  ##oxdna behaviour at T < 0 is uncertain
                        offset = -(DT*lr+T0-0.5)
                    #if DT*ur+T0>100 :
                    #    offset = - (DT*ur+T0-100) - 0.5

                    for l in range(lr,ur) :
                        rew_Ts.append( DT*l+T0  + offset)

                    cfun.rew_Ts_n5.append(rew_Ts)

                    if ws_n5_read == False:
                        #read weights

                        wfile = open("n5/Seq"+str(Seq_counter_n5)+"/Rep0/"+vals[8],'r')

                        weights = []

                        for line in wfile.readlines() :
                            vals_w = line.split()
                            if len(vals_w) == 0:
                                continue
                            if vals_w[0][0] == '#':
                                continue

                            weights.append(float(vals_w[1]))

                        cfun.weights_n5 = weights
                        wfile.close()
                        ws_n5_read = True


                #w parallel tempering!
                else:
                    Nrs = int(vals[4])
                    #print(Nrs)
                    cg.N_PT_reps_n5 = Nrs
                    sTs = []
                    npars = 4
                    for i in range(Nrs):
                        npars += 1
                        sT = float(vals[npars])
                        sTs.append(sT)

                    npars += 1
                    DT = float(vals[npars])
                    NT = int((sTs[len(sTs)-1]-sTs[0])/DT+1)
                    if NT < 15:
                        NT+=1

                    #print(NT)
                    cfun.sim_Ts_n5.append(sTs)
                    cfun.n_Ts_n5.append(NT)

                    #create T range around T0

                    #T0 = (SimT-SLt)/2.
                    T0 = sTs[0]

                    rew_Ts = []

                    for i in range(NT) :
                        rew_Ts.append(T0+DT*i)

                    #print(rew_Ts)
                    cfun.rew_Ts_n5.append(rew_Ts)


                    npars += 1 #next par is opfile, we don't need it

                    #read weights

                    if ws_n5_read == False:
                        weights = []
                        for i in range(Nrs):
                            weights_1rep = []
                            npars += 1
                            print(vals)
                            print(i, vals[npars])
                            name = "n5/Seq"+str(Seq_counter_n5)+"/Rep0/"+vals[npars]
                            print(name)
                            wfile = open(name,'r')

                            for line in wfile.readlines() :
                                vals_w = line.split()
                                if len(vals_w) == 0:
                                    continue
                                if vals_w[0][0] == '#':
                                    continue

                                weights_1rep.append(float(vals_w[1]))

                            wfile.close()
                            weights.append(weights_1rep)

                        cfun.weights_n5.append(weights)
                        #ws_n5_read = True

            elif nbp == 8:

                cg.good_n8.append(True)

                Seq_counter_n8 += 1
                Ts_index_n8 += 1

                cg.seq_n8.append(vals[1])
                cg.Njuns_n8.append(len(vals[1])-1)

                box = np.zeros(3)
                box[0] = float(vals[2])
                box[1] = float(vals[2])
                box[2] = float(vals[2])

                cg.boxes_n8.append(box)

                Ct = 2./pow(float(vals[2]),3.)*2.6868
                Cs = float(vals[3])

                cg.Ct_n8.append(Ct)
                cg.Cs_n8.append(Cs)

                SLt = SantaLucia.melting_temperature(vals[1],Ct,Cs)
                cfun.target_Tms_n8.append(SLt)

                if cg.parallel_tempering == False:
                    SimT = 0
                    if cg.first_step_flag == True:
                        SimT = float(vals[4])
                    else:
                        SimT = float(Ts_lines[Ts_index])

                    NT = int(vals[5])
                    DT = float(vals[6])

                    cfun.sim_Ts_n8.append(SimT)
                    cfun.n_Ts_n8.append(NT)

                    #DT = DT_new
                    cfun.D_Ts_n8.append(DT)

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
                    if DT*lr+T0<0 :  #oxdna behaviour at T < 0 is uncertain
                        offset = -(DT*lr+T0-0.5)
                    #if DT*ur+T0>100 :
                    #    offset = - (DT*ur+T0-100) - 0.5

                    for l in range(lr,ur) :
                        rew_Ts.append( DT*l+T0  + offset)

                    cfun.rew_Ts_n8.append(rew_Ts)

                    if ws_n8_read == False:
                        #read weights

                        wfile = open("n8/Seq"+str(Seq_counter_n8)+"/Rep0/"+vals[8],'r')

                        weights = []

                        for line in wfile.readlines() :
                            vals_w = line.split()
                            if len(vals_w) == 0:
                                continue
                            if vals_w[0][0] == '#':
                                continue

                            weights.append(float(vals_w[1]))

                        cfun.weights_n8 = weights
                        wfile.close()
                        ws_n8_read = True


                #w parallel tempering!
                else:
                    Nrs = int(vals[4])
                    cg.N_PT_reps_n8 = Nrs
                    sTs = []
                    npars = 4
                    for i in range(Nrs):
                        npars += 1
                        sT = float(vals[npars])
                        sTs.append(sT)

                    npars += 1
                    DT = float(vals[npars])
                    NT = int((sTs[len(sTs)-1]-sTs[0])/DT+1)
                    cfun.sim_Ts_n8.append(sTs)
                    cfun.n_Ts_n8.append(NT)

                    #create T range around T0

                    #T0 = (SimT-SLt)/2.
                    T0 = sTs[0]

                    rew_Ts = []

                    for i in range(NT) :
                        rew_Ts.append(T0+DT*i)

                    cfun.rew_Ts_n8.append(rew_Ts)


                    npars += 1 #opfile, we don't need it

                    #read weights

                    if ws_n8_read == False:
                        weights = []
                        for i in range(Nrs):
                            weights_1rep = []
                            npars += 1
                            name = "n8/Seq"+str(Seq_counter_n8)+"/Rep0/"+vals[npars]
                            wfile = open(name,'r')

                            for line in wfile.readlines() :
                                vals_w = line.split()
                                if len(vals_w) == 0:
                                    continue
                                if vals_w[0][0] == '#':
                                    continue

                                weights_1rep.append(float(vals_w[1]))

                            wfile.close()
                            weights.append(weights_1rep)

                        cfun.weights_n8.append(weights)
                        #ws_n8_read = True


            elif nbp == 15:

                cg.good_n15.append(True)
                Seq_counter_n15 += 1
                Ts_index_n15 += 1

                cg.seq_n15.append(vals[1])
                cg.Njuns_n15.append(len(vals[1])-1)

                box = np.zeros(3)
                box[0] = float(vals[2])
                box[1] = float(vals[2])
                box[2] = float(vals[2])

                cg.boxes_n15.append(box)

                Ct = 2./pow(float(vals[2]),3.)*2.6868
                Cs = float(vals[3])

                cg.Ct_n15.append(Ct)
                cg.Cs_n15.append(Cs)

                SLt = SantaLucia.melting_temperature(vals[1],Ct,Cs)
                cfun.target_Tms_n15.append(SLt)

                if cg.parallel_tempering == False:
                    SimT = 0
                    if cg.first_step_flag == True:
                        SimT = float(vals[4])
                    else:
                        SimT = float(Ts_lines[Ts_index])

                    NT = int(vals[5])
                    DT = float(vals[6])

                    cfun.sim_Ts_n15.append(SimT)
                    cfun.n_Ts_n15.append(NT)

                    #DT = DT_new
                    cfun.D_Ts_n15.append(DT)

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
                    if DT*lr+T0<0 :  #oxdna behaviour at T < 0 is uncertain
                        offset = -(DT*lr+T0-0.5)
                    #if DT*ur+T0>100 :
                    #    offset = - (DT*ur+T0-100) - 0.5

                    for l in range(lr,ur) :
                        rew_Ts.append( DT*l+T0  + offset)

                    cfun.rew_Ts_n15.append(rew_Ts)

                    if ws_n15_read == False:
                        #read weights

                        wfile = open("n15/Seq"+str(Seq_counter_n15)+"/Rep0/"+vals[8],'r')

                        weights = []

                        for line in wfile.readlines() :
                            vals_w = line.split()
                            if len(vals_w) == 0:
                                continue
                            if vals_w[0][0] == '#':
                                continue

                            weights.append(float(vals_w[1]))

                        cfun.weights_n15 = weights
                        wfile.close()
                        ws_n15_read = True


                #w parallel tempering!
                else:
                    Nrs = int(vals[4])
                    cg.N_PT_reps_n15 = Nrs
                    sTs = []
                    npars = 4
                    for i in range(Nrs):
                        npars += 1
                        sT = float(vals[npars])
                        sTs.append(sT)

                    npars += 1
                    DT = float(vals[npars])
                    NT = int((sTs[len(sTs)-1]-sTs[0])/DT+1)
                    cfun.sim_Ts_n15.append(sTs)
                    cfun.n_Ts_n15.append(NT)

                    #create T range around T0

                    #T0 = (SimT-SLt)/2.
                    T0 = sTs[0]

                    rew_Ts = []

                    for i in range(NT) :
                        rew_Ts.append(T0+DT*i)

                    cfun.rew_Ts_n15.append(rew_Ts)


                    npars += 1 #opfile, we don't need it

                    #read weights

                    if ws_n15_read == False:
                        weights = []
                        for i in range(Nrs):
                            weights_1rep = []
                            npars += 1
                            name = "n15/Seq"+str(Seq_counter_n15)+"/Rep0/"+vals[npars]
                            wfile = open(name,'r')

                            for line in wfile.readlines() :
                                vals_w = line.split()
                                if len(vals_w) == 0:
                                    continue
                                if vals_w[0][0] == '#':
                                    continue

                                weights_1rep.append(float(vals_w[1]))

                            wfile.close()
                            weights.append(weights_1rep)

                        cfun.weights_n15.append(weights)
                        #ws_n15_read = True

        if(vals[0] == "SKIP") :
            if vals[1] == "n5" :
                for k in range(2,len(vals)): cg.good_n5[int(vals[k])] = False
            if vals[1] == "n8" :
                for k in range(2,len(vals)): cg.good_n8[int(vals[k])] = False
            if vals[1] == "n15" :
                for k in range(2,len(vals)): cg.good_n15[int(vals[k])] = False

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
            cfun.add_opti_par(vals[1])
            checklist[3] = 1

        if(vals[0] == "SIM_TIME") :
            if vals[1] == "n5" :
                cg.tot_time_n5 = int(float(vals[2]))
                stf[0] = True
            if vals[1] == "n8" :
                cg.tot_time_n8 = int(float(vals[2]))
                stf[1] = True
            if vals[1] == "n15" :
                cg.tot_time_n15 = int(float(vals[2]))
                stf[2] = True
            if stf[0] and stf[1] and stf[2] : checklist[14] = 1

        if(vals[0] == "DELTA_TIME") :
            if vals[1] == "n5":
                cg.delta_time_n5 = int(float(vals[2]))
                dtf[0] = True
            if vals[1] == "n8":
                cg.delta_time_n8 = int(float(vals[2]))
                dtf[1] = True
            if vals[1] == "n15":
                cg.delta_time_n15 = int(float(vals[2]))
                dtf[2] = True
            if dtf[0] and dtf[1] and dtf[2] : checklist[15] = 1

        if(vals[0] == "DELTA_PRINT_EN") :
            if vals[1] == "n5":
                cg.delta_print_en_n5 = int(float(vals[2]))
                dpf[0] = True
            if vals[1] == "n8":
                cg.delta_print_en_n8 = int(float(vals[2]))
                dpf[1] = True
            if vals[1] == "n15":
                cg.delta_print_en_n15 = int(float(vals[2]))
                dpf[2] = True
            if dpf[0] and dpf[1] and dpf[2] : checklist[16] = 1

        if(vals[0] == "DELTA_PRINT_SPLIT_EN") :
            if vals[1] == "n5":
                cg.delta_split_en_n5 = int(float(vals[2]))
                dsf[0] = True
            if vals[1] == "n8":
                cg.delta_split_en_n8 = int(float(vals[2]))
                dsf[1] = True
            if vals[1] == "n15":
                cg.delta_split_en_n15 = int(float(vals[2]))
                dsf[2] = True
            if dsf[0] and dsf[1] and dsf[2] : checklist[17] = 1

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

        if(vals[0] == "read_energy_from_file") :
            if vals[1] == "True" :
                cg.read_energy_from_file = True
            elif vals[1] == "False" :
                cg.read_energy_from_file = False
            else:
                print("Could not understand option for read_energy_from_file.")
                print("Alternatives are True or False")
                print("Setting it to "+ str(cg.read_energy_from_file))
            checklist[10] = 1

        if(vals[0] == "read_coords_from_file") :
            if vals[1] == "True" :
                cg.read_coords_from_file = True
            elif vals[1] == "False" :
                cg.read_coords_from_file = False
            else:
                print("Could not understand option for read_coords_from_file.")
                print("Alternatives are True or False")
                print("Setting it to "+ str(cg.read_coords_from_file))
            checklist[11] = 1

        if(vals[0] == "print_energy_to_file") :
            if vals[1] == "True" :
                cg.print_energy_to_file = True
            elif vals[1] == "False" :
                cg.print_energy_to_file = False
            else:
                print("Could not understand option for print_energy_to_file.")
                print("Alternatives are True or False")
                print("Setting it to "+ str(cg.print_energy_to_file))
            checklist[12] = 1

        if(vals[0] == "print_coords_to_file") :
            if vals[1] == "True" :
                cg.print_coords_to_file = True
            elif vals[1] == "False" :
                cg.print_coords_to_file = False
            else:
                print("Could not understand option for print_coords_to_file.")
                print("Alternatives are True or False")
                print("Setting it to "+ str(cg.print_coords_to_file))
            checklist[13] = 1


    cfile.close()
    if cg.first_step_flag == False:
        tfile.close()

    #CHECK AND PRINT

    if checklist[0] == 1:

        cg.Nseq_n5 = len(cg.seq_n5)

        print("NBPS = 5")
        for i in range(cg.Nseq_n5):
                print("SEQUENCE " + str(i) + ": "+cg.seq_n5[i])
                print("Njunctions Seq" + str(i) + ": "+str(cg.Njuns_n5[i]))
                print("Single strand concentration: " + str(cg.Ct_n5[i]) + " M")
                print("Salt concentration: " + str(cg.Cs_n5[i]) + " M")
                print("Sim temperature: " + str(cfun.sim_Ts_n5[i]) + " M")
                if cg.good_n5[i] == False : print("EXCLUDED")


                string = ""

                for l in range(len(cfun.rew_Ts_n5[i])) :

                    string = string + str(cfun.rew_Ts_n5[i][l]) + " "

                print("Reweighted temperatures: " + string)

                string = ""
                """
                for l in range(len(cfun.weights_n5)) :
                    string = string + str(cfun.weights_n5[l]) + " "

                print("Weights: " +  string)
                """

        cg.Nseq_n8 = len(cg.seq_n8)

        print("NBPS = 8")
        for i in range(cg.Nseq_n8):
                print("SEQUENCE " + str(i) + ": "+cg.seq_n8[i])
                print("Njunctions Seq" + str(i) + ": "+str(cg.Njuns_n8[i]))
                print("Single strand concentration: " + str(cg.Ct_n8[i]) + " M")
                print("Salt concentration: " + str(cg.Cs_n8[i]) + " M")
                print("Sim temperature: " + str(cfun.sim_Ts_n8[i]) + " M")
                if cg.good_n8[i] == False : print("EXCLUDED")

                string = ""

                for l in range(len(cfun.rew_Ts_n8[i])) :

                    string = string + str(cfun.rew_Ts_n8[i][l]) + " "

                print("Reweighted temperatures: " + string)

                string = ""

                """
                for l in range(len(cfun.weights_n8)) :
                    string = string + str(cfun.weights_n8[l]) + " "

                print("Weights: " +  string)
                """

        cg.Nseq_n15 = len(cg.seq_n15)

        print("NBPS = 15")
        for i in range(cg.Nseq_n15):
                print("SEQUENCE " + str(i) + ": "+cg.seq_n15[i])
                print("Njunctions Seq" + str(i) + ": "+str(cg.Njuns_n15[i]))
                print("Single strand concentration: " + str(cg.Ct_n15[i]) + " M")
                print("Salt concentration: " + str(cg.Cs_n15[i]) + " M")
                print("Sim temperature: " + str(cfun.sim_Ts_n15[i]) + " M")
                if cg.good_n15[i] == False : print("EXCLUDED")

                string = ""

                for l in range(len(cfun.rew_Ts_n15[i])) :

                    string = string + str(cfun.rew_Ts_n15[i][l]) + " "

                print("Reweighted temperatures: " + string)

                string = ""

                """
                for l in range(len(cfun.weights_n15)) :
                    string = string + str(cfun.weights_n15[l]) + " "

                print("Weights: " +  string)
                """

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

        print("PARAMETERS used in the optimisation: ")
        print("# id type (base 10)")
        for l in range(len(cfun.OPT_PAR_LIST)) :
            print(cfun.OPT_PAR_LIST[l])

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

    if checklist[10] == 1:
        print("Read energy from file = "+str(cg.read_energy_from_file))
    else :
        print("OPTION. Reading energy from trajectory")
        print("If you want to read from file:")
        print("Usage:")
        print("read_energy_from_file True")

    if checklist[11] == 1:
        print("Read coords from file = "+str(cg.read_coords_from_file))
    else :
        print("OPTION. Reading coords from trajectory")
        print("If you want to read from file:")
        print("Usage:")
        print("read_coords_from_file True")

    if checklist[12] == 1:
        print("Print energy to file = "+str(cg.print_energy_to_file))
    else :
        print("OPTION. Printing energy to file disabled")
        print("If you want to print to file:")
        print("Usage:")
        print("print_energy_to_file True")

    if checklist[13] == 1:
        print("Print coords to file = "+str(cg.print_coords_to_file))
    else :
        print("OPTION. Printing coordinates to file disabled")
        print("If you want to print to file:")
        print("Usage:")
        print("print_coords_to_file True")

        if checklist[14] == 1:
            print("Total simulation time n5 = "+str(cg.tot_time_n5))
            print("Total simulation time n8 = "+str(cg.tot_time_n8))
            print("Total simulation time n15 = "+str(cg.tot_time_n15))
        else :
            print("MANDATORY. Total simulation time")
            print("Usage:")
            print("SIM_TIME nx tot_time")
            return False

    if checklist[15] == 1:
        print("n5")
        print("Delta t between snapshots = "+str(cg.delta_time_n5))
        cg.Nconfs_per_pt_rep_n5 = int(cg.tot_time_n5/cg.delta_time_n5)-cg.in_snap
        cg.tot_Nconfs_per_pt_rep_n5 = cg.Nreps*cg.Nconfs_per_pt_rep_n5
        print("Number of sampled configurations per replica = "+str(cg.Nconfs_per_pt_rep_n5))
        print("Total number of sampled configurations per replica (accounting for Nreps simulations)= "+str(cg.tot_Nconfs_per_pt_rep_n5))
        print("n8")
        print("Delta t between snapshots = "+str(cg.delta_time_n8))
        cg.Nconfs_per_pt_rep_n8 = int(cg.tot_time_n8/cg.delta_time_n8)-cg.in_snap
        cg.tot_Nconfs_per_pt_rep_n8 = cg.Nreps*cg.Nconfs_per_pt_rep_n8
        print("Number of sampled configurations per replica = "+str(cg.Nconfs_per_pt_rep_n8))
        print("Total number of sampled configurations per replica (accounting for Nreps simulations)= "+str(cg.tot_Nconfs_per_pt_rep_n8))
        print("n15")
        print("Delta t between snapshots = "+str(cg.delta_time_n15))
        cg.Nconfs_per_pt_rep_n15 = int(cg.tot_time_n15/cg.delta_time_n15)-cg.in_snap
        cg.tot_Nconfs_per_pt_rep_n15 = cg.Nreps*cg.Nconfs_per_pt_rep_n15
        print("Number of sampled configurations per replica = "+str(cg.Nconfs_per_pt_rep_n15))
        print("Total number of sampled configurations per replica (accounting for Nreps simulations)= "+str(cg.tot_Nconfs_per_pt_rep_n15))
    else :
        print("MANDATORY. Delta t between snapshots")
        print("Usage:")
        print("DELTA_TIME nx snap_time")
        return False

    if checklist[16] == 1:
        print("Delta t between printing energy n5 = "+str(cg.delta_print_en_n5))
        print("Delta t between printing energy n8 = "+str(cg.delta_print_en_n8))
        print("Delta t between printing energy n15 = "+str(cg.delta_print_en_n15))
    else :
        print("MANDATORY. Delta t between printing energy")
        print("Usage:")
        print("DELTA_PRINT_EN nx en_print_time")
        return False

    if checklist[17] == 1:
        print("Delta t between printing split energy n5 = "+str(cg.delta_print_en_n5))
        print("Delta t between printing split energy n8 = "+str(cg.delta_print_en_n8))
        print("Delta t between printing split energy n15 = "+str(cg.delta_print_en_n15))
    else :
        print("MANDATORY. Delta t between printing split energy")
        print("Usage:")
        print("DELTA_PRINT_SPLIT_EN nx split_en_print_time")
        return False


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
            for i in range(5):
                for j in range(5):
                    ty=i+base_to_id(vals_cn[1])*5+base_to_id(vals_cn[2])*5*5+j*5*5*5    #from base 4 to base 10
                    oi = [PARS_LIST.index("STCK_EPS"),ty]
                    over_indices.append(oi)
                    over_vals.append(float(vals[2]))
        elif vals_cn[0] == "HYDR" and len(vals_cn) == 3:
            for i in range(5):
                for j in range(5):
                    ty=i+base_to_id(vals_cn[1])*5+base_to_id(vals_cn[2])*5*5+j*5*5*5    #from base 4 to base 10
                    oi = [PARS_LIST.index("HYDR_EPS"),ty]
                    over_indices.append(oi)
                    over_vals.append(float(vals[2]))
        elif vals_cn[0] == "HYDR" or vals_cn[0] == "CRST" :
            par_name = vals_cn[0]
            if vals_cn[0] == "HYDR" or vals_cn[1] == "K":
                for i in range(1,len(vals_cn)-2):
                    par_name += "_"+vals_cn[i]
                for i in range(5):
                    for j in range(5):
                        ty=i+base_to_id(vals_cn[len(vals_cn)-2])*5+base_to_id(vals_cn[len(vals_cn)-1])*5*5+j*5*5*5    #from base 4 to base 10
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
                    for i in range(5):
                        for j in range(5):
                            ty=i+base_to_id(vals_cn[len(vals_cn)-2])*5+base_to_id(vals_cn[len(vals_cn)-1])*5*5+j*5*5*5    #from base 4 to base 10
                            oi = [PARS_LIST.index(par_name),ty]
                            over_indices.append(oi)
                            over_vals.append(float(vals[2]))
            else :
                for i in range(1,len(vals_cn)-4):
                    par_name += "_"+vals_cn[i]
                ty=base_to_id(vals_cn[len(vals_cn)-4])+base_to_id(vals_cn[len(vals_cn)-3])*5+base_to_id(vals_cn[len(vals_cn)-2])*5*5+base_to_id(vals_cn[len(vals_cn)-1])*5*5*5
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
                         ty=base_to_id(vals_cn[len(vals_cn)-4])+base_to_id(vals_cn[len(vals_cn)-3])*5+base_to_id(vals_cn[len(vals_cn)-2])*5*+base_to_id(vals_cn[len(vals_cn)-1])*5*5*5
                         oi = [PARS_LIST.index(par_name),ty]
                         over_indices.append(oi)
                         over_vals.append(float(vals[2]))

        elif vals_cn[0] == "STCK" or vals_cn[0] == "FENE" or vals_cn[0] == "EXCL":
            par_name = vals_cn[0]
            for i in range(1,len(vals_cn)-4):
                par_name += "_"+vals_cn[i]
            ty=base_to_id(vals_cn[len(vals_cn)-4])+base_to_id(vals_cn[len(vals_cn)-3])*5+base_to_id(vals_cn[len(vals_cn)-2])*5*5+base_to_id(vals_cn[len(vals_cn)-1])*5*5*5    #from base 4 to base 10
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

    OXPS_zero = np.zeros((len(PARS_LIST),625),dtype=float)
    shifts = np.zeros((2,625),dtype=float) #0 = hydr, 1 = stck

    for i in range(len(PARS_LIST)) :

        if PARS_LIST[i] == "STCK_EPS":
            for j in range(625) :
                OXPS_zero[i][j] = 1.
            continue
        if PARS_LIST[i] == "HYDR_EPS":
            for j in range(625) :
                OXPS_zero[i][j] = 0.
            continue

        index = pars_mh.index(PARS_LIST[i])
        val = vals_mh[index]

        for j in range(625) :
            OXPS_zero[i][j] = val

    #here we use the initial custom parameters, must include stck_x_y and hydr_x_y
    if [4,75] not in over_indices :
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
    for j in range(625) :
        #hydr
        shifts[0][j] = Morse(OXPS_zero[par_index[7]][j],OXPS_zero[par_index[4]][j],OXPS_zero[par_index[5]][j],OXPS_zero[par_index[6]][j])
        #stacking
        OXPS_zero[par_index[44]][j] =  OXPS_zero[par_index[44]][j]* (1.0 - stck_fact_eps + (T * 9.0 * stck_fact_eps))
        shifts[1][j] = Morse(OXPS_zero[par_index[47]][j],OXPS_zero[par_index[44]][j],OXPS_zero[par_index[45]][j],OXPS_zero[par_index[46]][j])

    return OXPS_zero, shifts

###################################################################################################
############## READ TRAJECTORY AND COMPUTE OXDNA COORDINATES (i.e angles and distances) ###########
###################################################################################################

#NOTE: FOR THE MELTING TEMPERATURE< WE NEED THIESE ONLY TO COMPUTE THE INITIAL ENERGY!

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

    #we are not changing crst, so we don't have to increase the range
    rcut_high = rcut_high
    rcut_low = rcut_low

    rcut_dh_n5 = []

    for i in range(len(cfun.dh_l_n5)):
        if cg.parallel_tempering:
           rcut_dh_n5.append(4.5*max(cfun.dh_l_n5[i])*1.04) #1.04 accounts for how much lamda would change for an extra 20K
           #0.962 ~ 2*sqrt((POS_MM_BACK1) * (POS_MM_BACK1) + (POS_MM_BACK2) * (POS_MM_BACK2))

        else: rcut_dh_n5.append(4.5*cfun.dh_l_n5[i]*1.07)

    rcut_dh_n8 = []

    for i in range(len(cfun.dh_l_n8)):
        if cg.parallel_tempering:
           rcut_dh_n8.append(4.5*max(cfun.dh_l_n8[i])*1.04) #1.04 accounts for how much lamda would change for an extra 20K
           #0.962 ~ 2*sqrt((POS_MM_BACK1) * (POS_MM_BACK1) + (POS_MM_BACK2) * (POS_MM_BACK2))

        else: rcut_dh_n8.append(4.5*cfun.dh_l_n8[i]*1.07)

    rcut_dh_n15 = []

    for i in range(len(cfun.dh_l_n15)):
        if cg.parallel_tempering:
           rcut_dh_n15.append(4.5*max(cfun.dh_l_n15[i])*1.04) #1.04 accounts for how much lamda would change for an extra 20K
           #0.962 ~ 2*sqrt((POS_MM_BACK1) * (POS_MM_BACK1) + (POS_MM_BACK2) * (POS_MM_BACK2))

        else: rcut_dh_n15.append(4.5*cfun.dh_l_n15[i]*1.07)


    return rcut_low , rcut_high, rcut_dh_n5, rcut_dh_n8, rcut_dh_n15

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
#pos_stck = [0.34,0.34,0.34,0.34]
#pos_hydr = [0.4,0.4,0.4,0.4]
#note: oxdna3
pos_stck = [0.37,0.37,0.37,0.37]
pos_hydr = [0.43,0.37,0.43,0.37]
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


def read_oxdna_trajectory_dist_and_angles(rcut_low, rcut_high, rcut_dh, tr_file, topo_file, box):

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
    types_unbn_33 = []
    types_unbn_55 = []

    dh_r = []
    dh_ty = []
    dh_chcut = []

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
    counts_b = -1
    config = []
    nid = 0
    for line in tr_file.readlines():
        a = line.strip()[0]
        if counts <= cg.in_snap :
            if a == 't': counts += 1
            continue
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
                        ty0 = 4
                        ty1 = 4
                        ty2 = 4
                        ty3 = 4
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
                                   break
                               else:
                                   for z in range(len(topology)) :
                                       if topology[z].id == topology[j].up_id:
                                          ty3 = base_to_id(topology[z].base_type)
                                          break
                               break

                        ty = ty0+ty1*5+ty2*5*5+ty3*5*5*5 #tetramer type in base 10

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
                    types_unbn_1conf_33 = []
                    types_unbn_1conf_55 = []

                    dh_r_1conf = []
                    dh_ty_1conf = []
                    dh_chcut_1conf = []


                    #UNDBONDED
                    #out_string = ""
                    for i in range(len(topology)) :
                        #out_string += " " + str(topology[i].id)
                        ty0_33 = 4
                        ty0_55 = 4
                        ty1 = 4
                        ty2 = 4
                        ty3_33 = 4
                        ty3_55 = 4
                        n1 = config[i]
                        ty1 = base_to_id(topology[i].base_type)
                        for z in range(len(topology)) :
                           if topology[z].id == topology[i].down_id:
                              ty0_55 = base_to_id(topology[z].base_type)
                           if topology[z].id == topology[i].up_id:
                              ty0_33 = base_to_id(topology[z].base_type)

                        for j in range(len(topology)) :

                            #min image to take care of duplexes with strands on different images
                            #this can happen in melting simulations since a duplex can melt and then reform
                            min_im_v = np.zeros(3)

                            min_im_v[0] = round((n1.c[0]-n2.c[0])/box[0])*box[0]
                            min_im_v[1] = round((n1.c[1]-n2.c[1])/box[1])*box[1]
                            min_im_v[2] = round((n1.c[2]-n2.c[2])/box[2])*box[2]

                            n2 = config[j]

                            out_of_range = False

                            hyv = n1.hydr - n2.hydr - min_im_v
                            hydist = np.linalg.norm(hyv)
                            hydist2 = hydist*hydist

                            bbv = n1.bb - n2.bb - min_im_v
                            bbdist = np.linalg.norm(bbv)
                            bbdist2 = bbdist*bbdist

                            if bbdist > rcut_dh and rcut_dh > rcut_high: continue #if we are out of debye huckle range
                            if hydist2 < rcut_sq_low: out_of_range = True #verlet cutoff
                            if hydist2 > rcut_sq_high: out_of_range = True #verlet cutoff
                            if topology[j].id <= topology[i].id: continue #ordered ids (avoid pairs repetition)
                            if topology[j].id == topology[i].down_id or topology[j].id == topology[i].up_id: continue #no bonded pairs

                            #out_string += str(topology[j].id)

                            if out_of_range and cg.debye_huckel == False: continue #if we dont't use debye huckel, continue.

                            ty2 = base_to_id(topology[j].base_type)

                            # if topology[j].up_id == -1:
                            #   ty3_33 = 0
                            #if topology[j].down_id == -1:
                               #ty3_55 = 0

                            if topology[j].down_id != -1:
                                for z in range(len(topology)) :
                                   if topology[z].id == topology[j].down_id:
                                      ty3_55 = base_to_id(topology[z].base_type)
                                      break
                            if topology[j].up_id != -1:
                                for z in range(len(topology)) :
                                   if topology[z].id == topology[j].up_id:
                                      ty3_33 = base_to_id(topology[z].base_type)
                                      break

                            #if topology[i].down_id == -1:
                            #   ty0_55 = 0
                            #if topology[i].up_id == -1:
                               #ty0_33 = 0

                            ty_33 = ty0_33+ty1*5+ty2*5*5+ty3_33*5*5*5 #tetramer type in base 10
                            ty_55 = ty0_55+ty1*5+ty2*5*5+ty3_55*5*5*5 #tetramer type in base 10


                            if cg.debye_huckel:
                                chcut = 1
                                if ty0_33 == 4 or ty0_55 == 4: chcut*=2
                                if ty3_33 == 4 or ty3_55 == 4: chcut*=2

                                dh_r_1conf.append(bbdist)
                                dh_ty_1conf.append(ty_33)
                                dh_chcut_1conf.append(chcut)

                            if out_of_range == False:
                                #compute unbnd pair coordinates
                                hydr_r_1conf.append(hydist)
                                rhydr = hyv/hydist

                                th1_1conf.append(np.arccos(-np.dot(n1.bv,n2.bv)))
                                th3_1conf.append(np.arccos(-np.dot(n1.bv,rhydr)))
                                th2_1conf.append(np.arccos(np.dot(n2.bv,rhydr)))

                                th4_unbn_1conf.append(np.arccos(np.dot(n1.n,n2.n)))
                                th8_1conf.append(np.arccos(-np.dot(n1.n,rhydr)))
                                th7_1conf.append(np.arccos(np.dot(n2.n,rhydr)))

                                types_unbn_1conf_33.append(ty_33)
                                types_unbn_1conf_55.append(ty_55)


                    #print(len(types_unbn_1conf_33), out_string)
                    types_unbn_33.append(types_unbn_1conf_33)
                    types_unbn_55.append(types_unbn_1conf_55)

                    hydr_r.append(hydr_r_1conf)

                    th1.append(th1_1conf)
                    th2.append(th2_1conf)
                    th3.append(th3_1conf)

                    th4_unbn.append(th4_unbn_1conf)
                    th7.append(th7_1conf)
                    th8.append(th8_1conf)

                    dh_r.append(dh_r_1conf)
                    dh_ty.append(dh_ty_1conf)
                    dh_chcut.append(dh_chcut_1conf)

    return fene_r, stck_r, th4_bn, th5, th6, cosphi1, cosphi2, types_bn, hydr_r, th1, th2, th3, th4_unbn, th7, th8, types_unbn_33, types_unbn_55, dh_r, dh_ty, dh_chcut

###################################################################################################
############## CALLBACK FUNCTIONS #############################################################
###################################################################################################

#change temperature to reweighting target in input2.an
"""
def update_T_input_file(T,file_name) :

    strold1 = r"(T) = .*(K)$"  #could be K or
    strold2 = r"(T) = .*(C)$"  #could be C
    strnew = "T = "+str(T)+"C"

    pysed.replace(strold1,strnew,file_name)
    pysed.replace(strold2,strnew,file_name)

    return
"""

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



###################################################################################################
############## PRINT FINAL PARAMETERS FILES ########################################################
###################################################################################################

def print_final_melting_temperatures(mTs) :

    ofile = open("final_rew_mTs.txt",'w')
    for i in range(len(mTs)) :
        if i%cg.Nreps == 0:
            print(mTs[i], file=ofile)
    ofile.close()

    return


#takes a list with the optimised parameters and the initial SD file, and produces the final SD dep file.
def print_final_pfile(FOPARS,infile) :

    ofile = open("oxDNA_sequence_dependent_parameters_fin.txt",'w')

    ids = np.array(cfun.OPT_PAR_LIST)[:,0]
    if (45 in ids) or (1 in ids) :
        ids = np.append(ids,78)
        ids = np.append(ids,117)


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

    CURR_PARS = torch.tensor(cfun.PAR0,device=cfun.device)
    PARS_OPTI = torch.tensor(FOPARS,device=cfun.device)

    CURR_PARS.put_(cfun.UPDATE_MAP, PARS_OPTI)

    #impose symmetries
    VALS = torch.gather( torch.reshape(PARS_OPTI,(-1,)),0,cfun.SYMM_LIST )
    CURR_PARS.put_(cfun.SYMM_LIST_SYMM,VALS)

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
        if (vals_cn[0] == "STCK" and len(vals_cn) > 3) or vals_cn[0] == "FENE" or vals_cn[0] == "EXCL" or (vals_cn[0] == "CRST" and len(vals_cn) > 6):
            par_name = vals_cn[0]
            for i in range(1,len(vals_cn)-4):
                par_name+="_"+vals_cn[i]
            index = PARS_LIST.index(par_name)

            if index in ids_to_update:
                ty3 = base_to_id(vals_cn[len(vals_cn)-1])
                ty2 = base_to_id(vals_cn[len(vals_cn)-2])
                ty1 = base_to_id(vals_cn[len(vals_cn)-3])
                ty0 = base_to_id(vals_cn[len(vals_cn)-4])

                ty = ty0+ty1*5+ty2*25+ty3*125

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

                ty = ty1*5+ty2*25

                print(vals[0] + " = " + str(float(FIN_PARS[index,ty])),file=ofile)
            else:
                print(line.strip(),file=ofile)

    ofile.close()

#takes a list with the optimised parameters and the initial SD file, and produces the final SD dep file.
def print_final_pfile_AllFromOpt(FOPARS,infile) :

    ofile = open("oxDNA_sequence_dependent_parameters_fin_opt.txt",'w')

    #CREATE TENSOR WITH FINAL VALUES OF ALL PARAMETERS
    #we do that on the cpu and copy it to the cpu

    CURR_PARS = torch.tensor(cfun.PAR0,device=cfun.device)
    PARS_OPTI = torch.tensor(FOPARS,device=cfun.device)

    CURR_PARS.put_(cfun.UPDATE_MAP, PARS_OPTI)

    #impose symmetries
    VALS = torch.gather( torch.reshape(PARS_OPTI,(-1,)),0,cfun.SYMM_LIST )
    CURR_PARS.put_(cfun.SYMM_LIST_SYMM,VALS)

    FIN_PARS = torch.tensor(CURR_PARS,device='cpu')

    #PARSE SD IN FILE UPDATE EVERYTHING

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
        if (vals_cn[0] == "STCK" and len(vals_cn) > 3) or vals_cn[0] == "FENE" or vals_cn[0] == "EXCL" or (vals_cn[0] == "CRST" and len(vals_cn) > 6):
            par_name = vals_cn[0]
            for i in range(1,len(vals_cn)-4):
                par_name+="_"+vals_cn[i]
            index = PARS_LIST.index(par_name)

            ty3 = base_to_id(vals_cn[len(vals_cn)-1])
            ty2 = base_to_id(vals_cn[len(vals_cn)-2])
            ty1 = base_to_id(vals_cn[len(vals_cn)-3])
            ty0 = base_to_id(vals_cn[len(vals_cn)-4])

            ty = ty0+ty1*5+ty2*25+ty3*125

            print(vals[0] + " = " + str(float(FIN_PARS[index,ty])),file=ofile)

        #2D parameters
        else:
            par_name = vals_cn[0]
            for i in range(1,len(vals_cn)-2):
                par_name+="_"+vals_cn[i]

            if par_name == "STCK" or par_name == "HYDR":
                par_name += "_EPS"

            index = PARS_LIST.index(par_name)

            ty2 = base_to_id(vals_cn[len(vals_cn)-1])
            ty1 = base_to_id(vals_cn[len(vals_cn)-2])

            ty = ty1*5+ty2*25

            print(vals[0] + " = " + str(float(FIN_PARS[index,ty])),file=ofile)


    ofile.close()



###################################################################################################
############## OPTIMISE ONLY AVERAGE (i.e. no SD) FOR GIVEN COORDINATE ############################
###################################################################################################

#makes gs of specified coordinates (e.g. propeller) flat (i.e. average). Must be called before initialising GPU tensors
def make_it_flat_GS(ids) :

    for l in range(cg.Nseq) :
        ave_sampled = []
        ave_target = []
        Ncoords = int(len(cfun.internal_coords[l][0])/len(cg.ids))
        #print("check counts: ", Ncoords)

        #make sampled specified coordinates flat (average)
        for i in range(len(ids)):
            ave_sampled.append(0.)
            ave_target.append(0.)

        for n in range(len(cfun.internal_coords[l])):
           for i in range(len(ids)):
               ave_sampled[i] = 0.
               ave_target[i] = 0.

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
        for i in range(len(cfun.target_mu[l])) :
           coord_type = cg.ids[i%len(cg.ids)]
           if coord_type in ids:
               idx = ids.index(coord_type)
               ave_target[idx] += cfun.internal_coords[l][n][i]/Ncoords

        for i in range(len(cfun.target_mu[l])) :
           coord_type = cg.ids[i%len(cg.ids)]
           if coord_type in ids:
               idx = ids.index(coord_type)
               cfun.target_mu[l][i] = ave_target[idx]
