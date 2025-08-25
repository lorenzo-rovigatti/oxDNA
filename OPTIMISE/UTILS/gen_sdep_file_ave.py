#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:29:49 2023

@author: yqb22156
"""


import numpy as np
import sys
import os
program_path=os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, program_path+'/..')
from functions import impose_continuity
import config as cg
current_path=os.getcwd()

def gen_sd_files_ave(par, opti) :
    
    string = "cp "+program_path+"/oxDNA_sequence_dependent_parameters.txt "+current_path
    os.system(string)
      
    ifile = open('oxDNA_sequence_dependent_parameters.txt','r')
    ofile = open('oxDNA_sequence_dependent_parameters_in.txt','w') 

    #ofile = open('oxDNA_sequence_dependent_parameters_in.txt','w')

    for line in ifile.readlines() :
        print(line.strip('\n'), file=ofile)
        
    print('\n', file=ofile)
        
    ifile.close()
    ofile.close()
    
    ofile1 = open('oxDNA_sequence_dependent_parameters.txt','a')    #file with fixed sequence dependent parameters
    ofile2 = open('oxDNA_sequence_dependent_parameters_in.txt','a') 
    
    for i in range(len(cg.par_codename)) :
        
        if cg.par_codename[i] == "FENE_DELTA":
            if opti[i] == False:
                print(cg.par_codename[i]+" = "+str(par[i]),file=ofile1)
                print(cg.par_codename[i]+"2 = "+str(par[i]*par[i]),file=ofile1)
                print('\n',file=ofile1)
            print(cg.par_codename[i]+" = "+str(par[i]),file=ofile2)
            print(cg.par_codename[i]+"2 = "+str(par[i]*par[i]),file=ofile2)
            print('\n',file=ofile2)
            continue
        
        for l in range(4):
            for m in range(4) :
                if opti[i] == False:
                    print(cg.par_codename[i]+"_"+cg.bases[l]+"_"+cg.bases[m]+" = "+str(par[i]),file=ofile1)
                print(cg.par_codename[i]+"_"+cg.bases[l]+"_"+cg.bases[m]+" = "+str(par[i]),file=ofile2)
        if opti[i] == False:
            print('\n',file=ofile1)
        print('\n',file=ofile2)
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
                        if opti == False:
                            print(names[k]+"_"+cg.bases[l]+"_"+cg.bases[m]+" = "+str(output[k+1]),file=ofile1) 
                        print(names[k]+"_"+cg.bases[l]+"_"+cg.bases[m]+" = "+str(output[k+1]),file=ofile2)
                if opti == False:
                    print('\n',file=ofile1)
                print('\n',file=ofile2)
                    
    ofile1.close()
    ofile2.close()
    
    return

def write_to_cnames(vals,ofile) :
    
    if vals[0] != 'cname' :
        return    
    
    line = vals[0]
    
    nvals = vals[1].split("_")
    
    app = " " +vals[2]
    
    if(vals[2]=='CONTINUITY') :
        app += " " + vals[3]
        
    if nvals[0] == 'STCK' or nvals[0] == 'FENE':
            line = vals[0]+" " + vals[1] + "_A_A" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_A_T" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_A_G" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_A_C" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_T_A" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_T_G" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_T_C" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_G_C" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_C_G" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_C_C" + app
            print(line,file=ofile)            
            
    elif nvals[0] == 'HYDR' :
            line = vals[0] + " " + vals[1] + "_A_T" + app
            print(line,file=ofile)
            line = vals[0] + " " + vals[1] + "_C_G" + app
            print(line,file=ofile)   
            
    elif nvals[0] == 'CRST':
            line = vals[0]+" " + vals[1] + "_A_A" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_A_T" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_A_G" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_A_C" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_T_T"+ app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_T_G" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_T_C" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_G_G" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_G_C" + app
            print(line,file=ofile)
            line = vals[0]+" " + vals[1] + "_C_C" + app
            print(line,file=ofile)
            
     
     
    return

if len(sys.argv) != 2:
    print("Unknown argument format.")
    print("Usage: " + sys.argv[0] +" gen_codenames_file")
    sys.exit(1)

iname = sys.argv[1]

ifile = open(iname,'r')
pfile = open("SD_par_codenames.txt",'w')

opti = []
par = []

for line in ifile.readlines() :
    vals = line.split()
    if len(vals) == 0:
        continue
    if vals[0][0] == '#':
        continue        
    #print(vals)
    
    if len(vals) < 4 :
        print("Not enough entries for parameter " + vals[1])
        print("Value missing?")
        print("Usage: ")
        print("cname/sname "+str(vals[0])+" FIXED/OPTIMISE/CONTINUITY value")
        print("cname for OPTIMISE")
        print("sname for FIXED")
        sys.exit(1)

    #read parameters to use for optimisation
    if(vals[0] == "cname") :
        #print(vals)
        if vals[2] == 'OPTIMISE' :
            cg.par_codename.append(vals[1])
            par.append(float(vals[3]))
            opti.append(True)            
        else :
            cg.continuity_par_codenames.append(vals[1])
            cg.continuity_par_values.append(float(vals[3]))
            
    if(vals[0] == "sname") :
        #print(vals)
        if vals[2] == 'FIXED' :
            cg.par_codename.append(vals[1])
            par.append(float(vals[3]))
            opti.append(False)
        else :
            cg.continuity_par_codenames.append(vals[1])
            cg.continuity_par_values.append(float(vals[3]))
            
    write_to_cnames(vals,pfile)
            
ifile.close()
pfile.close()
 
cg.par_dimension = len(cg.par_codename)
for i in range(cg.par_dimension) :
    cg.used.append(False)

#place fixed par on top
opti_ord = []
par_ord = []
for i in range(len(opti)) :
    if opti[i] == False:
        opti_ord.append(opti[i])
        par_ord.append(par[i])
        
for i in range(len(opti)) :
    if opti[i] == True:
        opti_ord.append(opti[i])
        par_ord.append(par[i])
        
gen_sd_files_ave(par_ord, opti_ord)