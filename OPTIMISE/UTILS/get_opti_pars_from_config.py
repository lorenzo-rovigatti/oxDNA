#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:45:37 2024

@author: AB

Given a par file and a conf file, this program loops through the par file and extracts the optimised parameters
Note that symmetries and continuity are accounted for. 
"""

import sys

symm_stck = True
bases = ['A','C','G','T']

all_fixed_pars = []

#given a parameter, check if continuity constraints should be imposed
def impose_continuity(par) :
    vals = par.split('_')
    output = []
    
    f1 = False
    f2 = False
    f4 = False
    
    
    #this is the scaling factor of the HYDRO and STCK (e.g. HYDR_A_T and STCK_G_A)
    if len(vals) == 3 and (vals[0] == "HYDR" or vals[0] == "STCK") and (vals[1] in bases) and (vals[2] in bases) :
        output.append('No')
        return output


    if vals[1] == 'R0' and vals[0] != 'FENE' and vals[0] != 'CRST': #FENE and CRST have different modulations or the radial part
        f1 = True

    elif vals[1] == 'RC' :
        f1 = True

    elif vals[1] == 'A':   
        f1 = True
   
    if f1:
       output.append('f1')
       return output
    
    
    #check if the parameter is in f2 modulation
    if len(vals) >= 2:
        if vals[1] == 'R0' and vals[0] == "CRST" :
            print("f2")
            f2 = True

    if f2 :        
        output.append('f2')
        
        return output
       

    #check if the parameter is in f4 modulation
    if len(vals) > 2 :
        if vals[2] == 'A' and (vals[1] != 'HYDR' or vals[1] != 'STCK') and  vals[0] != 'FENE':
            f4 = True

    if f4 :
       
        output.append('f4')
        
        return output
    
    if f1 == False and f2 == False and f4 == False:
        output.append('No')
        return output



def get_associated_pars(par) :
    
    all_fixed_pars.append(par)
    
    vals = par.split('_')
    name = vals[0]
    for k in range(1,len(vals)-2) :
        name = name + "_"+vals[k]
        
    if vals[1] == 'DELTA' :
        all_fixed_pars.append(name+"2"+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1])
        

    #symmetries
    if vals[0] == 'STCK':
        if vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'G' :
            all_fixed_pars.append(name+"_C_C")                      
        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'C' :
            all_fixed_pars.append(name+"_G_G")
            
        elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'A' :
            all_fixed_pars.append(name+"_T_C")
        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'C' :
            all_fixed_pars.append(name+"_G_A")
            
        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'G' :
            all_fixed_pars.append(name+"_C_T")
        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'T' :
            all_fixed_pars.append(name+"_A_G")
            
        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'G' :
            all_fixed_pars.append(name+"_C_A")
        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'A' :
            all_fixed_pars.append(name+"_T_G")
            
        elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'T' :
            all_fixed_pars.append(name+"_A_C")
        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'C' :
            all_fixed_pars.append(name+"_G_T")
        
        if symm_stck:
            if vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'A' :
                all_fixed_pars.append(name+"_T_T")
            elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'T' :
                all_fixed_pars.append(name+"_A_A")
        
    elif vals[0] == 'FENE':
        
        if vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'G' :
            all_fixed_pars.append(name+"_C_C")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_C_C")                    
        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'C' :
            all_fixed_pars.append(name+"_G_G")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_G_G")    
            
        elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'A' :
            all_fixed_pars.append(name+"_T_C")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_T_C")   
        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'C' :
            all_fixed_pars.append(name+"_G_A")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_G_A")    
            
        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'G' :
            all_fixed_pars.append(name+"_C_T")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_C_T")    
        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'T' :
            all_fixed_pars.append(name+"_A_G")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_A_G")  
            
        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'G' :
            all_fixed_pars.append(name+"_C_A")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_C_A")
        elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'A' :
            all_fixed_pars.append(name+"_T_G")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_T_G")   
            
        elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'T' :
            all_fixed_pars.append(name+"_A_C")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_A_C")  
        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'C' :
            all_fixed_pars.append(name+"_G_T")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_G_T")
            
        elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'A' :
            all_fixed_pars.append(name+"_T_T")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_T_T")    
        elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'T' :
            all_fixed_pars.append(name+"_A_A")
            if vals[1] == 'DELTA' :       
                all_fixed_pars.append(name+"2"+"_A_A")    
            
    elif vals[0] == 'CRST' or vals[0] == 'HYDR':
        if vals[len(vals)-2] !=  vals[len(vals)-1]  :
            all_fixed_pars.append(name+"_"+vals[len(vals)-1]+"_"+vals[len(vals)-2])
            
            
    #continuity
    output = impose_continuity(par)
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
        return
        
    #update seq dep file + impose symmetries
    for k in range(len(names)) :
        if vals[0] == 'STCK' :
            all_fixed_pars.append(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1])
            if vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'G' :
                all_fixed_pars.append(names[k]+"_C_C")
            elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'C' :
                all_fixed_pars.append(names[k]+"_G_G")
                
            elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'A' :
                all_fixed_pars.append(names[k]+"_T_C")
            elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'C' :
                all_fixed_pars.append(names[k]+"_G_A")
                
            elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'G' :
                all_fixed_pars.append(names[k]+"_C_T")
            elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'T' :
                all_fixed_pars.append(names[k]+"_A_G")
                
            elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'G' :
                all_fixed_pars.append(names[k]+"_C_A")
            elif vals[len(vals)-2] == 'C' and vals[len(vals)-1] == 'A' :
                all_fixed_pars.append(names[k]+"_T_G")
                
            elif vals[len(vals)-2] == 'G' and vals[len(vals)-1] == 'T' :
                all_fixed_pars.append(names[k]+"_A_C")
            elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'C' :
                all_fixed_pars.append(names[k]+"_G_T")
                
            elif vals[len(vals)-2] == 'A' and vals[len(vals)-1] == 'A' :
                all_fixed_pars.append(names[k]+"_T_T")
            elif vals[len(vals)-2] == 'T' and vals[len(vals)-1] == 'T' :
                all_fixed_pars.append(names[k]+"_A_A")
                
        elif vals[0] == 'CRST' or vals[0] == 'HYDR':
            all_fixed_pars.append(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1])
            if vals[len(vals)-2] !=  vals[len(vals)-1]  :
                all_fixed_pars.append(names[k]+"_"+vals[len(vals)-1]+"_"+vals[len(vals)-2])
        else :
            all_fixed_pars.append(names[k]+"_"+vals[len(vals)-2]+"_"+vals[len(vals)-1])
    
    
    
    return


if len(sys.argv) != 3 :
    print("Unknown argument format.")
    print("Usage: python3 get_opti_pars_from_config.py config_file par_file")
    sys.exit(1)
    
    
config = open(sys.argv[1],'r')

for line in config.readlines() :
    vals = line.split()
    
    if len(vals) == 0:
        continue
    elif vals[0][0] == '#':
        continue   
    
    #symm stacking default = True
    if(vals[0] == "SYMM_STCK") :
        if int(vals[1]) == 0 :
            symm_stck = True
        else :
            symm_stck = False
            
config.close()

config = open(sys.argv[1],'r')

for line in config.readlines() :
    vals = line.split()
    if len(vals) == 0:
        continue
    elif vals[0][0] == '#':
        continue  
        
    if(vals[0] == "cname") :
       get_associated_pars(vals[1])
       
config.close()

par_file = open(sys.argv[2],'r')
fixed_par_file = open("optimised_parameters.txt",'w')


for line in par_file.readlines() :
    vals = line.split()
    
    if len(vals) == 0:
        #print(line,file=fixed_par_file,end="")
        continue
    elif vals[0][0] == '#':
        #print(line,file=fixed_par_file,end="")
        continue
    
    if vals[0] in all_fixed_pars :
        print(line,file=fixed_par_file,end="")
    else:
        continue
     
par_file.close()
fixed_par_file.close()
    
    
