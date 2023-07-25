#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:44:00 2023

@author: andrea

Use oxpy to compute the energy of a trajectory. 
Note: it requires the sequence dependent version of oxDNA2 or oxDNA3!

In this example, we compute the energy of a trajectory (obtained by running standared oxDNA2) with two models: oxDNA2 and oxDNA2 with modified stacking th5_t0 (see oxDNA2_custom_sequence_dependent_parameters.txt), and compare the results.
"""

#import numpy as np
import oxpy

#compute energy with same trajectory, but two different models (needed for reweighting!)
energy1 = []
energy2 = []
time = []

with oxpy.Context():
    
    #read input script specifying sequence dependent file 1 (i.e. model1)
    inp = oxpy.InputFile()
    inp.init_from_filename("input1.an")
    
    #create analysys backend
    backend = oxpy.analysis.AnalysisBackend(inp)

    obs = backend.config_info().observables
    
    while backend.read_next_configuration() :
        
        energy1.append(float(obs[0].get_output_string(backend.conf_step)))
        time.append(backend.conf_step)
        
        #print(obs2[0].get_output_string(backend1.conf_step))

        #print(ob.get_potential_energy())
        print("bk1: "+str(backend.conf_step))
        
with oxpy.Context():
    
    #read input script specifying sequence dependent file 2
    inp = oxpy.InputFile()
    inp.init_from_filename("input2.an")
    
    #create analysys backend
    backend = oxpy.analysis.AnalysisBackend(inp)

    obs = backend.config_info().observables
    
    while backend.read_next_configuration() :
        
        energy2.append(float(obs[0].get_output_string(backend.conf_step)))
        
        #print(obs2[0].get_output_string(backend1.conf_step))

        #print(ob.get_potential_energy())
        print("bk2: "+str(backend.conf_step))    
        
   
for i in range(len(time)) :
    print(str(time[i])+" "+str(energy1[i])+" "+str(energy2[i]))
        
        
