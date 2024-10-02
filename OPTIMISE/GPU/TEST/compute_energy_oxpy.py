#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:11:37 2024

@author: yqb22156
"""

import sys
import oxpy


#compute energy of trajectory by using oxpy


inp = []
backend = []
obs = []

counts = 0

for i in range(100) :

	energy_sampled = []


	with oxpy.Context():
	    
		#read input script specifying sequence dependent file
		inp.append(oxpy.InputFile())

		inp[i].init_from_filename('input1.an')
		#create analysys backend
		backend.append(oxpy.analysis.AnalysisBackend(inp[i]))

		obs.append(backend[i].config_info().observables)

		read = False
		"""
		while backend[i].read_next_configuration() :

		counts+=1

		if(counts < in_snap) :
		    continue

		a = float(obs[i][0].get_output_string(backend[i].conf_step).split()[0])
		print(a)
		energy_sampled.append((Njuns+1)*10*a)
		"""
		
		counts = -1
		while 1==1 : 
			try:
				read =  backend[i].read_next_configuration()
			except:
				counts+=1
				energy_sampled.append(999)
				print("Warning: exception in oxpy energy computation; conf" +str(counts))
				continue
			if read == False :
				break
			counts+=1

			a = float(obs[i][0].get_output_string(backend[i].conf_step).split()[0])
			energy_sampled.append(48*10*a)
                
