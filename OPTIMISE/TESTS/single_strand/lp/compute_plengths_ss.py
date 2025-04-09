#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 12:46:13 2023

@author: yqb22156
"""


import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad

from oxdna_to_internal_triadII import read_oxdna_trajectory_standard_order #Carlon style (see 2017 J. chem. phys)
#from oxdna_to_internal_cgna import read_oxdna_trajectory_standard_order #Maddocks style (see theses)

#ids of internal coordinates we use to compute stiffness matrices. 9,10,11 = tilt, roll, twist
ids = [8,9,10,11]

Njuns = 99 #total number of junctions

in_j = 3 #ignore ends
fin_j = Njuns -3 #ignore ends
in_snap = 100 #ignore first in_snap snapshots (equilibration)

dimension = (Njuns-in_j-(Njuns-fin_j))*len(ids)

N_fit =  (Njuns-in_j-(Njuns-fin_j))

internal_coords = [] # internal coordinates

#average internal coordinates
mu_sampled = np.zeros(dimension, dtype = float)
mu_global_sampled = np.zeros(len(ids), dtype = float)
#cov_sampled = np.zeros((dimension,dimension), dtype = float)

e3s = []
#store normals
def store_normals(traj,in_j,fin_j,in_snap,overwrite=True) :

    global e3s

    e3s_loc = []

    Nsnaps = len(traj)
    Njuns = len(traj[0])

    for i in range(in_snap,Nsnaps) :
        e3 = []
        for j in range(Njuns) :
            if j < in_j or j >= fin_j :
                continue
            #print(traj[i][j].base_pair1.frame.orientation[:,2])
            e3.append(traj[i][j].base1.frame.orientation[:,2])

        if overwrite == False :
            e3s.append(e3)
        else :
            e3s_loc.append(e3)
           
    if overwrite == True :
        e3s = e3s_loc
                    
    return

#store chosen set of internal coordinates
def store_internal_coord(traj,ids,in_j,fin_j,in_snap,overwrite=True) :
    
    global internal_coords
     
    coords = []
    
    Nsnaps = len(traj)
    Njuns = len(traj[0])
    
    for i in range(in_snap,Nsnaps) :

        coord = []
        for j in range(Njuns) :
            #print("snap:"+str(i)+" junn:"+str(j))
            if j < in_j or j >= fin_j :
                continue
            for z in range(len(ids)) :
                #no intra for single strand
                """
                if ids[z] == 0 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[0])
                elif ids[z] == 1 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[1])
                elif ids[z] == 2 :
                    coord.append(traj[i][j].base_pair1.intra_coord.tran[2])
                elif ids[z] == 3 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[0])
                elif ids[z] == 4 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[1])
                elif ids[z] == 5 :
                    coord.append(traj[i][j].base_pair1.intra_coord.rot[2])
                """
                 
                elif ids[z] == 6 :
                    coord.append(traj[i][j].inter_coord.tran[0])
                elif ids[z] == 7 :
                    coord.append(traj[i][j].inter_coord.tran[1])
                elif ids[z] == 8 :
                    coord.append(traj[i][j].inter_coord.tran[2])
                elif ids[z] == 9 :
                    coord.append(traj[i][j].inter_coord.rot[0])
                elif ids[z] == 10 :
                    coord.append(traj[i][j].inter_coord.rot[1])
                elif ids[z] == 11 :
                    coord.append(traj[i][j].inter_coord.rot[2])
                    
        if overwrite == False :
            internal_coords.append(coord)
        else :
            coords.append(coord)
           
    if overwrite == True :
        internal_coords = coords
                    
    return

#compute averege internal coordinates
def ave_stored() :
    
    global mu_sampled
   # global cov_sampled
    
    Nsnaps = len(internal_coords)
    Ncoords = len(internal_coords[0])
    
    for i in range(Ncoords) :
        mu_sampled[i] = 0.
        #for j in range(Ncoords) :
        #    cov_sampled[i][j] = 0.
    
    #sequence dependent
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            mu_sampled[j] += internal_coords[i][j]/Nsnaps
            
    #everaged over sequence
    for i in range(Ncoords)        :
        mu_global_sampled[i%len(ids)]+=mu_sampled[i]/(Ncoords/len(ids))
    
    """
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            for z in range(j,Ncoords) :
                cov_sampled[j][z] += (internal_coords[i][j] - mu_sampled[j])*(internal_coords[i][z] - mu_sampled[z])/Nsnaps
    
    for j in range(Ncoords) :
        for z in range(j+1,Ncoords) :
            cov_sampled[z][j] = cov_sampled[j][z]    
    """    
    return


#Compute persistence length


#read oxdna trajectory and compute internal coordinates
iname = 'trajectory.dat'
tname = 'generated.top'

ifile = open(iname,'r')
tfile = open(tname,'r')


traj = read_oxdna_trajectory_standard_order(ifile, tfile)

ifile.close()
tfile.close()

# True = overwrite, False = append
# append is for averaging over multiple trajectories

store_internal_coord(traj,ids,in_j,fin_j,in_snap,True)
store_normals(traj,in_j,fin_j,in_snap,overwrite=True)
 
#compute average coordinates
ave_stored()


#compute stiffness matrices
N = int(len(internal_coords[0])/len(ids))

#print(mu_global_sampled[0])
av_rise = mu_global_sampled[0]*8.518/10.
print(av_rise)

def lb(m) :
    costheta = 0.
    counts = 0
    for i in range(len(e3s)):	#loop sampled confs
    	for j in range(N-m):
    	    #print(i,j)
    	    #print(np.dot(e3s[i][j],e3s[i][j+m]))
    	    costheta += np.dot(e3s[i][j],e3s[i][j+m])
    	    counts+=1
    #print(costheta)
    #print(math.log(costheta/counts))
    lb = -m*av_rise/math.log(costheta/counts)
    
    return lb
    
print(mu_global_sampled)
def lt(m) :
    cosomega3 = 0.
    counts = 0
    for i in range(len(internal_coords)):	#loop sampled confs
    	for j in range(N-m):
    	    omega3 = 0.
    	    for z in range(j,j+m):
    	        omega3 += (internal_coords[i][z*4+3]-mu_global_sampled[3])
    	    #print(omega3)
    	    cosomega3 += math.cos(omega3/180*math.pi)
    	    #print(i,j)
    	    #print(np.dot(e3s[i][j],e3s[i][j+m]))
    	    counts+=1
    #print(cosomega3)
    #print(math.log(cosomega3/counts))
    lt = -m*av_rise/math.log(cosomega3/counts)
    
    return lt
   
lb(11)

ofile = open("lps_test.txt","w")
for m in range(1,90) :
    print(m,lb(m),lt(m),file=ofile)
    
ofile.close()
