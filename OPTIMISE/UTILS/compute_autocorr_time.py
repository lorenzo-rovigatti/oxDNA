#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 14:07:04 2022

@author: andrea bonato

code to map oxdna coordinates (center of mass + base vectors)
to internal coordinates (slide, twist, etc.)

it takes a oxdna trajectory and topology as input 
and outputs a file with the average internal coordinates


"""
import numpy as np
import math
from scipy import linalg
import sys
import os
import copy
#from oxdna_to_internal_triadII import read_oxdna_trajectory_standard_order #Carlon style (see 2017 J. chem. phys)
from oxdna_to_internal_cgna import read_oxdna_trajectory_standard_order #Maddocks style (see theses)

program_path=os.path.dirname(os.path.realpath(__file__))
current_path=os.getcwd()

print("Options for angles: deg or cay")
print("deg = degrees")
print("cay = radiants/5")

unit = 'cay'

alpha = 5

if len(sys.argv) != 5:
    print("Unknown argument format.")
    print("Usage: " + sys.argv[0] +" trajectory_file topology_file snap0 M")
    print("M is max tau for autocorrelation function (to avoid noise we don't sum over all trajectory when computing tau)")
    print("M >= C*tau(M), C ~ 5")
    sys.exit(1)

iname = sys.argv[1]
tname = sys.argv[2]
in_snap = int(sys.argv[3])
M = int(sys.argv[4])

print("Trajectory file: "+ iname)
print("Topology file: "+ tname)
print("Discarding first "+ str(in_snap-1)+ " snapshots.")

internal_coords = []

#given a junction trajectory (read_oxdna_trajectory_standard_order), store specific internal coordinates in
#global variable internal_coords
def store_internal_coord(traj,ids,in_j,fin_j,in_snap,overwrite=True) :

    global internal_coords
    coords = []

    Nsnaps = len(traj)
    Njuns = len(traj[0])

    alpha = 1. #5.*math.pi/180 #angles in cgna are in radiants/5

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
            internal_coords.append(coord)
        else :
            coords.append(coord)

    if overwrite == True :
        internal_coords = coords

    return

sigma2_sampled = None
mu_sampled = None

#compute averege internal coordinates
def ave_stored() :
    
    global mu_sampled
    global sigma2_sampled
    
    Nsnaps = len(internal_coords)
    Ncoords = len(internal_coords[0])
    
    mu_sampled = np.zeros(Ncoords)
    sigma2_sampled = np.zeros(Ncoords)
    
    for i in range(Ncoords) :
        mu_sampled[i] = 0.
        sigma2_sampled[i] = 0.
    
    #sequence dependent
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            mu_sampled[j] += internal_coords[i][j]/Nsnaps
            sigma2_sampled[j] += internal_coords[i][j]*internal_coords[i][j]/Nsnaps
            
    sigma2_sampled = sigma2_sampled - mu_sampled*mu_sampled  
    
    return


#read trajectory

#iname = 'trajectory.dat'
#tname = 'generated.top'

    
ifile = open(iname,'r')
tfile = open(tname,'r')

tr_data = read_oxdna_trajectory_standard_order(ifile, tfile)

Nj = len(tr_data[0])

ifile.close()
tfile.close()

ids = [7,8,11]

in_j = 3
fin_j = Nj - 3

store_internal_coord(tr_data,ids,in_j,fin_j,in_snap,True)
ave_stored()

print(mu_sampled)
print(sigma2_sampled)
    
Ncoords = len(internal_coords[0])
Nsnaps = len(internal_coords)
greater_0 = np.zeros(Ncoords)

tau = np.ones(Ncoords)

for i in range(1,M) :
    counts = 0
    ACorr = np.zeros(Ncoords)
    for j in range(Nsnaps-i-1) :
        counts += 1
        ACorr += (internal_coords[j]-mu_sampled)*(internal_coords[j+i+1]-mu_sampled)
    ACorr /= counts
    for j in range(len(ACorr)) : 
        if ACorr[j] < 0 :
            greater_0[j] = 1
    for j in range(len(tau)) :
        if greater_0[j] < 1:
            tau[j] += 2*ACorr[j]/sigma2_sampled[j]
    
print("ids:", ids)
print("Integrated autocorrelation time [# of snapshots]: ")
print(tau)
    
