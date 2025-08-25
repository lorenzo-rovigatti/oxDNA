#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:06:08 2024

@author: yqb22156

plot internal coordinates and comparison with cgna+

"""
import sys
import numpy as np
import math
from oxdna_to_internal_cgna import read_oxdna_trajectory_standard_order
import get_cgdna_pars
import matplotlib.pyplot as plt
import os
import "../../colours.py"
path_style=os.path.dirname(os.path.realpath(__file__))

#style

plt.style.use("../../style.sty")

internal_coords = []
mu_sampled = []
cov_sampled = []
target_mu = []
target_cov = []


#given a junction trajectory (read_oxdna_trajectory_standard_order), store specific internal coordinates in
#global variable internal_coords
def store_internal_coord(traj,ids,in_j,jfe,in_snap,overwrite=True) :
       
    global internal_coords
    coords = []
    
    Nsnaps = len(traj)
    Njuns = len(traj[0])
    fin_j = Njuns-jfe-1
    print("Nsnap: " + str(Nsnaps))
    print("Njuns: " + str(Njuns))
    
    alpha = 5.*math.pi/180 #angles in cgna are in radiants/5
    
    print("js: "+str(in_j) + " "+ str(fin_j))
    print(ids)
    
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
            internal_coords.append(coord)
        else :
            coords.append(coord)
           
    if overwrite == True :
        internal_coords = coords
        
    
                    
    return

def ave_and_cov_stored() :
      
    global internal_coords
    global mu_sampled
    global cov_sampled

    Nsnaps = len(internal_coords)
    Ncoords = len(internal_coords[0])
    
    Nave = Nsnaps
    
    print("Computing ground state")
    
    for i in range(Ncoords) :
        mu_sampled.append(0.)
        cov_line = []
        for j in range(Ncoords) :
            cov_line.append(0)
        cov_sampled.append(cov_line)
        
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            mu_sampled[j] += internal_coords[i][j]/Nave
    
    print("Computing covariance")
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            for z in range(j,Ncoords) :
               cov_sampled[j][z] += (internal_coords[i][j] - mu_sampled[j])*(internal_coords[i][z] - mu_sampled[z])/Nave

    for j in range(Ncoords) :
        for z in range(j+1,Ncoords) :
            cov_sampled[z][j] = cov_sampled[j][z]
    
    return



def unscrumble_gs(gs) :
    
    unsc_gs_all = []
   
        
    for j in range(len(ids)) :
        
        unsc_gs = []
        
        for z in range(len(gs)) :
            if z%len(ids) == j :
                unsc_gs.append(gs[z])
                
        unsc_gs_all.append(unsc_gs)
        
    return unsc_gs_all


def plot_gs_sampled(ids, seq_id, seq) :
    
    global target_mu
    global mu_sampled
    global target_cov
    global cov_sampled
    
    unscr_gs_sampled = unscrumble_gs(mu_sampled)
    unscr_gs_target = unscrumble_gs(target_mu)
    
    for j in range(len(unscr_gs_sampled)) : #Nids
    
        coord_name = ""
        
        if ids[j] == 0:
            coord_name = "buckle"
        if ids[j] == 1:
            coord_name = "propeller"
        if ids[j] == 2:
            coord_name = "opening"
        if ids[j] == 3:
            coord_name = "shear"
        if ids[j] == 4:
            coord_name = "stretch"
        if ids[j] == 5:
            coord_name = "stagger"
        if ids[j] == 6:
            coord_name = "tilt"
        if ids[j] == 7:
            coord_name = "roll"
        if ids[j] == 8:
            coord_name = "twist"
        if ids[j] == 9:
            coord_name = "shift"
        if ids[j] == 10:
            coord_name = "slide"
        if ids[j] == 11:
            coord_name = "rise"  
    
        x = []            
        for z in range(len(unscr_gs_sampled[j])) :
            x.append(z)
            
        ys = unscr_gs_sampled[j]
        yt = unscr_gs_target[j]
        
        if (ids[j] == 0) or (ids[j] == 1) or(ids[j] == 2) or(ids[j] == 6) or (ids[j] == 7) or(ids[j] == 8):
            for i in range(len(unscr_gs_sampled[j])):
                unscr_gs_sampled[j][i]*=0.2*180/math.pi
                unscr_gs_target[j][i]*=0.2*180/math.pi
                
        ave = 0.
        counts = 0
        for i in range(len(unscr_gs_sampled[j])) : 
            ave += unscr_gs_sampled[j][i] 
            counts+=1
        ave /= counts
        print("SD ave GS "+coord_name+" = "+str(ave) )
        
        #print(ys)
        #print(yt)
        
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        # plt.title(r"All: $R=7$ $\phi=0.364$")
        ax.title.set_fontsize(20)
        ticks = []
        ticks_labels = []
        for i in range(len(unscr_gs_sampled[j])) :
            ticks.append(i)
            if ids[j] < 6: ticks_labels.append(seq[i])
            else: ticks_labels.append(seq[i]+"/"+seq[i+1])
        
        plt.xticks(ticks,ticks_labels)
        if ids[j] < 6: ax.set_xlabel(r"Sequence [base]",fontsize=20)
        else: ax.set_xlabel(r"Sequence [base step]",fontsize=20)
        if (ids[j] == 0) or (ids[j] == 1) or(ids[j] == 2) or(ids[j] == 6) or (ids[j] == 7) or(ids[j] == 8): ax.set_ylabel(coord_name + " [deg]",fontsize=20)
        else: ax.set_ylabel(coord_name + " [nm]",fontsize=20)
        #ax.set_ylim(0,160)
        #ax.set_xlim(-1.2,1.2)
        if ids[j] < 6: ax.tick_params(axis='both', which='major', labelsize=20)
        else:
            plt.xticks(rotation=60)
            ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        ax.plot(x, yt, color=CCS_G[1], label="cgna+")
        ax.plot(x, ys, color=CCS_G[7], label="oxdna")
        ax.legend(fontsize = 20)
        plt.savefig(coord_name+"_"+str(seq_id)+".pdf",bbox_inches='tight',pad_inches=0.05)
        
        plt.close()        
        
    return


def unscrumble_cov_diag(cov) :
    
    unsc_cov_all = []
   
        
    for j in range(len(ids)) :
        
        unsc_cov = []
        
        for z in range(len(cov)) :
            if z%len(ids) == j :
                unsc_cov.append(cov[z][z])
                
        unsc_cov_all.append(unsc_cov)
        
    return unsc_cov_all


def plot_cov_sampled_diag(ids, seq_id, seq) :
    
    global target_mu
    global mu_sampled
    global target_cov
    global cov_sampled
    
    unscr_cov_sampled = unscrumble_cov_diag(cov_sampled)
    unscr_cov_target = unscrumble_cov_diag(target_cov)
    
    for j in range(len(unscr_cov_sampled)) : #Nids
    
        coord_name = ""
        
        if ids[j] == 0:
            coord_name = "buckle"
        if ids[j] == 1:
            coord_name = "propeller"
        if ids[j] == 2:
            coord_name = "opening"
        if ids[j] == 3:
            coord_name = "shear"
        if ids[j] == 4:
            coord_name = "stretch"
        if ids[j] == 5:
            coord_name = "stagger"
        if ids[j] == 6:
            coord_name = "tilt"
        if ids[j] == 7:
            coord_name = "roll"
        if ids[j] == 8:
            coord_name = "twist"
        if ids[j] == 9:
            coord_name = "shift"
        if ids[j] == 10:
            coord_name = "slide"
        if ids[j] == 11:
            coord_name = "rise"  
    
        x = []            
        for z in range(len(unscr_cov_sampled[j])) :
            x.append(z)
            
        ys = unscr_cov_sampled[j]
        yt = unscr_cov_target[j]
        
        if (ids[j] == 0) or (ids[j] == 1) or(ids[j] == 2) or(ids[j] == 6) or (ids[j] == 7) or(ids[j] == 8):
            for i in range(len(unscr_cov_sampled[j])):
                unscr_cov_sampled[j][i]=math.sqrt(unscr_cov_sampled[j][i])*0.2*180/math.pi
                unscr_cov_target[j][i]=math.sqrt(unscr_cov_target[j][i])*0.2*180/math.pi
                
                
        ave = 0.
        counts = 0
        for i in range(len(unscr_cov_sampled[j])) : 
            ave += unscr_cov_sampled[j][i] 
            counts+=1
        ave /= counts
        print("SD ave STD "+coord_name+" = "+str(ave) )
        
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        # plt.title(r"All: $R=7$ $\phi=0.364$")
        ax.title.set_fontsize(20)
        ticks = []
        ticks_labels = []
        for i in range(len(unscr_cov_sampled[j])) :
            ticks.append(i)
            if ids[j] < 6: ticks_labels.append(seq[i])
            else: ticks_labels.append(seq[i]+"/"+seq[i+1])
        
        plt.xticks(ticks,ticks_labels)
        if ids[j] < 6: ax.set_xlabel(r"Sequence [base]",fontsize=20)
        else: ax.set_xlabel(r"Sequence [base step]",fontsize=20)
        if (ids[j] == 0) or (ids[j] == 1) or(ids[j] == 2) or(ids[j] == 6) or (ids[j] == 7) or(ids[j] == 8): ax.set_ylabel("Std "+ coord_name + " [deg]",fontsize=20)
        else: ax.set_ylabel("Std "+ coord_name + " [nm]",fontsize=20)
        #ax.set_ylim(0,160)
        #ax.set_xlim(-1.2,1.2)
        if ids[j] < 6: ax.tick_params(axis='both', which='major', labelsize=20)
        else:
            plt.xticks(rotation=60)
            ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        ax.plot(x, yt, color=CCS_G[1], label="cgna+")
        ax.plot(x, ys, color=CCS_G[7], label="oxdna")
        ax.legend(fontsize = 20)
        plt.savefig("Std_"+coord_name+"_"+str(seq_id)+".pdf",bbox_inches='tight',pad_inches=0.05)  
        
        plt.close()        
        
    return


###################################################
############### MAIN #################################
###################################################

if len(sys.argv) != 4 :
    print("Unknown argument format.")
    print("Usage: python3 "+sys.argv[0]+" trj_file topo_file seq_id")
    sys.exit(1)


#iname = './Seq'+str(l)+'/Rep'+str(i)+'/trajectory.dat'
#tname = './Seq'+str(l)+'/Rep'+str(i)+'/generated.top'

#read trajectory and compute internal coordinates

ifile = open(sys.argv[1],'r')
tfile = open(sys.argv[2],'r')
seq_id = int(sys.argv[3])

traj = read_oxdna_trajectory_standard_order(ifile, tfile)

ifile.close()
tfile.close()


seq = ""

#read sequence from topo file
counts = -1

tfile = open(sys.argv[2],'r')

for line in tfile.readlines() :
    
    counts += 1
    vals = line.split(" ")
    if counts == 0:
        continue
    if vals[0] == "1":
        seq += vals[1]
    
tfile.close()

print("Sequence: "+seq)

ids = [1,6,7,8,11]

in_snap = 20
in_j = 3
jfe = 3

#store selected internal coordinates and compute gs and cov

store_internal_coord(traj,ids,in_j,jfe,in_snap,overwrite=True)
print("len: " + str(len(internal_coords[0])))
ave_and_cov_stored()


#compute gs and cov from cgna+

target_mu, target_cov = get_cgdna_pars.get_target_mean_and_covariance(seq, ids, in_j, jfe)
target_mu, target_cov = get_cgdna_pars.get_target_mean_and_covariance(seq, ids, in_j, jfe)

seq_stripped = ""
for i in range(in_j+1,len(seq)) : seq_stripped += seq[i]

#plot
plot_gs_sampled(ids,seq_id,seq_stripped)
plot_cov_sampled_diag(ids,seq_id,seq_stripped)
