#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:50:09 2024

@author: yqb22156
"""

import config_multi as cg
import matplotlib.pyplot as plt



def unscrumble_gs(gs) :
    
    unsc_gs_all = []
   
        
    for j in range(len(cg.ids)) :
        
        unsc_gs = []
        
        for z in range(len(gs)) :
            if z%len(cg.ids) == j :
                unsc_gs.append(gs[z])
                
        unsc_gs_all.append(unsc_gs)
        
    return unsc_gs_all


def plot_gs_sampled() :
    
    unscr_gs_sampled = unscrumble_gs(cg.mu_sampled)
    unscr_gs_target = unscrumble_gs(cg.target_mu[cg.seq_id])
    
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
        plt.savefig("Seq"+str(cg.seq_id)+"_"+coord_name+".pdf",bbox_inches='tight',pad_inches=0.05)  
        
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


def plot_cov_sampled_diag() :
    
    unscr_cov_sampled = unscrumble_cov_diag(cg.cov_sampled)
    unscr_cov_target = unscrumble_cov_diag(cg.target_cov[cg.seq_id])
    
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
        plt.savefig("Seq"+str(cg.seq_id)+"_Cov_diag_"+coord_name+".pdf",bbox_inches='tight',pad_inches=0.05)  
        
        plt.close()        
        
    return
