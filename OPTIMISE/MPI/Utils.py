#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 14:50:09 2024

@author: yqb22156
"""

import config_multi as cg
import matplotlib.pyplot as plt



def unscrumble_gs(gs) :
    
    unscrumbled_gs_all_seqs = []
    
    for i in range(cg.Nseq) :
        
        unscrumbled_gs = []        
        
        for j in range(len(cg.ids)) :
            
            unscr_gs_single = []
            
            for z in range(len(gs[i])) :
                if z%len(cg.ids) == j :
                    unscr_gs_single.append(gs[i][z])
                    
            unscrumbled_gs.append(unscr_gs_single)
            
        unscrumbled_gs_all_seqs.append(unscrumbled_gs)
        
    return unscrumbled_gs_all_seqs


def plot_gs_sampled() :
    
    unscr_gs_sampled = unscrumble_gs(cg.mu_sampled)
    unscr_gs_target = unscrumble_gs(cg.target_mu)
    
    for i in range(len(unscr_gs_sampled)) :    #NSeq
    
        for j in range(len(unscr_gs_sampled[i])) : #Nids
        
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
            for z in range(len(unscr_gs_sampled[i][j])) :
                x.append(z)
                
            ys = unscr_gs_sampled[i][j]
            yt = unscr_gs_target[i][j]
                
            
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
            plt.savefig("Seq"+str(i)+"_"+coord_name+".pdf",bbox_inches='tight',pad_inches=0.05)  
            
            plt.close()

        
        
    return