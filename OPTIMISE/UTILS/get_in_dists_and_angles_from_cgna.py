#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:11:15 2024

@author: yqb22156
"""

dst = 0.26  #stacking dist from center
dh = 0.2
dbb = 1.0

import sys
import math
import matplotlib.pyplot as plt
from get_cgdna_pars import get_target_mean_and_covariance
import numpy as np

ox_to_Ang = 8.518

def Vsmooth(x,b,xc):
    f = b*(xc-x)**2
    return f

def Vharmonic(x,x0):
    f = 0.5*(x0-x)**2
    return f

def f2(r,r0,rc,bl,rcl,bh,rch,rl,rh) :
    f = 0
    if r >= rl and r <= rh :
        f = Vharmonic(r,r0)-Vharmonic(rc,r0)
    elif r > rcl and r < rl :
        f = Vsmooth(r,bl,rcl)
    elif r > rh and r < rch :
        f = Vsmooth(r,bh,rch)
    else :
        f = 0
        
    return f

def check_continuity_f2(r0,rc,bl,rcl,bh,rch,rl,rh,eps) :
    continuous = True
    for t in range(-1200,1200) :
         r1 = r0+t/1000.
         r2 = r0+(t+1)/1000.
         d = 47.5*abs( f2(r1,r0,rc,bl,rcl,bh,rch,rl,rh)-f2(r2,r0,rc,bl,rcl,bh,rch,rl,rh) ) 
         if d > eps :
             continuous = False
             
    return continuous

def plot_modulation_f2(r0,rc,bl,rcl,bh,rch,rl,rh,out_file) :
    erres = []
    mod = []
    
    k = 47.5/abs(rc-r0)/abs(rc-r0)*0.1*0.1 #rescaled base value in ox2
    
    for t in range(-1200,1200) :
        r = r0+t/1000.
        erres.append(r)
        mod.append(k*f2(r,r0,rc,bl,rcl,bh,rch,rl,rh))
        
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    # plt.title(r"All: $R=7$ $\phi=0.364$")
    ax.title.set_fontsize(20)
    ax.set_ylabel(r"$f2(r)$",fontsize=20)
    ax.set_xlabel(r"$r$ [oxDNA units]",fontsize=20)
    #ax.set_ylim(-0.1,1.2)
    ax.set_xlim(r0-1.2,r0+1.2)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.scatter(erres,mod,color="black")
    plt.savefig(out_file,bbox_inches='tight',pad_inches=0.05)



#The amplitude of the modulation f2 is regulated by abs(r0-rc). This difference also affects the depth.
#to control depth and amplitude separately, we have to rescale the parameter k by multiplying it by 0.1*0.1/abs(r0-rc)/abs(r0-rc).
#Here we fix the amplitude by setting rc=r0+0.1, so that abs(r0-rc) = 0.1 always. In this way the depth is directly regulated by k.
def continuity_f2_fixed_width(r0) :
    
    #note: f2, contrary to f1, is symmetric
    rl = r0-0.08    #r low, the -0.08 is as in oxDNA2
    rh = r0+0.08      #r high, the +0.08 is as in oxDNA2
    rc = r0+0.1     #as in ox2
    
    while 1==1 :
        
        
        rcl = rl - 1./(rl-r0)*( (rl-r0)*(rl-r0)-(rc-r0)*(rc-r0) )
        rch = rh - 1./(rh-r0)*( (rh-r0)*(rh-r0)-(rc-r0)*(rc-r0) )
        
        bl = 0.5*(r0-rl)/(rcl-rl)
        bh = 0.5*(r0-rh)/(rch-rh)
        
        if Vsmooth(rl,bl,rcl) < 0 and Vsmooth(rh,bh,rcl) < 0 :
            break
        
        #note: f2, contrary to f1, is symmetric
        if Vsmooth(rl,bl,rcl) > 0 :
            rl += 0.005
        if Vsmooth(rh,bh,rcl) > 0 :
            rh -= 0.005
    
    
    continuous = check_continuity_f2(r0,rc,bl,rcl,bh,rch,rl,rh,0.05)

    return rc,rl,rh,bl,bh,rcl,rch,continuous

#continuity constraints for f4

def get_pars(rise,twist):    #twist in rad*5, rise in Ang [cgdna units]
    twist = twist/5.
    rise = rise/ox_to_Ang
    r0 = math.sqrt(rise*rise+2*dst*dst*(1-math.cos(twist)))
    r0 = r0 -0.03 #oxdna correction   
    R0 = math.sqrt(rise*rise+2*dbb*dbb*(1-math.cos(twist)))
    R0 = R0 + 0.006
    
    crst_r0 = math.sqrt(rise*rise+2*dh*dh*(1+math.cos(twist)))
    th0_1 = twist
    th0_1 = th0_1 + 0.14 #oxdna correction    
    th0_2 = math.acos(dh/crst_r0*(1+math.cos(twist)))
    th0_2 = th0_2 + 0.14 #oxdna correction 
    th0_7 = math.acos(rise/crst_r0)
    th0_7 = th0_7 + 0.14 #oxdna correction 
    crst_r0 = crst_r0 + 0.015 #oxdna correction 
    
    return r0,R0,crst_r0,th0_1,th0_2,th0_7


bases = ['A','C','G','T']

def base_to_id(b):
    if b == 'A':
        return 0
    elif b == 'C':
        return 1
    elif b == 'G':
        return 2
    elif b == 'T':
        return 3
    else:
        return -1


"""
rise = 3.41
twist = 0.646*5

r0,R0,crst_r0,th0_1,th0_2,th0_7=get_pars(rise,twist)

rc,rl,rh,bl,bh,rcl,rch,continuous = continuity_f2_fixed_width(crst_r0)

#out_file = "f2.pdf"

#plot_modulation_f2(crst_r0,rc,bl,rcl,bh,rch,rl,rh,out_file) 
"""

"""
if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 get_in_dists....py sequences_file")
    sys.exit(1)

ifile = open(sys.argv[1],'r')
"""

ifile = open('seq.txt','r')


seqs = []

rise = np.zeros((4,4,4,4),dtype=float )
twist = np.zeros((4,4,4,4),dtype=float )
counts = np.zeros((4,4,4,4),dtype=float )

            

for line in ifile.readlines() :
    seqs.append(line)
    
ifile.close()
ids = [8,11]
jx_in = 2
jx_from_end = 2
    
for seq in seqs:
    
    red_gs, red_cov = get_target_mean_and_covariance(seq, ids, jx_in, jx_from_end)
    cnt = 0
    for i in range(jx_in, len(seq)-jx_from_end):
        print(seq[i])
        twist[base_to_id(seq[i-1])][base_to_id(seq[i])][base_to_id(seq[i+1])][base_to_id(seq[i+1])] += red_gs[(i-jx_in)*2]
        rise[base_to_id(seq[i-1])][base_to_id(seq[i])][base_to_id(seq[i+1])][base_to_id(seq[i+1])] += red_gs[(i-jx_in)*2+1]
        counts[base_to_id(seq[i-1])][base_to_id(seq[i])][base_to_id(seq[i+1])][base_to_id(seq[i+1])]+=1


for i in range(4):
    for j in range(4):
        for l in range(4):
            for m in range(4):                
                if counts[i][j][l][m] == 0:
                    print("warning. untrained tetramer")
                    print(i,j,l,m)
                    rise[i][j][l][m] = 3.35
                    twist[i][j][l][m] = 3.2
                else :
                    rise[i][j][l][m] /= counts[i][j][l][m] 
                    twist[i][j][l][m] /= counts[i][j][l][m] 
                    
                    print(bases[i]+bases[j]+[bases[l]+bases[m]]+" "+str(rise[i][j][l][m])+" "+str(twist[i][j][l][m]))
                   
                    

