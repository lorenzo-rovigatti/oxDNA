#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 17:13:31 2023

@author: andrea


Includes functions to impose continuity constraints to the modulation coefficients.
There are also functions to plot the modulations.

"""

import math
import matplotlib.pyplot as plt


def Vmod(theta, a, theta0) :
    f = 1-a*(theta-theta0)**2
    return f

def Vsmooth(x,b,xc):
    f = b*(xc-x)**2
    return f

def Morse(r,epsilon,r0,a) :
    f = epsilon*(1-math.exp(-(r-r0)*a))*(1-math.exp(-(r-r0)*a))
    return f

def Vharmonic(x,x0):
    f = 0.5*(x0-x)**2
    return f

def f1(r,epsilon,r0,a,rc,bl,rcl,bh,rch,rl,rh) :
    f = 0
    if r >= rl and r <= rh :
        f = Morse(r,epsilon,r0,a) - Morse(rc,epsilon,r0,a)
    elif r > rcl and r < rl :
        f = epsilon*Vsmooth(r,bl,rcl)
    elif r > rh and r < rch :
        f = epsilon*Vsmooth(r,bh,rch)
    else :
        f = 0
        
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

def f4(theta, a, b, theta0, dtheta_s, dtheta_c) :
    f = 0
    if theta >= (theta0-dtheta_s) and theta <= (theta0+dtheta_s) :
        f = Vmod(theta,a,theta0)
    elif theta > (theta0-dtheta_c) and theta < (theta0-dtheta_s) :
        f = Vsmooth(theta,b,theta0-dtheta_c)
    elif theta > (theta0+dtheta_s) and theta < (theta0+dtheta_c) :
        f = Vsmooth(theta,b,theta0+dtheta_c)
    else :
        f = 0
    
    return f




#continuity constraints for f1


def plot_modulation_f1(r0,a,rc,bl,rcl,bh,rch,rl,rh,out_file) :
    erres = []
    mod = []
    epsilon = 1. #We set epsilon equal 1. for convenienece; epsilon diff from 1 just rescales the modulation.
    for t in range(-1200,1200) :
        r = r0+t/1000.
        erres.append(r)
        mod.append(f1(r,epsilon,r0,a,rc,bl,rcl,bh,rch,rl,rh))
        
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    # plt.title(r"All: $R=7$ $\phi=0.364$")
    ax.title.set_fontsize(20)
    ax.set_ylabel(r"$f1(r)$",fontsize=20)
    ax.set_xlabel(r"$r$ [oxDNA units]",fontsize=20)
    #ax.set_ylim(-0.1,1.2)
    ax.set_xlim(r0-1.2,r0+1.2)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.scatter(erres,mod,color="black")
    plt.savefig(out_file,bbox_inches='tight',pad_inches=0.05)
        


def check_continuity_f1(r0,a,rc,bl,rcl,bh,rch,rl,rh,eps) :
    continuous = True
    for t in range(-1200,1200) :
         r1 = r0+t/1000.
         r2 = r0+(t+1)/1000.
         epsilon = 1.
         d = abs( f1(r1,epsilon,r0,a,rc,bl,rcl,bh,rch,rl,rh)-f1(r2,epsilon,r0,a,rc,bl,rcl,bh,rch,rl,rh) ) 
         if d > eps :
             continuous = False
            
    return continuous

# given a value of the width a, compute other parameteres to impose continuity of f4 and its first derivative.
def continuity_f1(r0,a,rc) :
    
    rl = r0-0.08    #r low, the -0.08 is as in oxDNA2
    rh = r0+0.35      #r high, the +0.35 is as in oxDNA2
    
    while 1==1 :
    
        num = a*a*(math.exp(-(rl-r0)*a)-math.exp(-2*(rl-r0)*a))**2
        den = (1-math.exp(-(rl-r0)*a))**2-(1-math.exp(-(rc-r0)*a))**2
        
        bl = num/den
        
        num = a*a*(math.exp(-(rh-r0)*a)-math.exp(-2*(rh-r0)*a))**2
        den = (1-math.exp(-(rh-r0)*a))**2-(1-math.exp(-(rc-r0)*a))**2
        
        bh = num/den
        
        rcl = rl-a/bl*(math.exp(-(rl-r0)*a)-math.exp(-2*(rl-r0)*a))
        rch = rh-a/bh*(math.exp(-(rh-r0)*a)-math.exp(-2*(rh-r0)*a))
        
        if Vsmooth(rl,bl,rcl) < 0 and Vsmooth(rh,bh,rcl) < 0 :
            break
        
        if Vsmooth(rl,bl,rcl) > 0 :
            rl += 0.005
        if Vsmooth(rh,bh,rcl) > 0 :
            rh -= 0.01
    
    
    continuous = check_continuity_f1(r0,a,rc,bl,rcl,bh,rch,rl,rh,0.05)

    return rl,rh,bl,bh,rcl,rch,continuous



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


def check_continuity_f2(r0,rc,bl,rcl,bh,rch,rl,rh,eps) :
    continuous = True
    for t in range(-1200,1200) :
         r1 = r0+t/1000.
         r2 = r0+(t+1)/1000.
         d = 47.5*abs( f2(r1,r0,rc,bl,rcl,bh,rch,rl,rh)-f2(r2,r0,rc,bl,rcl,bh,rch,rl,rh) ) 
         if d > eps :
             continuous = False
             
    return continuous

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


def plot_modulation_f4(a,th0,dts,dtc,b,out_file) :
    angle = []
    mod = []
    for t in range(-10000,10000) :
        theta = math.pi*t/5000.
        angle.append(theta/math.pi)
        mod.append(f4(theta, a, b, th0, dts, dtc))
        
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    # plt.title(r"All: $R=7$ $\phi=0.364$")
    ax.title.set_fontsize(20)
    ax.set_ylabel(r"$f4(\theta)$",fontsize=20)
    ax.set_xlabel(r"$\theta$ [$\pi$]",fontsize=20)
    #ax.set_ylim(-0.1,1.2)
    ax.set_xlim(-1.2,1.2)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.scatter(angle,mod,color="black")
    plt.savefig(out_file,bbox_inches='tight',pad_inches=0.05)
        


def check_continuity_f4(a,th0,dts,dtc,b,eps) :
    continuous = True
    for t in range(-10000,10000) :
        theta1 = math.pi*t/5000.
        theta2 = math.pi*(t+1)/5000.        
        if abs( f4(theta1, a, b, th0, dts, dtc)-f4(theta2, a, b, th0, dts, dtc) ) > eps :
            continuous = False
    return continuous
    


# given a value of the width a, compute other parameteres to impose continuity of f4 and its first derivative.
def continuity_f4(a) :
    dts = math.sqrt(0.81225/a) #delta theta star; smoothing kicks in at fmod = 0.18775 (as in STCK oxdna2)
    dtc = 1./(a*dts) #delta theta c
    b = a*dts/(dtc-dts) #width b
    
    #continuous = check_continuity_f4(a,th0,dts,dtc,b,0.01)

    return dts,dtc,b




