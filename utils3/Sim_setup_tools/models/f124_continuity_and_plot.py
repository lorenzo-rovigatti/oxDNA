#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 17:13:31 2023

@author: andrea
"""

import numpy as np
import math
from scipy import linalg
import sys
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

def Harmonic(r,k,r0) :
    f = k/2.*(r-r0)*(r-r0)
    return f

#f1
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

#f2
def f2(r,k,r0,rc,bl,rcl,bh,rch,rl,rh) :
    f = 0
    if r >= rl and r <= rh :
        f = Harmonic(r,k,r0) - Harmonic(rc,k,r0)
    elif r > rcl and r < rl :
        f = k*Vsmooth(r,bl,rcl)
    elif r > rh and r < rch :
        f = k*Vsmooth(r,bh,rch)
    else :
        f = 0
        
    return f

#f4
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

#function to impose continuity of f4 and its derivative (see report). Delta theta* (dts) is set such that
#smoothing kicks in at fmod = 0.18775 (as in oxdna2 for theta5,6 stacking)
def continuity_f4_v1(a):
    
    dts = math.sqrt(0.81225/a)
    dtc = 1./(a*dts)
    b = a*dts/(dtc-dts)
    return b,dts,dtc

#function to impose continuity of f4 and its first derivative. delta theta* (dts) is choosen by the user
def continuity_f4_v2(dts,a):
    
    dtc = 1./(a*dts)
    b = a*dts/(dtc-dts)
    return b,dts,dtc


#function to impose continuity of f1 and its first derivative.
def continuity_f1(rl,r0,a,rc) :
    num = a*a*(math.exp(-(rl-r0)*a)-math.exp(-2*(rl-r0)*a))**2
    den = (1-math.exp(-(rl-r0)*a))**2-(1-math.exp(-(rc-r0)*a))**2
    b = num/den
    rcl = rl-a/b*(math.exp(-(rl-r0)*a)-math.exp(-2*(rl-r0)*a))
    
    return b,rcl


#Impose continuity of f2 and its first derivative.
def continuity_f2(rl,r0,rc) :
    rcl = rl-((rl-r0)*(rl-r0)-(rc-r0)*(rc-r0))/(rl-r0)
    b = -0.5*(rl-r0)/(rcl-rl)
    return b,rcl


###########################
#continuity and plot: f4

#OX2 values
STCK_THETA5_A = 0.9
STCK_THETA5_B = 3.89361
STCK_THETA5_T0 = 0.
STCK_THETA5_TS = 0.95
STCK_THETA5_TC = 1.16959
STCK_THETA6_A = 0.9
STCK_THETA6_B = 3.89361
STCK_THETA6_T0 = 0.
STCK_THETA6_TS = 0.95
STCK_THETA6_TC = 1.16959

# given a value of STCK_THETA5_A, compute other parameteres to impose continuity of f4 and its first derivative.

print("############################\n")
print("f4, first run")


print("th0 = ")
print(str(STCK_THETA5_T0)+"\n")
print("a = ")
print(str(STCK_THETA5_A)+"\n")

b,dts,dtc = continuity_f4_v1(STCK_THETA5_A)

print("dtheta_s =")
print(str(dts)+"\n")
print("dtheta_c =")
print(str(dtc)+"\n")
print("b =")
print(str(b)+"\n")

STCK_THETA5_TS = dts
STCK_THETA5_TC =dtc
STCK_THETA5_B = b


#save

V1 = []
VM = []
VS1 = []
VS2 = []
angle = []

for t in range(-5000,5000) :
    theta = math.pi*t/5000.
    angle.append(theta/math.pi)
    V1.append(f4(theta, STCK_THETA5_A, STCK_THETA5_B,STCK_THETA5_T0,STCK_THETA5_TS,STCK_THETA5_TC))
    VM.append(Vmod(theta,STCK_THETA5_A,STCK_THETA5_T0))
    VS1.append(Vsmooth(theta,STCK_THETA5_B,STCK_THETA5_T0-STCK_THETA5_TC))
    VS2.append(Vsmooth(theta,STCK_THETA5_B,STCK_THETA5_T0+STCK_THETA5_TC))

###################################
#second run

print("############################\n")
print("f4, second run")

STCK_THETA5_T0 = 0 #rad
STCK_THETA5_A = 0.8
print("th0 = ")
print(str(STCK_THETA5_T0)+"\n")
print("a = ")
print(str(STCK_THETA5_A)+"\n")

b,dts,dtc = continuity_f4_v1(STCK_THETA5_A)

print("dtheta_s =")
print(str(dts)+"\n")
print("dtheta_c =")
print(str(dtc)+"\n")
print("b =")
print(str(b)+"\n")

STCK_THETA5_TS = dts
STCK_THETA5_TC =dtc
STCK_THETA5_B = b


V2 = []    
for t in range(-5000,5000) :
     theta = math.pi*t/5000.
     V2.append(f4(theta, STCK_THETA5_A, STCK_THETA5_B,STCK_THETA5_T0,STCK_THETA5_TS,STCK_THETA5_TC))   
    
    
###################################   
#third run

print("############################\n")
print("f4, third run")

STCK_THETA5_T0 = 0.0 #rad
STCK_THETA5_A = 0.7
print("th0 = ")
print(str(STCK_THETA5_T0)+"\n")
print("a = ")
print(str(STCK_THETA5_A)+"\n")


b,dts,dtc = continuity_f4_v1(STCK_THETA5_A)

print("dtheta_s =")
print(str(dts)+"\n")
print("dtheta_c =")
print(str(dtc)+"\n")
print("b =")
print(str(b)+"\n")

STCK_THETA5_TS = dts
STCK_THETA5_TC =dtc
STCK_THETA5_B = b

  
V3 = []    
for t in range(-5000,5000) :
     theta = math.pi*t/5000.
     V3.append(f4(theta, STCK_THETA5_A, STCK_THETA5_B,STCK_THETA5_T0,STCK_THETA5_TS,STCK_THETA5_TC))   
    
###################################     
#forth run
print("f4, forth run")

STCK_THETA5_T0 = 0 #rad
STCK_THETA5_A = 0.3
print("th0 = ")
print(str(STCK_THETA5_T0)+"\n")
print("a = ")
print(str(STCK_THETA5_A)+"\n")

b,dts,dtc = continuity_f4_v1(STCK_THETA5_A)

print("dtheta_s =")
print(str(dts)+"\n")
print("dtheta_c =")
print(str(dtc)+"\n")
print("b =")
print(str(b)+"\n")

STCK_THETA5_TS = dts
STCK_THETA5_TC =dtc
STCK_THETA5_B = b
  
V4 = []    
for t in range(-5000,5000) :
     theta = math.pi*t/5000.
     V4.append(f4(theta, STCK_THETA5_A, STCK_THETA5_B,STCK_THETA5_T0,STCK_THETA5_TS,STCK_THETA5_TC))    

     
###############################
#plot


fig = plt.figure(figsize=(5,5))
ax =fig.add_axes([0., 0., 1.0, 1.0])
#plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)


ax.set_ylabel(r"$f4(\theta)$",fontsize=20)
ax.set_xlabel(r"$\theta \, [\pi]$",fontsize=20)
ax.set_ylim(-0.1,1.1)
ax.set_xlim(-1,1)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)

ax.scatter(angle,V1,color="black",label=r"$a=0.9$")
ax.scatter(angle,V2,color="red",label=r"$a=0.8$")
ax.scatter(angle,V3,color="blue",label=r"$a=0.7$")
ax.scatter(angle,V4,color="silver",label=r"$a=0.3$")
#ax.scatter(angle,V1,label="f4(theta5)")
#ax.scatter(angle,VM,label="fmod(theta5)")
#ax.scatter(angle,VS1,label="fsm1(theta5)")
#ax.plot(angle,VS1,label=r"smooth")
#ax.scatter(angle,VS2,label="fsm2(theta5)")

ax.legend(bbox_to_anchor=(0.75, 1.05), prop={'size':24})
     
     
#####################################
#continuity and plot: f1

print("############################\n")
print("f1, first run")

#ox2 values
    
STCK_F1 = 1
STCK_BASE_EPS_OXDNA = 1.3448
STCK_BASE_EPS_OXDNA2 = 1.3523
STCK_FACT_EPS_OXDNA = 2.6568
STCK_FACT_EPS_OXDNA2 = 2.6717
STCK_A = 6.
STCK_RC = 0.9
STCK_R0 = 0.40
STCK_BLOW = -68.1857
STCK_BHIGH = -3.12992
STCK_RLOW = 0.32
STCK_RHIGH = 0.75
STCK_RCLOW = 0.23239
STCK_RCHIGH = 0.956




#continuity for f1
    
bl,rcl = continuity_f1(STCK_RLOW,STCK_R0,STCK_A,STCK_RC)
#run1
print("b low: ")
print(bl)
print("rc low: ")
print(rcl)

STCK_BLOW = bl
STCK_RCLOW = rcl

bh,rch = continuity_f1(STCK_RHIGH,STCK_R0,STCK_A,STCK_RC)

print("b high: ")
print(bh)
print("rc high: ")
print(rch)

STCK_BHIGH = bh
STCK_RCHIGH = rch


V_f1_1 = []
erre = []
for t in range(0,2200) :
     r = t/2000.
     erre.append(r)
     V_f1_1.append(f1(r, STCK_BASE_EPS_OXDNA2, STCK_R0,STCK_A,STCK_RC,STCK_BLOW,STCK_RCLOW,STCK_BHIGH,STCK_RCHIGH,STCK_RLOW,STCK_RHIGH))    
 
####################################
#print minimum! 
MIN = - Morse(STCK_RC,STCK_BASE_EPS_OXDNA2,STCK_R0,STCK_A) 
print("Minimum: ")
print(MIN)
 
#####################################
# run2

print("############################\n")
print("f1, second run")

STCK_A = 4.
STCK_RC = 0.9
STCK_R0 = 0.40
STCK_RLOW = 0.32
STCK_RHIGH = 0.75


bl,rcl = continuity_f1(STCK_RLOW,STCK_R0,STCK_A,STCK_RC)

print("b low: ")
print(bl)
print("rc low: ")
print(rcl)

STCK_BLOW = bl
STCK_RCLOW = rcl

bh,rch = continuity_f1(STCK_RHIGH,STCK_R0,STCK_A,STCK_RC)

print("b high: ")
print(bh)
print("rc high: ")
print(rch)

STCK_BHIGH = bh
STCK_RCHIGH = rch

 
V_f1_2 = []
for t in range(0,2200) :
     r = t/2000.
     V_f1_2.append(f1(r, STCK_BASE_EPS_OXDNA2, STCK_R0,STCK_A,STCK_RC,STCK_BLOW,STCK_RCLOW,STCK_BHIGH,STCK_RCHIGH,STCK_RLOW,STCK_RHIGH))  
     
MIN = - Morse(STCK_RC,STCK_BASE_EPS_OXDNA2,STCK_R0,STCK_A) 
print("Minimum: ")
print(MIN)
 
 
#####################################
# run3

print("############################\n")
print("f1, third run")

STCK_A = 5.
STCK_RC = 0.9
STCK_R0 = 0.40
STCK_RLOW = 0.32
STCK_RHIGH = 0.75

bl,rcl = continuity_f1(STCK_RLOW,STCK_R0,STCK_A,STCK_RC)

print("b low: ")
print(bl)
print("rc low: ")
print(rcl)

STCK_BLOW = bl
STCK_RCLOW = rcl

bh,rch = continuity_f1(STCK_RHIGH,STCK_R0,STCK_A,STCK_RC)

print("b high: ")
print(bh)
print("rc high: ")
print(rch)

STCK_BHIGH = bh
STCK_RCHIGH = rch
 
V_f1_3 = []
for t in range(0,2200) :
     r = t/2000.
     V_f1_3.append(f1(r, STCK_BASE_EPS_OXDNA2, STCK_R0,STCK_A,STCK_RC,STCK_BLOW,STCK_RCLOW,STCK_BHIGH,STCK_RCHIGH,STCK_RLOW,STCK_RHIGH))     

MIN = - Morse(STCK_RC,STCK_BASE_EPS_OXDNA2,STCK_R0,STCK_A) 
print("Minimum: ")
print(MIN)
 
#####################################
# run4

print("############################\n")
print("f1, forth run")

STCK_A = 7.
STCK_RC = 0.9
STCK_R0 = 0.4
STCK_RLOW = 0.32
STCK_RHIGH = 0.75

bl,rcl = continuity_f1(STCK_RLOW,STCK_R0,STCK_A,STCK_RC)

print("b low: ")
print(bl)
print("rc low: ")
print(rcl)

STCK_BLOW = bl
STCK_RCLOW = rcl

bh,rch = continuity_f1(STCK_RHIGH,STCK_R0,STCK_A,STCK_RC)

print("b high: ")
print(bh)
print("rc high: ")
print(rch)

STCK_BHIGH = bh
STCK_RCHIGH = rch
 
V_f1_4 = []
for t in range(0,2200) :
     r = t/2000.
     V_f1_4.append(f1(r, STCK_BASE_EPS_OXDNA2, STCK_R0,STCK_A,STCK_RC,STCK_BLOW,STCK_RCLOW,STCK_BHIGH,STCK_RCHIGH,STCK_RLOW,STCK_RHIGH)) 
     
     
MIN = - Morse(STCK_RC,STCK_BASE_EPS_OXDNA2,STCK_R0,STCK_A) 
print("Minimum: ")
print(MIN)
 

#####################################
# run5

print("############################\n")
print("f1, fifth run")

STCK_A = 8.
STCK_RC = 0.9
STCK_R0 = 0.44
STCK_RLOW = 0.36
STCK_RHIGH = 0.75

bl,rcl = continuity_f1(STCK_RLOW,STCK_R0,STCK_A,STCK_RC)

print("b low: ")
print(bl)
print("rc low: ")
print(rcl)

STCK_BLOW = bl
STCK_RCLOW = rcl

bh,rch = continuity_f1(STCK_RHIGH,STCK_R0,STCK_A,STCK_RC)

print("b high: ")
print(bh)
print("rc high: ")
print(rch)

STCK_BHIGH = bh
STCK_RCHIGH = rch

 
V_f1_5 = []
for t in range(0,2200) :
     r = t/2000.
     V_f1_5.append(f1(r, STCK_BASE_EPS_OXDNA2, STCK_R0,STCK_A,STCK_RC,STCK_BLOW,STCK_RCLOW,STCK_BHIGH,STCK_RCHIGH,STCK_RLOW,STCK_RHIGH))  
     
     
MIN = - Morse(STCK_RC,STCK_BASE_EPS_OXDNA2,STCK_R0,STCK_A) 
print("Minimum: ")
print(MIN)
 
     
     
######################################
# plot 


fig = plt.figure(figsize=(5,5))
ax =fig.add_axes([0., 0., 1.0, 1.0])
#plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)

ax.set_ylabel(r"$f1(r)$",fontsize=20)
ax.set_xlabel(r"$r$ [oxDNA units]",fontsize=20)
#ax.set_ylim(-0.1,1.2)
ax.set_xlim(-0.1,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
    
ax.scatter(erre,V_f1_2,color="black",label=r"$a=4$")
ax.scatter(erre,V_f1_3,color="red",label=r"$a=5$")
ax.scatter(erre,V_f1_1,color="silver",label=r"$a=6$")
ax.scatter(erre,V_f1_4,color="blue",label=r"$a=7$")
#ax.plot(erre,V_f1_5,color="magenta",label=r"$a=8$")
      
     
ax.legend(bbox_to_anchor=(0.75, 0.7), prop={'size':24})
   

#####################################
#continuity for f2


CRST_R0 = 0.575
CRST_RC = 0.675
CRST_K = 47.5
CRST_BLOW = -0.888889
CRST_RLOW = 0.495
CRST_RCLOW = 0.45
CRST_BHIGH = -0.888889
CRST_RHIGH = 0.655
CRST_RCHIGH = 0.7
    
####################################
#run1

print("############################\n")
print("f2, first run")

CRST_R0 = 0.61
CRST_RC = 0.71
CRST_RLOW = 0.525
CRST_RHIGH = 0.695
CRST_K = 47.5

bl,rcl = continuity_f2(CRST_RLOW,CRST_R0,CRST_RC)

CRST_BLOW = bl
CRST_RCLOW = rcl

print("b low: ")
print(bl)
print("rc low: ")
print(rcl)

bh,rch = continuity_f2(CRST_RHIGH,CRST_R0,CRST_RC)

CRST_BHIGH = bh
CRST_RCHIGH = rch

print("b high: ")
print(bh)
print("rc high: ")
print(rch)


V_f2_1 = []
erre = []
for t in range(0,2200) :
     r = t/2000.
     erre.append(r)
     V_f2_1.append(f2(r,CRST_K,CRST_R0,CRST_RC,CRST_BLOW,CRST_RCLOW,CRST_BHIGH,CRST_RCHIGH,CRST_RLOW,CRST_RHIGH))

 
####################################
#print 
   
fig = plt.figure(figsize=(5,5))
ax =fig.add_axes([0., 0., 1.0, 1.0])
#plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)

ax.set_ylabel(r"$f2(r)$",fontsize=20)
ax.set_xlabel(r"$r$ [oxDNA units]",fontsize=20)
#ax.set_ylim(-0.1,1.2)
ax.set_xlim(-0.1,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
    
ax.scatter(erre,V_f2_1,color="black")
#ax.plot(erre,V_f2_3,color="red",label=r"$a=5$")
#ax.plot(erre,V_f2_1,color="silver",label=r"$a=6$")
#ax.plot(erre,V_f2_4,color="blue",label=r"$a=7$")
#ax.plot(erre,V_f2_5,color="magenta",label=r"$a=8$")

#ax.legend(bbox_to_anchor=(0.75, 0.7), prop={'size':24})
