#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:53:18 2024

@author: yqb22156
"""

import math
import numpy as np
import sys

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


def check_types_number(pars, dim, hydr=False) :
    
    
    allgood = 0
    number = 4*4
    if dim == 4:
        number = number *4*4
        
    if dim != 2 and dim != 4:
        print("Dimension is not 2 nor 4")
        allgood = 1
    
    if hydr: 
        number = 4
    
    if len(pars) != number :
        print("Wrong number of parameters; expected ", number)
        print("Got: ", len(pars))
        allgood = 1

    return allgood


bases = ['A', 'C', 'G', 'T']

def get_compl_type(ty):
    id_ty = []
    for t in ty:
        if t == 'A':
            id_ty.append(0)
        elif t == 'C':
            id_ty.append(1)
        elif t == 'G':
            id_ty.append(2)
        elif t == 'T':
            id_ty.append(3)
            
    ty_c = ""
    for l in range(len(id_ty)): 
        ty_c += bases[3-id_ty[len(id_ty)-l-1]]
    return ty_c

#symm_type = 0 : no symmetry
#symm_type = 1 : reverse symmetry (i.e. AT = TA)
#symm_type = 2 : complementary symmetry (i.e. AG = CT)
def check_symmetry(pars, symm_type) :
    
   if symm_type == 0:
       return 0
   
   allgood = 0 
   
   for par in pars :
       if symm_type == 1:
           par_s = par[::-1]
           #print(par,par_s)
           #print(pars[par],pars[par_s])
           if pars[par] > pars[par_s]+1e-12 or pars[par] < pars[par_s]-1e-12 :
               print("Something weird for type: ", par)
               print("Expected reverse symmetry: "+par+" = "+par_s)
               print("Got : "+str(pars[par])+", "+str(pars[par_s]))
               allgood = 1
       if symm_type == 2:
            par_s = get_compl_type(par)
            if pars[par] > pars[par_s]+1e-12 or pars[par] < pars[par_s]-1e-12 :
                print("Something weird for type: ", par)
                print("Expected complementary symmetry: "+par+" = "+par_s)
                print("Got : "+str(pars[par])+", "+str(pars[par_s]))
                allgood = 1
   return allgood

    
def check_continuity_f1(pars_r0, pars_a, pars_rc, pars_bl, pars_bh, pars_rl, pars_rh, pars_rcl, pars_rch ) :
    
    allgood = 0
    
    for par in pars_r0 :
        
        r0 =  pars_r0[par]
        a = pars_a[par]
        
        rl = r0-0.08
        rh = r0+0.35
        rc = r0+0.5
        
        num = a*a*(math.exp(-(rl-r0)*a)-math.exp(-2*(rl-r0)*a))**2
        den = (1-math.exp(-(rl-r0)*a))**2-(1-math.exp(-(rc-r0)*a))**2
        
        bl = num/den
        
        num = a*a*(math.exp(-(rh-r0)*a)-math.exp(-2*(rh-r0)*a))**2
        den = (1-math.exp(-(rh-r0)*a))**2-(1-math.exp(-(rc-r0)*a))**2
        
        bh = num/den
        
        rcl = rl-a/bl*(math.exp(-(rl-r0)*a)-math.exp(-2*(rl-r0)*a))
        rch = rh-a/bh*(math.exp(-(rh-r0)*a)-math.exp(-2*(rh-r0)*a))
        
        if  pars_rc[par] > rc+1e-6 or pars_rc[par] < rc-1e-6 :
               print("Something weird for rc, type: ", par)
               print("Expected: ",rc)
               print("Got : ",pars_rc[par])
               allgood = 1
        
        if  pars_bl[par] > bl+1e-6 or pars_bl[par] < bl-1e-6 :
               print("Something weird for bl, type: ", par)
               print("Expected: ",bl)
               print("Got : ",pars_bl[par])
               allgood = 1
               
        if  pars_bh[par] > bh+1e-6 or pars_bh[par] < bh-1e-6 :
               print("Something weird for bh, type: ", par)
               print("Expected: ",bh)
               print("Got : ",pars_bh[par])
               allgood = 1
               
        if  pars_rl[par] > rl+1e-6 or pars_rl[par] < rl-1e-6 :
               print("Something weird for rl, type: ", par)
               print("Expected: ",rl)
               print("Got : ",pars_rl[par])
               allgood = 1
               
        if  pars_rh[par] > rh+1e-6 or pars_rh[par] < rh-1e-6 :
               print("Something weird for rh, type: ", par)
               print("Expected: ",rh)
               print("Got : ",pars_rh[par])
               allgood = 1

        if  pars_rcl[par] > rcl+1e-6 or pars_rcl[par] < rcl-1e-6 :
               print("Something weird for rcl, type: ", par)
               print("Expected: ",rcl)
               print("Got : ",pars_rcl[par])
               allgood = 1

        if  pars_rch[par] > rch+1e-6 or pars_rch[par] < rch-1e-6 :
               print("Something weird for rch, type: ", par)
               print("Expected: ",rch)
               print("Got : ",pars_rch[par])
               allgood = 1
               
    return allgood

def check_continuity_f2(pars_r0, pars_rc, pars_bl, pars_bh, pars_rl, pars_rh, pars_rcl, pars_rch ) :
    
    allgood = 0
    
    for par in pars_r0 :
        
        r0 =  pars_r0[par]
        
        rl = r0-0.08
        rh = r0+0.08
        rc = r0+0.1
        
        rcl = rl - 1./(rl-r0)*( (rl-r0)*(rl-r0)-(rc-r0)*(rc-r0) )
        rch = rh - 1./(rh-r0)*( (rh-r0)*(rh-r0)-(rc-r0)*(rc-r0) )
        
        bl = 0.5*(r0-rl)/(rcl-rl)
        bh = 0.5*(r0-rh)/(rch-rh)
        
        if  pars_rc[par] > rc+1e-6 or pars_rc[par] < rc-1e-6 :
               print("Something weird for rc, type: ", par)
               print("Expected: ",rc)
               print("Got : ",pars_rc[par])
               allgood = 1
        
        if  pars_bl[par] > bl+1e-6 or pars_bl[par] < bl-1e-6 :
               print("Something weird for bl, type: ", par)
               print("Expected: ",bl)
               print("Got : ",pars_bl[par])
               allgood = 1
               
        if  pars_bh[par] > bh+1e-6 or pars_bh[par] < bh-1e-6 :
               print("Something weird for bh, type: ", par)
               print("Expected: ",bh)
               print("Got : ",pars_bh[par])
               allgood = 1
               
        if  pars_rl[par] > rl+1e-6 or pars_rl[par] < rl-1e-6 :
               print("Something weird for rl, type: ", par)
               print("Expected: ",rl)
               print("Got : ",pars_rl[par])
               allgood = 1
               
        if  pars_rh[par] > rh+1e-6 or pars_rh[par] < rh-1e-6 :
               print("Something weird for rh, type: ", par)
               print("Expected: ",rh)
               print("Got : ",pars_rh[par])
               allgood = 1

        if  pars_rcl[par] > rcl+1e-6 or pars_rcl[par] < rcl-1e-6 :
               print("Something weird for rcl, type: ", par)
               print("Expected: ",rcl)
               print("Got : ",pars_rcl[par])
               allgood = 1

        if  pars_rch[par] > rch+1e-6 or pars_rch[par] < rch-1e-6 :
               print("Something weird for rch, type: ", par)
               print("Expected: ",rch)
               print("Got : ",pars_rch[par])
               allgood = 1
               
    return allgood
        
    
def check_continuity_f4(pars_a, pars_b, pars_dts, pars_dtc ) :
    
    allgood = 0
    
    for par in pars_a :
        
        a =  pars_a[par]
        
        dts = math.sqrt(0.81225/a) #delta theta star; smoothing kicks in at fmod = 0.18775 (as in STCK oxdna2)
        dtc = 1./(a*dts) #delta theta c
        b = a*dts/(dtc-dts) #width b
        
        if  pars_b[par] > b+1e-6 or pars_b[par] < b-1e-6 :
               print("Something weird for b, type: ", par)
               print("Expected: ",b)
               print("Got : ",pars_b[par])
               allgood = 1
               
        if  pars_dts[par] > dts+1e-6 or pars_dts[par] < dts-1e-6 :
               print("Something weird for dts, type: ", par)
               print("Expected: ",dts)
               print("Got : ",pars_dts[par])
               allgood = 1
               
        if  pars_dtc[par] > dtc+1e-6 or pars_dtc[par] < dtc-1e-6 :
               print("Something weird for dtc, type: ", par)
               print("Expected: ",dtc)
               print("Got : ",pars_dtc[par])
               allgood = 1
        
        return allgood
    

if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 "+sys.argv[0]+" par_file")
    sys.exit()

par_file = open(sys.argv[1], 'r')


params_fene = {}
params_stck = {}
params_hydr = {}
params_crst_33 = {}
params_crst_55 = {}

#read parameters file

for line in par_file.readlines() :
    vals = line.strip().split()
    if len(vals) == 0:
        continue
    if vals[0][0] == '#':
        continue
    name_vals = vals[0].split('_')
    if name_vals[0] == "FENE":
        ty = ""
        pname = ""
        for l in range(4):
            ty += name_vals[2+l]
        pname = name_vals[0] + "_" + name_vals[1]
        if pname not in params_fene.keys() :
            params_fene[pname] = {}
        
        params_fene[pname][ty] = float(vals[2])
    
    if name_vals[0] == "STCK":
        ty = ""
        pname = ""
        if len(name_vals) == 3 :
            if name_vals[1] == "FACT" :
                continue
            else:
                for l in range(2):
                    ty += name_vals[1+l]
                pname = name_vals[0]
                if pname not in params_stck.keys() :
                    params_stck[pname] = {}
        else:     
            offset = 2
            if len(name_vals) == 7: offset += 1
            for l in range(4):
                ty += name_vals[offset+l]
            pname = name_vals[0]
            for l in range(1,offset) :
                pname += "_"+name_vals[l]
        if pname not in params_stck.keys() :
            params_stck[pname] = {}
        
        params_stck[pname][ty] = float(vals[2])
        
    if name_vals[0] == "HYDR":
        ty = ""
        pname = ""
        if len(name_vals) == 3 :
            for l in range(2):
                ty += name_vals[1+l]
            pname = name_vals[0]
            if pname not in params_hydr.keys() :
                params_hydr[pname] = {}
        else:     
            offset = 2
            if len(name_vals) == 5: offset += 1
            for l in range(2):
                ty += name_vals[offset+l]
            pname = name_vals[0]
            for l in range(1,offset) :
                pname += "_"+name_vals[l]
        if pname not in params_hydr.keys() :
            params_hydr[pname] = {}
        
        params_hydr[pname][ty] = float(vals[2])
        
        
    if name_vals[0] == "CRST":
        ty = ""
        pname = ""
        if len(name_vals) == 5 :
            for l in range(2):
                ty += name_vals[3+l]
            pname = name_vals[0]+"_"+name_vals[1]+"_"+name_vals[2]
            if name_vals[2] == "33" :
                if pname not in params_crst_33.keys() :
                    params_crst_33[pname] = {}
                params_crst_33[pname][ty] = float(vals[2])
            elif name_vals[2] == "55" :
                if pname not in params_crst_55.keys() :
                    params_crst_55[pname] = {}
                params_crst_55[pname][ty] = float(vals[2])
        else:     
            offset = 3
            if len(name_vals) == 8: offset += 1
            for l in range(4):
                ty += name_vals[offset+l]
            pname = name_vals[0]
            for l in range(1,offset) :
                pname += "_"+name_vals[l]
            if name_vals[offset-1] == "33":
                if pname not in params_crst_33.keys() :
                    params_crst_33[pname] = {}      
                params_crst_33[pname][ty] = float(vals[2])
            elif name_vals[offset-1] == "55":
                if pname not in params_crst_55.keys() :
                    params_crst_55[pname] = {}      
                params_crst_55[pname][ty] = float(vals[2])


par_file.close()
    
#print(params_fene)
#print(params_stck)
#print(params_hydr)
#print(params_crst_33)
#print(params_crst_55)

#checking that all parameters are included (256 for 4d pars, 16 for 2d pars, 4 for 2d hydro)

print("CHECKING NUMBERS")

for par in params_fene :
    isgood = check_types_number(params_fene[par], len(next(iter(params_fene[par]))) )
    if isgood != 0 :
        print("Something weird for parameter ", par)
for par in params_hydr :
    isgood = check_types_number(params_hydr[par], len(next(iter(params_hydr[par]))), True )
    if isgood != 0 :
        print("Something weird for parameter ", par)
for par in params_stck :
    isgood = check_types_number(params_stck[par], len(next(iter(params_stck[par]))) )
    if isgood != 0 :
        print("Something weird for parameter ", par)
for par in params_crst_33 :
    isgood = check_types_number(params_crst_33[par], len(next(iter(params_crst_33[par]))) )
    if isgood != 0 :
        print("Something weird for parameter ", par)
for par in params_crst_55 :
    isgood = check_types_number(params_crst_55[par], len(next(iter(params_crst_55[par]))) )
    if isgood != 0 :
        print("Something weird for parameter ", par)


#checking symmetries
print("CHECKING SYMMETRIES")

for par in params_fene :
    isgood = check_symmetry(params_fene[par], 0 )
    if isgood != 0 :
        print("Something weird for parameter ", par)
for par in params_hydr :
    isgood = check_symmetry(params_hydr[par], 1 )
    if isgood != 0 :
        print("Something weird for parameter ", par)
for par in params_stck :
    isgood = check_symmetry(params_stck[par], 0 )
    if isgood != 0 :
        print("Something weird for parameter ", par)
for par in params_crst_33 :
    isgood = check_symmetry(params_crst_33[par], 1 )
    if isgood != 0 :
        print("Something weird for parameter ", par)

for par in params_crst_55 :
    isgood = check_symmetry(params_crst_55[par], 1 )
    if isgood != 0 :
        print("Something weird for parameter ", par)


#checking continuity

print("CHECKING CONTINUITY")

#f1
pastck = {}
if "STCK_A" in params_stck : 
    pastck = params_stck["STCK_A"]
else:
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for z in range(4):
                    ty = bases[i]+bases[j]+bases[k]+bases[z] 
                    pastck[ty] = 6.
                    
isgood = check_continuity_f1(params_stck["STCK_R0"],pastck,params_stck["STCK_RC"],params_stck["STCK_BLOW"],params_stck["STCK_BHIGH"], \
                             params_stck["STCK_RLOW"],params_stck["STCK_RHIGH"],params_stck["STCK_RCLOW"],params_stck["STCK_RCHIGH"])
if isgood != 0:
    print("Something weird with stacking f1 continuity")
#f2
isgood = check_continuity_f2(params_crst_33["CRST_R0_33"],params_crst_33["CRST_RC_33"],params_crst_33["CRST_BLOW_33"],params_crst_33["CRST_BHIGH_33"], \
                             params_crst_33["CRST_RLOW_33"],params_crst_33["CRST_RHIGH_33"],params_crst_33["CRST_RCLOW_33"],params_crst_33["CRST_RCHIGH_33"])
if isgood != 0:
    print("Something weird with crst_33 f2 continuity")

isgood = check_continuity_f2(params_crst_55["CRST_R0_55"],params_crst_55["CRST_RC_55"],params_crst_55["CRST_BLOW_55"],params_crst_55["CRST_BHIGH_55"], \
                             params_crst_55["CRST_RLOW_55"],params_crst_55["CRST_RHIGH_55"],params_crst_55["CRST_RCLOW_55"],params_crst_55["CRST_RCHIGH_55"])
if isgood != 0:
    print("Something weird with crst_55 f2 continuity")
    
#f4
isgood = check_continuity_f4(params_hydr["HYDR_THETA4_A"],params_hydr["HYDR_THETA4_B"],params_hydr["HYDR_THETA4_TS"],params_hydr["HYDR_THETA4_TC"])
if isgood != 0:
    print("Something weird with hydr theta4 f4 continuity")
    
isgood = check_continuity_f4(params_stck["STCK_THETA4_A"],params_stck["STCK_THETA4_B"],params_stck["STCK_THETA4_TS"],params_stck["STCK_THETA4_TC"])
if isgood != 0:
    print("Something weird with stck theta4 f4 continuity")
    
    isgood = check_continuity_f4(params_stck["STCK_THETA5_A"],params_stck["STCK_THETA5_B"],params_stck["STCK_THETA5_TS"],params_stck["STCK_THETA5_TC"])
if isgood != 0:
    print("Something weird with stck theta5 f4 continuity")
    
    isgood = check_continuity_f4(params_crst_33["CRST_THETA4_A_33"],params_crst_33["CRST_THETA4_B_33"],params_crst_33["CRST_THETA4_TS_33"],params_crst_33["CRST_THETA4_TC_33"])
if isgood != 0:
    print("Something weird with crst_33 theta4 f4 continuity")
    
isgood = check_continuity_f4(params_crst_55["CRST_THETA4_A_55"],params_crst_55["CRST_THETA4_B_55"],params_crst_55["CRST_THETA4_TS_55"],params_crst_55["CRST_THETA4_TC_55"])
if isgood != 0:
    print("Something weird with crst_55 theta4 f4 continuity")