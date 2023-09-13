#!/usr/bin/env python
# coding: utf-8

import itertools
import numpy as np
import math 
import re
import sys
from os import path
import os
import subprocess

# handle the arguments
if len(sys.argv) < 3:
	print 'USAGE:',sys.argv[0],'<trajectory file> <number of nucleotides>  <alpha>'
	sys.exit(-1)

#PLEASE Note:ALPHA is com correction. If you are not interested in applying that, please use alpha = 0.00. The optimium value to get correct com, alpha= 0.06(units of sim length).

n_nucls = int(sys.argv[2]) #no. of nucleotides
N_bp = n_nucls/2 #no of basepairs

# Fetch nucleotide indices plane by plane
ids_nucls = [] 
for i in range(0,N_bp):
    id1=[i] if i < N_bp else None
    id2=[n_nucls-1-i] if i < N_bp else None
    ids_nucls.append(zip(id1,id2))
    
n_com = len(ids_nucls) #no. of centre of masses
print ("no of basepairs =" , n_com)

bad_words = ['t','b','E']

#This can be done elegantly and faster. I wanted to remove lines containing t,b and E to read the trajectory. So I just wrote a new file without those.

with open(sys.argv[1]) as oldfile, open('trajectory_final.dat', 'w') as newfile:
    for line in oldfile:
        if not any(bad_word in line for bad_word in bad_words):
            newfile.write(line)

#loading as array
traject = np.loadtxt('trajectory_final.dat')
s = (len(traject)/n_nucls)

cms = np.zeros(( s, n_nucls, 3), dtype=np.float64)
b = np.zeros(( s, n_nucls, 3), dtype=np.float64)
n = np.zeros(( s, n_nucls, 3), dtype=np.float64)

#reading the r , b and n values from the trajectory file
for i in range (s):
    cms[i] = traject[(i*n_nucls):(i+1)*n_nucls, :3]
    b[i] = traject[(i*n_nucls):(i+1)*n_nucls, 3:6]
    n[i] = traject[(i*n_nucls):(i+1)*n_nucls, 6:9]

      
#get R_duples =(r_nuc1 + r_nuc2)/2 
def get_centres(cms):
    R_dup = np.zeros((s,n_com,3), dtype=np.float64)
    for id_com,ids_pairs in enumerate(ids_nucls):
        n_pairs = 0
        alpha = float(sys.argv[3])
        for i in range (s):
            for ids in ids_pairs:
                id1,id2  = ids
                r_nuc_1 = cms[i][id1]
                r_nuc_2 = cms[i][id2]
                n_nuc_1 = n[i][id1]
                n_nuc_2 = n[i][id2]
                b_nuc_1 = b[i][id1]
                b_nuc_2 = b[i][id2]
                n_pairs += 1
                R_dup[i][id_com] += 0.5*(r_nuc_1 + r_nuc_2) + alpha *0.5*(np.cross(b_nuc_1,n_nuc_1) + np.cross(b_nuc_2,n_nuc_2))
    return R_dup
R_dup  = get_centres(cms)

#a(helical rise)
def a(R_dup):
    a = np.zeros((s,n_com-1), dtype=np.float64)
    for i in range (s):
		for j in range (n_com-1):
			a[i][j] = 0.8518*np.linalg.norm(R_dup[i][j+1]- R_dup[i][j])  #now nm after 0.8518
    return a
a = a(R_dup)


#y = r_nuc_1-r_nuc_2/||r_nuc_1-r_nuc_2||
def y(cms):
    y = np.zeros((s,n_com,3), dtype=np.float64)
    for id_com,ids_pairs in enumerate(ids_nucls):
        n_pairs = 0
        for i in range (s):
            for ids in ids_pairs:
                id1,id2  = ids
                r_nuc_1 = cms[i][id1]
                r_nuc_2 = cms[i][id2]
                n_pairs += 1
                y[i][id_com] += (r_nuc_2 - r_nuc_1) / np.linalg.norm(r_nuc_2 - r_nuc_1)
    return y
y = y(cms)

#Just random checks to see if things are being taken properly
R_bp = np.zeros((s,n_com-2,3), dtype=np.float64)
for i in range (s):
    R_bp[i] += R_dup[i][1:N_bp-1]

y_f = np.zeros((s,n_com-2,3), dtype=np.float64)
for i in range (s):
     y_f[i] += y[i][1:N_bp-1]

id_nucls = [] 
for i in range(0,N_bp-2):
    i1=[i] if i < N_bp else None
    i2=[i+2] if i < N_bp else None
    id_nucls.append(zip(i1,i2))    


def e3(R_dup):
    e3 = np.zeros((s,n_com-2,3), dtype=np.float64)
    for id_com1,id_pair in enumerate(id_nucls):
        n_pairs = 0
        for i in range (s):
            for ids2 in id_pair:
                i1,i2 = ids2
                R1 = R_dup[i][i1]
                R2 = R_dup[i][i2]
                n_pairs += 1
                e3[i][id_com1] += (R2-R1)/np.linalg.norm(R2-R1)
    return e3
e3 = e3(R_dup)

#e2= (y -(y.e3)*e3)/Mod of numerator
def e2(y_f,e3):
    e2 = np.zeros((s,n_com-2,3), dtype=np.float64)
    for j in range (s):
        for i in range(N_bp-2):
            yf = y_f[j][i]
            e3f = e3[j][i]
            e2[j][i] += (yf - ((np.dot(yf,e3f))*e3f))/(np.linalg.norm(yf - ((np.dot(yf,e3f))*e3f)))
    return e2
e2 = e2(y_f,e3)

e1 = np.cross(e2,e3) #e1=e2Xe3

def T(e1,e2,e3):
    T = np.zeros((s,N_bp-2,3,3), dtype=np.float64)
    for j in range (s):
        for i in range (0,N_bp-2):
            T[j][i][:,0] = e1[j][i]
            T[j][i][:,1] = e2[j][i]
            T[j][i][:,2] = e3[j][i]
    return T
T = T(e1,e2,e3)

def TT(e1,e2,e3):
    TT = np.zeros((s,N_bp-2,3,3), dtype=np.float64)
    for j in range (s):
        for i in range (0,N_bp-2):
            TT[j][i][0] = e1[j][i]
            TT[j][i][1] = e2[j][i]
            TT[j][i][2] = e3[j][i]
    return TT
TT = TT(e1,e2,e3)

def R(T,TT):
    R = np.zeros((s,N_bp-3,3,3), dtype=np.float64)
    for j in range (s):
        for i in range (0,N_bp-3):
            tt = TT[j][i]
            t = T[j][i+1]
            R[j][i] = np.mat(tt)*np.mat(t)
    return R
R = R(T,TT) 

def Trace(R):
    Trace = np.zeros((s,N_bp-3), dtype=np.float64)
    for j in range(s):
        for i in range(0,N_bp-3):
            Trace[j][i] = np.trace(R[j][i]) 
    return Trace
Trace = Trace(R)

def theta(Trace):
    theta = np.zeros((s,N_bp-3), dtype=np.float64)
    for j in range (s):
        for i in range(0,N_bp-3):
            theta[j][i] = (math.acos((Trace[j][i]-1)/2 ))
    return theta
theta = theta(Trace)

def theta11(R,theta):
    theta11 = np.zeros((s,N_bp-3), dtype=np.float64)
    for j in range (s):
        for i in range(0,N_bp-3):
            theta11[j][i] =(0.5*(theta[j][i]/np.sin(theta[j][i]))*(R[j][i][2,1] -R[j][i][1,2]))
    return theta11
theta11 = theta11(R,theta)

def theta22(R,theta):
    theta22 = np.zeros((s, N_bp-3), dtype=np.float64)
    for j in range (s):
        for i in range(0,N_bp-3):
            theta22[j][i]= (0.5*(theta[j][i]/np.sin(theta[j][i]))*(R[j][i][0,2] -R[j][i][2,0]))
    return theta22
theta22 = theta22(R,theta)

def theta33(R,theta):
    theta33 = np.zeros((s, N_bp-3), dtype=np.float64)
    for j in range (s):
        for i in range(0,N_bp-3):
            theta33[j][i]= (0.5*(theta[j][i]/np.sin(theta[j][i]))*(R[j][i][1,0] -R[j][i][0,1]))
    return theta33
theta33 = theta33(R,theta)


a1 = (np.mean(theta11))
a2 = (np.mean(theta22))
a3 = (np.mean(theta33))


a4 = np.mean(theta22 , axis = 0)
a5 = np.mean(theta33 , axis = 0)



omega1 = (theta11 -a1)/(np.mean(a[0:][10:-10]))
omega2 = (theta22 -a2)/(np.mean(a[0:][10:-10]))
omega3 = (theta33 -a3)/(np.mean(a[0:][10:-10]))

def omegas(omega1,omega2,omega3):
    omegas = np.zeros((s,N_bp-3, 3),dtype=np.float64)
    for i in range (s):
        for j in range (N_bp-3):
            omegas[i][j][0] = omega1[i][j]
            omegas[i][j][1] = omega2[i][j]
            omegas[i][j][2] = omega3[i][j]
    return omegas
omegas = omegas(omega1,omega2,omega3)

with open ("theta_deg" , "a+") as f:
                print >> f , '{:^18}{:^18}{:^18}{:^18}{:^18}'.format(("alpha"),("n_nucls"),("t1"),("t2"),("t3"))
        	print >> f , '{:^18}{:^18}{:^18}{:^18}{:^18}'.format((sys.argv[3]),(n_nucls),round(np.rad2deg(a1),4),round(np.rad2deg(a2),4),round(np.rad2deg(a3),4))

print (np.mean(a[0:][10:-10]))
np.save('theta2_alpha_'+ sys.argv[3], a4) #rad
np.save('theta3_alpha_'+ sys.argv[3], a5) #rad
np.save('Omegas_alpha_'+ sys.argv[3], omegas) #rad
np.save('e3_alpha_'+ sys.argv[3], e3) #rad
np.save('Omega3_alpha_'+ sys.argv[3], omega3) #rad
np.save('a_alpha_'+sys.argv[3] , a)  #mean of this file gives exact value of "a" in nm
