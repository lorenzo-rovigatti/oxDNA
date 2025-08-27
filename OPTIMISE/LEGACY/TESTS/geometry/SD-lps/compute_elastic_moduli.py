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
ids = [9,10,11]

Njuns = 99 #total number of junctions

in_j = 3 #ignore ends
fin_j = Njuns -3 #ignore ends
in_snap = 50 #ignore first in_snap snapshots (equilibration)

dimension = (Njuns-in_j-(Njuns-fin_j))*len(ids)

N_fit =  (Njuns-in_j-(Njuns-fin_j))

internal_coords = [] # internal coordinates

#average internal coordinates
mu_sampled = np.zeros(dimension, dtype = float)
mu_global_sampled = np.zeros(len(ids), dtype = float)
#cov_sampled = np.zeros((dimension,dimension), dtype = float)

Ms = [] # Stiffness matrices in Furier space

#fitting function parameters
X_At = []
X_At_err = []
X_Ar = []
X_Ar_err = []
X_C = []
X_C_err = []
X_G = []
X_G_err = []

#oxDNA2 pars
X_At_ox2 = [54,17,4.0,1.1]
X_Ar_ox2 = [38,2,0.8,0.2]
X_C_ox2 = [78,22,6.5,1.3]
X_G_ox2 = [23,6.0,1.9,0.4]


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


# Computes Delta(q) for a given snapshot and q, where Delta is the displacement of the chosen set of internal coordinates 
def Displacement_Fourier(q,snap_id) :
    
    disp_q = np.zeros(len(ids),dtype=complex)
    
    Ncoords = len(internal_coords[0])   
    
    #ave_mu[i%len(ids)] += mu[i]/(1.*len(internal_coords[0])/(1.*len(ids)))
    
    n = -1
    N = Ncoords/len(ids)
    
    for i in range(Ncoords) :
            
            if i%len(ids) == 0 :
                n+=1;
            disp_q[i%len(ids)]  += (internal_coords[snap_id][i]-mu_sampled[i])*math.pi/180/0.34*complex(math.cos(-2*math.pi*q*n/N),math.sin(-2*math.pi*q*n/N))
            #disp_q[i%len(ids)]  += (internal_coords[snap_id][i]-mu_global_sampled[i%len(ids)])*math.pi/180./0.34*complex(math.cos(-2*math.pi*q*n/N),math.sin(-2*math.pi*q*n/N))
    return disp_q


# Print the stiffness matrices to a file
def print_Ms(dim,Nqs,N) :
    
    global Ms
    
    q = -N/2
    
    oname = "Mq_vs_q.txt"
    
    ofile = open(oname,'w')
    
    for i in range(Nqs) :
        if i == 0:
            print("#q, Mq[0,0], Mq[0,1], Mq[0,2], Mq[1,1], Mq[1,2], Mq[2,2]", file=ofile)   
            
        string = str(q)
        for l in range(dim) :
            for m in range(l,dim) :
                if (l == 0 and m == 1) or (l == 0 and m == 1) :
                    string += " "+str(Ms[i][l,m].imag)
                else :
                    string += " "+str(Ms[i][l,m].real)
        print(string,file=ofile)
                
        q += N/Nqs
        
    ofile.close()
        
    return

#fitting function
def fit_f(x,X0,X1,X2,X3) :
    
    return X0+X1*np.cos(2*math.pi*x/N_fit)+X2*np.cos(2*math.pi*2*x/N_fit)+X3*np.cos(2*math.pi*3*x/N_fit)
    
# Print the persistence lengths to a file
def print_lp(ms,lbs,lts) :
    
   
    oname = "lps_vs_m.txt"
    
    ofile = open(oname,'w')
    
    for i in range(len(ms)) :
        if i == 0:
            print("#m, lb, lp", file=ofile)   
            
        string = str(ms[i]) + " " + str(lbs[i]) + " " + str(lts[i])

        print(string,file=ofile)
        
    ofile.close()
        
    return


#Factors used to compute the persistence lengths

m_int = 0

def integrand1(y) :
    
    #print("int1: " +str(m_int))
    
    q = N_fit/math.pi*y
    C = fit_f(q,X_C[0],X_C[1],X_C[2],X_C[3])
    G = fit_f(q,X_G[0],X_G[1],X_G[2],X_G[3])
    Ar = fit_f(q,X_Ar[0],X_Ar[1],X_Ar[2],X_Ar[3])
    if y == 0 :
        return 1./(2*m_int*math.pi)/(C-G*G/Ar)
    else :
        return 1./(2*m_int*math.pi)*math.sin(m_int*y)*math.sin(m_int*y)/math.sin(y)/math.sin(y)/(C-G*G/Ar)


def aphi_N(y) :
    
    aw0 = mu_global_sampled[2]/180*math.pi
    q = N_fit/math.pi*y
    
    C = fit_f(q,X_C[0],X_C[1],X_C[2],X_C[3])
    G = fit_f(q,X_G[0],X_G[1],X_G[2],X_G[3])
    Ar = fit_f(q,X_Ar[0],X_Ar[1],X_Ar[2],X_Ar[3])
    At = fit_f(q,X_At[0],X_At[1],X_At[2],X_At[3])
    
    return (1-math.cos(aw0))/2/aw0/aw0*(1./At + 1./(Ar-G*G/C))

def integrand2(y) :
    
    Dy = mu_global_sampled[2]/180*math.pi/2

    if y == 0 :
        return 1./math.pi/m_int*(aphi_N(y+Dy)+aphi_N(y-Dy))
    else :
        return 1./math.pi/m_int*math.sin(m_int*y)*math.sin(m_int*y)/math.sin(y)/math.sin(y)*(aphi_N(y+Dy)+aphi_N(y-Dy))

#Compute persistence length
def persistence_lengths(N,Delta) :
    
    global m_int   
    
    lbs = []
    lts = []
    ms = []
    
    for i in range(0,N) :
        m_int = 1+i*Delta
        print("m: "+str(m_int))
        ms.append(m_int)
        
        lt = 1./quad(integrand1,-math.pi/2.,math.pi/2.)[0]
        lb = 1./quad(integrand2,-math.pi/2.,math.pi/2.)[0]
        
        lts.append(lt)
        lbs.append(lb)  
    
    return ms,lbs,lts

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

 
#compute average coordinates
ave_stored()


#compute stiffness matrices
N = len(internal_coords[0])/len(ids)

q = -float(N/2)
Nqs = 100

for i in range(Nqs):
    
    print(q)
        
    M = np.zeros((len(ids),len(ids)),dtype=complex)
    
    #compute displacement covariance
    for snap in range(len(internal_coords)) :
    
        disp_q = Displacement_Fourier(q,snap)
    
        for i in range(len(ids)) :
            for j in range(i,len(ids)) :
                    M[i,j] += disp_q[i]*disp_q[j].conjugate()/len(internal_coords)
                    
        for i in range(len(ids)) :
            for j in range(i,len(ids)) :
                M[j,i] = M[i,j].conjugate()                
    
    #from covariance to stiffness matrix
    M = M *(0.34/N)
    M = np.linalg.inv(M)
    
    Ms.append(M)   
    
    q += 1.*float(N)/Nqs
    
    
print_Ms(len(ids),Nqs,N)

print(mu_global_sampled)


#fitting

x = np.zeros(Nqs,dtype=float)
q_rescaled = np.zeros(Nqs,dtype=float)
y0 = np.zeros(Nqs,dtype=float)
y1 = np.zeros(Nqs,dtype=float)
y2 = np.zeros(Nqs,dtype=float)
y3 = np.zeros(Nqs,dtype=float)
y4 = np.zeros(Nqs,dtype=float)
y5 = np.zeros(Nqs,dtype=float)
q = -float(N/2)

for i in range(Nqs) :
    x[i] = q
    q_rescaled[i] = q*math.pi/N
    y0[i] = Ms[i][0,0].real
    y1[i] = Ms[i][1,1].real
    y2[i] = Ms[i][2,2].real
    y3[i] = Ms[i][2,1].real
    y4[i] = Ms[i][0,1].imag
    y5[i] = Ms[i][0,2].imag
    q += 1.*float(N)/Nqs

    

X_At, X_At_err = curve_fit(fit_f, x, y0)
X_Ar, X_Ar_err = curve_fit(fit_f, x, y1)
X_C, X_C_err = curve_fit(fit_f, x, y2)
X_G, X_G_err = curve_fit(fit_f, x, y3)

q_rescaled_fit = np.zeros(4*Nqs,dtype=float)

q = -float(N/2)
for i in range(4*Nqs) :
    q_rescaled_fit[i] = q*math.pi/N
    q += 1.*float(N)/Nqs/4.
    

fig = plt.figure(figsize=(4, 6))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$\pi q/N$",fontsize=20)
ax.set_ylabel(r"[nm]",fontsize=20)
ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.scatter(q_rescaled,y0,label="At")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_At[0], X_At[1], X_At[2], X_At[3]), 'b--', label=None)
ax.scatter(q_rescaled,y1,label="Ar")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_Ar[0], X_Ar[1], X_Ar[2], X_Ar[3]), 'b--', label=None)
ax.scatter(q_rescaled,y2,label="C")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_C[0], X_C[1], X_C[2], X_C[3]), 'b--', label=None)
ax.scatter(q_rescaled,y3,label="G")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_G[0], X_G[1], X_G[2], X_G[3]), 'b--', label=None)
ax.legend(fontsize = 20)
plt.savefig("M_q_vs_q.pdf",bbox_inches='tight',pad_inches=0.05)

fig = plt.figure(figsize=(4, 6))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$\pi q/N$",fontsize=20)
ax.set_ylabel(r"[nm]",fontsize=20)
ax.set_ylim(0,160)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_At_ox2[0], X_At_ox2[1], X_At_ox2[2], X_At_ox2[3]), 'k-', label="ox2")
ax.scatter(q_rescaled,y0,label="At")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_At[0], X_At[1], X_At[2], X_At[3]), 'b--', label=None)
ax.scatter(q_rescaled,y1,label="Ar")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_Ar[0], X_Ar[1], X_Ar[2], X_Ar[3]), 'b--', label=None)
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_Ar_ox2[0], X_Ar_ox2[1], X_Ar_ox2[2], X_Ar_ox2[3]), 'k-', label=None)
ax.scatter(q_rescaled,y2,label="C")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_C[0], X_C[1], X_C[2], X_C[3]), 'b--', label=None)
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_C_ox2[0], X_C_ox2[1], X_C_ox2[2], X_C_ox2[3]), 'k-', label=None)
ax.scatter(q_rescaled,y3,label="G")
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_G[0], X_G[1], X_G[2], X_G[3]), 'b--', label=None)
ax.plot(q_rescaled_fit, fit_f(q_rescaled_fit*N/math.pi, X_G_ox2[0], X_G_ox2[1], X_G_ox2[2], X_G_ox2[3]), 'k-', label=None)
ax.legend(fontsize = 20)
plt.savefig("M_q_vs_q_vs_ox2.pdf",bbox_inches='tight',pad_inches=0.05)


#plot other two non-diagonal elements
fig = plt.figure(figsize=(4, 6))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$\pi q/N$",fontsize=20)
ax.set_ylabel(r"[nm]",fontsize=20)
#ax.set_ylim(0,120)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.scatter(q_rescaled,y4,label="B")
ax.scatter(q_rescaled,y5,label="Atr")
ax.legend(fontsize = 20)
plt.savefig("M_q_others_vs_q.pdf",bbox_inches='tight',pad_inches=0.05)

#compute persistence lengths (using fit functions for the element of the stiffness matrices)
ms,lbs,lts = persistence_lengths(100, 1)

print_lp(ms,lbs,lts)


lts_half = np.array(lts)
lts_half = lts_half*0.5

#Lps = np.loadtxt("ox2_lps_vs_m.txt")


CCS_G = []
CCS_G.append('#000000')
CCS_G.append('#006884')
CCS_G.append('#00909E')
CCS_G.append('#5B5B5B')
CCS_G.append('#6E006C')
CCS_G.append('#89DBEC')
CCS_G.append('#91278F')
CCS_G.append('#B00051')
CCS_G.append('#CF97D7')
CCS_G.append('#D4D4D4')
CCS_G.append('#ED0026')
CCS_G.append('#F68370')
CCS_G.append('#FA9D00')
CCS_G.append('#FEABB9')
CCS_G.append('#FFD08D')

#style

plt.style.use("../../style.sty")

fig = plt.figure(figsize=(4, 4))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
# plt.title(r"All: $R=7$ $\phi=0.364$")
ax.title.set_fontsize(20)
ax.set_xlabel(r"$m$",fontsize=20)
ax.set_ylabel(r"persistence length [nm]",fontsize=20)
ax.set_ylim(30,140)
#ax.set_xlim(-1.2,1.2)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.tick_params(axis='both', which='minor', labelsize=8)
#ax.plot(Lps[:,0],Lps[:,1], "k-", label="ox2")
#ax.plot(Lps[:,0],Lps[:,2]*0.5, "k-", label=None)
ax.plot(ms,lbs,color=CCS_G[1],label="$l_b$")
ax.plot(ms,lts_half,color=CCS_G[7],label="$l_t/2$")
ax.legend(fontsize = 20)
plt.savefig("lps.pdf",bbox_inches='tight',pad_inches=0.05)


#m to infinity

C = fit_f(0,X_C[0],X_C[1],X_C[2],X_C[3])
G = fit_f(0,X_G[0],X_G[1],X_G[2],X_G[3])
Ar = fit_f(0,X_Ar[0],X_Ar[1],X_Ar[2],X_Ar[3])
At = fit_f(0,X_At[0],X_At[1],X_At[2],X_At[3])

lt_half_inf = (C-G*G/Ar)
aw0 = mu_global_sampled[2]/180*math.pi
lb_inf = aw0*aw0/(1-math.cos(aw0))/(1./At + 1./(Ar-G*G/C))

print("m to infinity:")
print("lb = "+str(lb_inf))
print("lt/2 = "+str(lt_half_inf))
