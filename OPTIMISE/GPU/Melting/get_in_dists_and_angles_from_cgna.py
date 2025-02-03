#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:11:15 2024

@author: yqb22156
"""

dst = 0.26  #stacking dist from center
dh = 0.2
dbb = 1.0
target_ave_delta = 0.2

import sys
import math
import matplotlib.pyplot as plt
from get_cgdna_pars import get_target_mean_and_covariance
import numpy as np
import continuity_constraints as cc

ox_to_Ang = 8.518


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
    th0_3 = th0_2
    th0_7 = math.acos(rise/crst_r0)
    th0_7 = th0_7 + 0.14 #oxdna correction 
    th0_8 = th0_7
    crst_r0 = crst_r0 + 0.024 #oxdna correction 
    
    return r0,R0,crst_r0,th0_1,th0_2,th0_3,th0_7,th0_8


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

r0,R0,crst_r0,th0_1,th0_2,th0_3,th0_7,th0_8=get_pars(rise,twist)

rc,rl,rh,bl,bh,rcl,rch,continuous = cc.continuity_f2_fixed_width(crst_r0)

a = 6
rc = 0.9
r0 = 0.3815265815609129
rl = 0.3015265815609129
rh = 0.6715265815609128
bl = -118.70378868875504
bh = -1850.7809069256127
rcl = 0.2512018826328554
rch = 0.6719957229564483

rl,rh,bl,bh,rcl,rch,continuous = cc.continuity_f1(r0,a,rc)

CRST_R0_33_A_A_A_A = 0.5740150624614532
CRST_RC_33_A_A_A_A = 0.6740150624614532
CRST_RLOW_33_A_A_A_A = 0.49401506246145316
CRST_RHIGH_33_A_A_A_A = 0.6540150624614531
CRST_BLOW_33_A_A_A_A = -0.8888888888888916
CRST_BHIGH_33_A_A_A_A = -0.8888888888888876
CRST_RCLOW_33_A_A_A_A = 0.4490150624614533
CRST_RCHIGH_33_A_A_A_A = 0.6990150624614532


#crst_r0 = CRST_R0_55_A_A_A_A
#rc = CRST_RC_55_A_A_A_A
#rl = CRST_RLOW_55_A_A_A_A
#rh = CRST_RHIGH_55_A_A_A_A
#bl = CRST_BLOW_55_A_A_A_A
#bh = CRST_BHIGH_55_A_A_A_A
#rcl = CRST_RCLOW_55_A_A_A_A
#rch = CRST_RCHIGH_55_A_A_A_A


crst_r0 = CRST_R0_33_A_A_A_A
rc = CRST_RC_33_A_A_A_A
rl = CRST_RLOW_33_A_A_A_A
rh = CRST_RHIGH_33_A_A_A_A
bl = CRST_BLOW_33_A_A_A_A
bh = CRST_BHIGH_33_A_A_A_A
rcl = CRST_RCLOW_33_A_A_A_A
rch = CRST_RCHIGH_33_A_A_A_A


out_file = "f1.pdf"



#rc,rl,rh,bl,bh,rcl,rch,continuous = cc.continuity_f2_fixed_width(crst_r0)

cc.plot_modulation_f2(crst_r0,rc,bl,rcl,bh,rch,rl,rh,out_file) 
#cc.plot_modulation_f1(r0,a,rc,bl,rcl,bh,rch,rl,rh,out_file)



if len(sys.argv) != 2 :
    print("Unknown argument format.")
    print("Usage: python3 get_in_dists....py sequences_file")
    sys.exit(1)

ifile = open(sys.argv[1],'r')
"""

ifile = open('seq2.txt','r')


seqs = []

rise = np.zeros((4,4,4,4),dtype=float )
twist = np.zeros((4,4,4,4),dtype=float )
deltas = np.zeros((4,4,4,4),dtype=float )
counts = np.zeros((4,4,4,4),dtype=float )

            

for line in ifile.readlines() :
    seqs.append(line.rstrip())
    
ifile.close()
ids = [8,11]
jx_in = 2
jx_from_end = 2
tot_counts = 0
ave_delta = 0.
min_delta = 1.
max_delta = -0.3
resc = 1.
for seq in seqs:
    
    
    red_gs, red_cov = get_target_mean_and_covariance(seq, ids, jx_in, jx_from_end)

    for i in range(jx_in, len(seq)-jx_from_end-1):
        twist[base_to_id(seq[i-1])][base_to_id(seq[i])][base_to_id(seq[i+1])][base_to_id(seq[i+2])] += red_gs[(i-jx_in)*2]
        rise[base_to_id(seq[i-1])][base_to_id(seq[i])][base_to_id(seq[i+1])][base_to_id(seq[i+2])] += red_gs[(i-jx_in)*2+1]
        delta_tmp = (red_cov[(i-jx_in)*2][(i-jx_in)*2]-0.09)/0.53
        deltas[base_to_id(seq[i-1])][base_to_id(seq[i])][base_to_id(seq[i+1])][base_to_id(seq[i+2])] += delta_tmp
        ave_delta += delta_tmp
        tot_counts +=1
        if delta_tmp > max_delta :
            max_delta = delta_tmp
        if delta_tmp < min_delta :
            min_delta = delta_tmp
        counts[base_to_id(seq[i-1])][base_to_id(seq[i])][base_to_id(seq[i+1])][base_to_id(seq[i+2])]+=1
        
ave_delta /= (1.*tot_counts)
#0.1 and 0.45 are the min and max values for fene_delta. Don't want to push them more than that

resc1 = (0.1-target_ave_delta)/(min_delta-ave_delta)
resc2 = (0.45-target_ave_delta)/(max_delta-ave_delta)

if resc1 < resc2 and resc1 < 1.:
    resc = resc1
elif resc2 < resc1 and resc2 < 1.:
    resc = resc2


for i in range(4):
    for j in range(4):
        for l in range(4):
            for m in range(4):                
                if counts[i][j][l][m] == 0:
                    print("warning. untrained tetramer")
                    print(i,j,l,m)
                    rise[i][j][l][m] = 3.35
                    twist[i][j][l][m] = 3.2
                    deltas[i][j][l][m] = target_ave_delta
                else :
                    rise[i][j][l][m] /= counts[i][j][l][m] 
                    twist[i][j][l][m] /= counts[i][j][l][m] 
                    deltas[i][j][l][m] /= counts[i][j][l][m] 
                    #rescale and shift deltas so that average is target_ave_delta and fluctuations are rescaled so that 0.1<delta<0.45
                    deltas[i][j][l][m] = target_ave_delta + resc*(deltas[i][j][l][m]-ave_delta)
                    
                    
                    print(str(bases[i])+str(bases[j])+str(bases[l])+str(bases[m])+" "+str(rise[i][j][l][m])+" "+str(twist[i][j][l][m]))
                   
                    
ofile = open("4D_SD_par_file.txt",'w')

ofile2 = open("4D_SD_par_file_crst_only.txt",'w')

ofile3 = open("4D_SD_par_file_stck_only.txt",'w')

ofile4 = open("4D_SD_par_file_fene_r0_only.txt",'w')

ofile5 = open("4D_SD_par_file_fene_delta_only.txt",'w')

print("printing SD parameters to 4D_SD_par_file.txt")
print("Warning!")
print("We are using: STCK_A = 6; STCK_RC = 0.9")

for i in range(4):
    for j in range(4):
        for l in range(4):
            for m in range(4):  
                
                r0,R0,crst_r0,th0_1,th0_2,th0_3,th0_7,th0_8=get_pars(rise[i][j][l][m],twist[i][j][l][m])
                
                string = "STCK_R0_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(r0)                
                print(string, file=ofile)
                print(string, file=ofile3)                
                
                rl,rh,bl,bh,rcl,rch,continuous = cc.continuity_f1(r0,6.0,0.9)
                               
                string = "STCK_RLOW_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rl)                
                print(string, file=ofile)
                print(string, file=ofile3)
                string = "STCK_RHIGH_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rh)                
                print(string, file=ofile)
                print(string, file=ofile3)
                string = "STCK_BLOW_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(bl)                
                print(string, file=ofile)
                print(string, file=ofile3)
                string = "STCK_BHIGH_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(bh)                
                print(string, file=ofile)
                print(string, file=ofile3)
                string = "STCK_RCLOW_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rcl)                
                print(string, file=ofile)
                print(string, file=ofile3)
                string = "STCK_RCHIGH_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rch)                
                print(string, file=ofile)
                print(string, file=ofile3)
                
                string = "FENE_DELTA_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(deltas[i][j][l][m])                
                print(string, file=ofile)
                print(string, file=ofile5)
                string = "FENE_DELTA2_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(deltas[i][j][l][m]*deltas[i][j][l][m])                
                print(string, file=ofile)
                print(string, file=ofile5)

                string = "FENE_R0_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(R0)                
                print(string, file=ofile)
                print(string, file=ofile4)


                #33 direction
                r0,R0,crst_r0,th0_1,th0_2,th0_3,th0_7,th0_8=get_pars(rise[3-i][3-j][l][m],twist[3-i][3-j][l][m])
                rc,rl,rh,bl,bh,rcl,rch,continuous = cc.continuity_f2_fixed_width(crst_r0)

                string = "CRST_R0_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(crst_r0)                
                print(string, file=ofile)
                print(string, file=ofile2)                  
                
                string = "CRST_RC_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rc)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RLOW_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rl)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RHIGH_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rh)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_BLOW_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(bl)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_BHIGH_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(bh)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RCLOW_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rcl)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RCHIGH_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rch)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                
                """
                string = "CRST_THETA1_T0_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_1)                
                print(string, file=ofile)
                print(string, file=ofile2)    
                """
                string = "CRST_THETA2_T0_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_2)                
                print(string, file=ofile)
                print(string, file=ofile2)                  
                string = "CRST_THETA3_T0_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_3)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                
                string = "CRST_THETA7_T0_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_7)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_THETA8_T0_33_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_8)                
                print(string, file=ofile)
                print(string, file=ofile2) 
                
                
                #55 direction
                r0,R0,crst_r0,th0_1,th0_2,th0_3,th0_7,th0_8=get_pars(rise[i][j][3-l][3-m],twist[i][j][3-l][3-m])
                rc,rl,rh,bl,bh,rcl,rch,continuous = cc.continuity_f2_fixed_width(crst_r0)
                
                string = "CRST_R0_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(crst_r0)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RC_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rc)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RLOW_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rl)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RHIGH_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rh)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_BLOW_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(bl)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_BHIGH_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(bh)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RCLOW_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rcl)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_RCHIGH_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(rch)                
                print(string, file=ofile)
                print(string, file=ofile2)  
         
                """
                string = "CRST_THETA1_T0_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_1)                
                print(string, file=ofile)
                print(string, file=ofile2) 
                """
                string = "CRST_THETA2_T0_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_2)                
                print(string, file=ofile)
                print(string, file=ofile2)   
                string = "CRST_THETA3_T0_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(th0_3)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                string = "CRST_THETA7_T0_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(math.pi-th0_7)                
                print(string, file=ofile)
                print(string, file=ofile2)                   
                string = "CRST_THETA8_T0_55_"+str(bases[i])+"_"+str(bases[j])+"_"+str(bases[l])+"_"+str(bases[m])+" = "+str(math.pi-th0_8)                
                print(string, file=ofile)
                print(string, file=ofile2)  
                
                
ofile.close()
ofile2.close()
ofile3.close()
ofile4.close()

