#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 14:07:04 2022

@author: andrea bonato

code to map oxdna coordinates (center of mass + base vectors)
to internal coordinates (slide, twist, etc.)

it takes a oxdna trajectory and topology as input 
and outputs a trajectory with internal coordinates


"""
import numpy as np
import math
from scipy import linalg

###############################
#Parameteres to map oxdna coordinates to interaction centers
STACK_X = 0.34
HB_X = 0.4
BB_X = -0.34
BB_Y = -0.34

###############################
#parameter to map nucleotide centers to Euler translations
OX_TZU_DISP = [0.0,0.0,0.0]

###############################
#class for oxdna coordinates
#c (center),bv (base-backbone vector),bnv (base normal vector)
#and normal (bv x bnv) are np arrays
#int_centers returns the positions of the interaction centers
class oxdna_frame :
    def __init__(self, c,bv,n):
        self.center = c
        self.base_v = -bv #minus sign for Tsukuba convention
        self.normal = n
        self.base_norv = np.cross(bv,n)
    def int_centers(self):
        back_bone = self.center-BB_X*self.base_v-BB_Y*self.base_norv
        hbond = self.center+HB_X*self.base_v
        stack = self.center+STACK_X*self.base_v
        return(back_bone,hbond,stack)
    
##############################  
#class for euler coordinates
#pos (translation) is a np array, orientation is a np matrix   
class eframe :
    def __init__(self, p,ori):
        self.pos = p
        self.orientation = ori
        
##############################  
#class for interanl (intra or inter) coordinates
#tran (translations) and rot (rotations) are np array   
class int_coord :
    def __init__(self, tr,rot):
        self.tran = tr
        self.rot = rot
 
##############################
#base class i.e. 1 nucleotide and its frames
#members are oxframe, which is an oxdna_frame object
#and frame, which is a eframe object
class base :
    def __init__(self,oxc):
        self.oxframe = oxc
        
        #map oxdna coordinates to euler coordinates
        p = oxc.center - OX_TZU_DISP
        """
        A = np.zeros((3,3), dtype=float)
        
        #Euler angles
        #there's a minus sign because we want the inverse rotation (from rigid body to lab)
        alpha = -math.atan2(oxc.normal[0],-oxc.normal[1])
        betha = -math.acos(oxc.normal[2])
        gamma = -math.atan2(oxc.base_norv[2], oxc.base_v[2])
        
        #adeg = alpha/math.pi*180
        #bdeg = betha/math.pi*180
        #gdeg = gamma/math.pi*180
        #print(adeg,bdeg,gdeg)

        #Rotation matrix
        A[0][0] = math.cos(gamma)*math.cos(alpha)-math.cos(betha)*math.sin(gamma)*math.sin(alpha)
        A[0][1] = -math.cos(gamma)*math.sin(alpha)-math.cos(betha)*math.cos(alpha)*math.sin(gamma)
        A[0][2] = math.sin(gamma)*math.sin(betha)
        
        A[1][0] = math.cos(alpha)*math.sin(gamma)+math.cos(gamma)*math.cos(betha)*math.sin(alpha)
        A[1][1] = math.cos(gamma)*math.cos(betha)*math.cos(alpha)-math.sin(gamma)*math.sin(alpha)
        A[1][2] = -math.cos(gamma)*math.sin(betha)
        
        A[2][0] = math.sin(betha)*math.sin(alpha)
        A[2][1] = math.cos(alpha)*math.sin(betha)
        A[2][2] = math.cos(betha)
        
        self.frame.rot = A
        """
        ori = np.column_stack((oxc.base_norv,oxc.base_v,oxc.normal))
        self.frame = eframe(p,ori)



###########################
#Inverse (m1 stands for -1) of the Cayley transformation which maps a vector to a rotation (SO(3))
#a = 1 => variables are in radiants

def caym1(A,a) :
    c = 2*a/(1+np.trace(A))
    v = np.zeros(3,dtype=float)
    M = A -A.transpose()
    
    #note: v is vect(A-A.transposed()) 
    #vect(M), with M skew is defined as v = (M(2,1), M(0,2), M(1,0))
    #see, e.g., Daiva Petkevičiūtė thesis (Maddocks student)
    
    v[0] = M[2][1]
    v[1] = M[0][2]
    v[2] = M[1][0]
    t = c*v
    return t      

#############################
#base pair class. Stores two bases, a base pair frame (which is a eframe object) 
#and the intra coordinates
class base_pair :
    def __init__(self,b1,b2) :
        self.base_W = b1 #+
        self.base_C = b2 #-
      
        #flip the Crick base
        F = np.zeros((3,3), dtype=float)
        F[0][0] = 1.
        F[1][1] = -1.
        F[2][2] = -1.
        
        #compute average bp frame
        p = (b1.frame.pos+b2.frame.pos)*0.5
        DC = np.dot(b2.frame.orientation,F)
        A2 = np.dot(DC.transpose(),b1.frame.orientation)
        A = linalg.sqrtm(A2).real
        ori = np.dot(DC,A)
        self.frame = eframe(p,ori)
        
        #compute intra coordinates
        rot = caym1(A2,1)
        tr = np.dot(self.frame.orientation.transpose(),b1.frame.pos-b2.frame.pos)
        self.intra_coord = int_coord(tr, rot)

#############################
#junction class. Stores two base_pairs, a junction frame, and the inter coordinates (inter_coord)
class junction :
    def __init__(self,bp1,bp2) :
        self.base_pair1 = bp1 #bp n
        self.base_pair2 = bp2 #bp n+1

        #compute average junction frame        
        p = (bp1.frame.pos+bp2.frame.pos)*0.5
        A2 = np.dot(bp1.frame.orientation.transpose(),bp2.frame.orientation)
        A = linalg.sqrtm(A2).real
        ori = np.dot(bp1.frame.orientation,A)
        self.frame = eframe(p,ori)
        
        #compute inter coordinates

        rot = caym1(A2,1)
        tr = np.dot(self.frame.orientation.transpose(),bp2.frame.pos-bp1.frame.pos)
        self.inter_coord = int_coord(tr, rot)

##############################
#read oxdna trajectory (XXXTODO: pandas?)
#XXXTODO: read topology from topo file!
#XXXTODO: memory wise it's better to read AND print one snapshot at a time
def read_oxdna_trajectory(tr_file, Nb):
    #read Nb from topology!
    trajectory = []
    for line in tr_file.readlines():
        a = line.strip()[0]
        if a == 't':
            counts = 0
            #print(a)
            config = []
        else:
            if a != 'b' and a != 'E':
                counts += 1
                #print(counts)
                vals = line.split()
                #print(vals[0])
                c = np.array([float(vals[0]),float(vals[1]), float(vals[2])])
                #print(c)
                bv = np.array([float(vals[3]),float(vals[4]), float(vals[5])])
                #print(bv)
                n = np.array([float(vals[6]),float(vals[7]), float(vals[8])])
                #print(n)
                b = base(oxdna_frame(c, bv, n))
                config.append(b)
                if len(config) == Nb :
                    bps = []
                    juns = []
                    for k in range(0,int(Nb/2)) : #modify to account for general topology
                        bp =  base_pair(config[k],config[Nb-1-k])
                        bps.append(bp)
                    for k in range(0,len(bps)-1) :
                        jun = junction(bps[k],bps[k+1])
                        juns.append(jun)
                    trajectory.append(juns)
            
    return trajectory



#TEST
"""
c = np.array([8.51176430302203, 2.97670929267179, 26.2153512548549])
bv = np.array([0.38485191698477, 0.710495480645037, -0.589139350218692])
bnv = np.array([0.915651863281132, -0.37416748949958, 0.146902535959532])

rot = np.zeros((3,3), dtype=float)

oc = oxdna_frame(c,bv,bnv)

print("normal:")
print(oc.normal)
ec = eframe(c,rot)

#print(ec.tran)

print("rot matrix:")
print(ec.rot)


#the following should be 3x3 identity matrix
print("test: ")
print(np.dot(ec.rot,oc.base_norv))
print(np.dot(ec.rot,oc.base_v))
print(np.dot(ec.rot,oc.normal))
"""
iname = 'trajectory.dat'

ifile = open(iname,'r')

tr_data = read_oxdna_trajectory(ifile, 34)

print("reading test:")
#print(tr_data[0][0][0].center)
#print(tr_data[0][0][1].center)
print("JUNCTION 3, SNAP 0:")
print("###############\nNOTE1: translations are in oxdna units, rotations in radiants\nNOTE2: translations are between base/base pair barycenters\n##############")
print("BP1 INTRA: ")
print("translations:")
print(tr_data[0][3].base_pair1.intra_coord.tran)
print("rotations:")
print(tr_data[0][3].base_pair1.intra_coord.rot)
print("BP2 INTRA: ")
print("translations:")
print(tr_data[0][3].base_pair2.intra_coord.tran)
print("rotations:")
print(tr_data[0][3].base_pair2.intra_coord.rot)
print("JUNCTION INTER: ")
print("translations:")
print(tr_data[0][3].inter_coord.tran)
print("rotations:")
print(tr_data[0][3].inter_coord.rot)