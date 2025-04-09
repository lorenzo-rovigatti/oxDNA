# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 14:07:04 2022

@author: andrea bonato

code to map oxdna coordinates (center of mass + base vectors)
to internal coordinates (slide, twist, etc.)

it takes a oxdna trajectory and topology as input 
and outputs a file with the average internal coordinates


"""
import numpy as np
import math
from scipy import linalg
import sys



###############################
OX_TO_ANG = 8.518 #1 oxdna unit = 8.518 angstrongs
#Parameteres to map oxdna coordinates to interaction centers
#A = 0; C = 1, G = 2, T = 3
#Note: currently, only the COM is used to compute translations.
STACK_X = [0.405*OX_TO_ANG, 0.275*OX_TO_ANG, 0.405*OX_TO_ANG, 0.275*OX_TO_ANG]
HB_X = [0.465*OX_TO_ANG, 0.335*OX_TO_ANG, 0.465*OX_TO_ANG, 0.335*OX_TO_ANG]
BB_X = [-0.34*OX_TO_ANG, -0.34*OX_TO_ANG, -0.34*OX_TO_ANG, -0.34*OX_TO_ANG]
BB_Y = [-0.34*OX_TO_ANG, -0.34*OX_TO_ANG, -0.34*OX_TO_ANG, -0.34*OX_TO_ANG]

#HB_X = [0.405*OX_TO_ANG, 0.405*OX_TO_ANG, 0.405*OX_TO_ANG, 0.405*OX_TO_ANG]



###############################
#TZU_HB_PU = np.array([-0.684,0.5865,0.0])
#TZU_HB_PI = np.array([-0.3445,2.3755,0.0])
#TZU_HB_PU = np.array([-0.5143,0.5865,0.0])
#TZU_HB_PI = np.array([-0.5143,2.3755,0.0])
#TZU_HB_PU = np.array([-0.5143,0.795,0.0])
#TZU_HB_PI = np.array([-0.5143,2.5885,0.0])
TZU_HB_PU = np.array([0,0.795,0.0])
TZU_HB_PI = np.array([0,2.5885,0.0])

#TZU_HB_PU = np.array([0,0,0.0])
#TZU_HB_PI = np.array([0,0,0.0])

#parameter to map nucleotide centers to Euler translations
TZU_C_PU = np.array([0, 5.11, 0.0])
TZU_C_PI = np.array([0, 5.11, 0.0])

###############################
#class for oxdna coordinates
#c (center),bv (base-backbone vector),bnv (base normal vector)
#and normal (bv x bnv) are np arrays
#int_centers returns the positions of the interaction centers
class oxdna_frame :
    def __init__(self, c,bv,n):
        self.center = c*OX_TO_ANG
        self.base_v = bv
        self.normal = n
        self.base_norv = np.cross(n,bv)
    def int_centers(self,ty):
        tid = 0
        if ty == 'A' :
            tid = 0
        elif ty == 'C' :
            tid = 1
        elif ty == 'G' :
            tid = 2
        elif ty == "T" :
            tid = 3
        else :
            print("Unknown base type " + ty)
            quit()        
        
        back_bone = self.center+BB_X[tid]*self.base_v+BB_Y[tid]*self.base_norv
        stack = self.center+STACK_X[tid]*self.base_v
        hbond = self.center+HB_X[tid]*self.base_v
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
    def __init__(self,oxc,ty):
        self.oxframe = oxc
        self.type = ty
        
        
        #map oxdna coordinates to euler coordinates
        
        #some algebra to compute Tsukuba-like x,y vectors from oxdna x,y vectors
        #Note: this is not relevant anymore, since now Tsukuba y = -base_norv.
        #solution of equations:
            
        # ua+vb+wc = cos(alpha)
        # ud+ve+wf = 0
        # u^2+v^2+w^2 = 1
        #dot(new y,old x) > 0 for purines,
        #dot(new y,old x) < 0 for pirimidines
        
        #where new y vector = (u,v,w)
        #old y vector = (a,b,c)
        #z vector = (d,e,f)
        

        
        y = -oxc.base_v
        z = -oxc.normal
        x = -oxc.base_norv
        
        """
        
        cos_alpha = 0.99358 #alpha ~= 6.532°
        
        n = cos_alpha/(y[1]-z[1]*y[0]/z[0])
        m = -(y[2]-z[2]*y[0]/z[0])/(y[1]-z[1]*y[0]/z[0])
        
        p = -z[1]*n/z[0]
        q = -(z[2]+z[1]*m)/z[0]
        
        r = (q*q+m*m+1)
        s = (2*p*q+2*n*m)
        t = (p*p+n*n-1)
        
        
        ynew = np.zeros(3,dtype=float) 
        ynew_p = np.zeros(3,dtype=float) 
        ynew_m = np.zeros(3,dtype=float) 
        
        ynew_p[2] = (-s + math.sqrt(s*s-4*r*t))/(2*r)
        ynew_m[2] = (-s - math.sqrt(s*s-4*r*t))/(2*r)
        
        ynew_p[1] = n+m*ynew_p[2]
        ynew_m[1] = n+m*ynew_m[2]
        
        ynew_p[0] = p+q*ynew_p[2]
        ynew_m[0] = p+q*ynew_m[2]


        if np.dot(ynew_p,x) > 0 :
            if ty == 'A' or ty == 'G' :
                    ynew = ynew_p
            elif ty == 'C' or ty == 'T' :
                    ynew = ynew_m
        else:
            if ty == 'A' or ty == 'G' :
                ynew = ynew_m
            elif ty == 'C' or ty == 'T' :
                ynew = ynew_p
        """
 
        #ori = np.column_stack((np.cross(ynew,z),ynew,z))
        ori = np.column_stack((x,y,z))
        
        p = np.zeros(3,dtype=float)       

        #hb_pos = oxc.int_centers(ty)[1]       

      
        """  
        if ty == 'A' or ty == 'G' :
            p = hb_pos - np.dot(ori,TZU_HB_PU)
        elif ty == 'C' or ty == 'T' :
            p = hb_pos - np.dot(ori,TZU_HB_PI)
            #p = hb_pos - np.dot(ori,TZU_HB_PU)
  
        """
        #mapping to center of mass. Gives same result
       
        if ty == 'A' or ty == 'G' :
            p = oxc.center - np.dot(ori,TZU_C_PU)
        elif ty == 'C' or ty == 'T' :
            p = oxc.center - np.dot(ori,TZU_C_PI)
        
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
    
    v[0] = M[2][1].real
    v[1] = M[0][2].real
    v[2] = M[1][0].real
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
        #DC = np.dot(b1.frame.orientation,F)
        #A2 = np.dot(DC.transpose(),b2.frame.orientation)
        A = linalg.sqrtm(A2)
        ori = np.dot(DC,A)
        self.frame = eframe(p,ori)
        
        #compute intra coordinates
        rot = caym1(A2,1)*180./math.pi
        tr = np.dot(self.frame.orientation.transpose(),b1.frame.pos-b2.frame.pos)
        self.intra_coord = int_coord(tr.real, rot.real)

#############################
#junction class. Stores two base_pairs, a junction frame, and the inter coordinates (inter_coord)
class junction :
    def __init__(self,bp1,bp2) :
        self.base_pair1 = bp1 #bp n
        self.base_pair2 = bp2 #bp n+1

        #compute average junction frame        
        p = (bp1.frame.pos+bp2.frame.pos)*0.5
        A2 = np.dot(bp1.frame.orientation.transpose(),bp2.frame.orientation)
        A = linalg.sqrtm(A2)
        ori = np.dot(bp1.frame.orientation,A)
        self.frame = eframe(p,ori)
        
        #compute inter coordinates
        rot = caym1(A2,1)*180./math.pi
        tr = np.dot(self.frame.orientation.transpose(),bp2.frame.pos-bp1.frame.pos)
        self.inter_coord = int_coord(tr.real, rot.real)

##############################
#read oxdna trajectory
#XXXNOTE: There is no information on base pairs in topology file + 
#bp can potentially change runtime
#for now I'm assuming standard order: (A)Nbp-1,Nbp-2,...,0(B)0,1,..,Nbp-1 (Nbp = Nb/2 number of base pairs)  
#XXXTODO: memory wise it's better to read AND print one snapshot at a time

class topo :
    def __init__(self, nid, sid, bty, do, up) :
        self.id = nid
        self.strand_id = sid
        self.base_type = bty
        self.down_id = do
        self.up_id = up

def read_oxdna_trajectory_standard_order(tr_file, topo_file):
    
    Nb = 0
    Ns = 0
    nid = 0
    
    topology = []
    counts = 0
    for line in topo_file.readlines() :
        counts+=1
        vals = line.split()
        if len(vals) == 2 :
            Nb = int(vals[0])
            Ns = int(vals[1])
            if Ns != 2 :
                print("Number of strands in not 2.")
                quit()
        else :
           to = topo(nid, int(vals[0]), vals[1], int(vals[2]), int(vals[3]))
           topology.append(to)
           nid += 1
           
    trajectory = []
    counts = 0
    for line in tr_file.readlines():
        a = line.strip()[0]
        if a == 't':
            nid = 0
            #print(a)
            config = []
            counts+=1
            print(counts)
            #if counts == 501:
            #   break
        else:
            if a != 'b' and a != 'E':
                vals = line.split()
                #print(vals[0])
                c = np.array([float(vals[0]),float(vals[1]), float(vals[2])])
                #print(c)
                bv = np.array([float(vals[3]),float(vals[4]), float(vals[5])])
                #print(bv)
                n = np.array([float(vals[6]),float(vals[7]), float(vals[8])])
                #print(n)
                b = base(oxdna_frame(c, bv, n),topology[nid].base_type)
                config.append(b)
                nid += 1
                if len(config) == Nb :
                    bps = []
                    juns = []
                    for k in range(0,int(Nb/2)) : #modify to account for general topology
                        #bp =  base_pair(config[k],config[Nb-1-k])
                        bp =  base_pair(config[int(Nb/2)-1-k],config[int(Nb/2)+k])
                        bps.append(bp)
                    for k in range(0,len(bps)-1) :
                        jun = junction(bps[k],bps[k+1])
                        juns.append(jun)
                    trajectory.append(juns)
            
        if counts > 1000 :
            break
            
    return trajectory

def average_and_cov_internal_coord_over_trajectory(traj) :
    Nsn = len(traj)
    Nj = len(traj[0])
    
    Nbp = Nj+1
    av_intra_tr = np.zeros((Nbp,3),dtype = float)
    av_intra_rot = np.zeros((Nbp,3),dtype = float)
    av_inter_tr = np.zeros((Nj,3),dtype = float)
    av_inter_rot = np.zeros((Nj,3),dtype = float)
    unified_av = np.zeros((Nj,12),dtype = float)
    unified_coord1 = np.zeros(12,dtype = float)
    unified_coord2 = np.zeros(12,dtype = float)
    cov = np.zeros((Nj*12,Nj*12),dtype=float)
    topo = []
    
    for i in range (0, Nj) :
        topo.append(traj[0][i].base_pair1.base_W.type)
        if i == Nj-1 :
            topo.append(traj[0][i].base_pair1.base_C.type)       
    
    for i in range (0,Nsn) :
        for j in range (0, Nj) :
            av_intra_tr[j] += traj[i][j].base_pair1.intra_coord.tran/(float(Nsn))
            av_intra_rot[j] += traj[i][j].base_pair1.intra_coord.rot/(float(Nsn))
            av_inter_tr[j] += traj[i][j].inter_coord.tran/(float(Nsn))
            av_inter_rot[j] += traj[i][j].inter_coord.rot/(float(Nsn))
            if j == Nj-1 :
                av_intra_tr[j+1] += traj[i][j].base_pair2.intra_coord.tran/(float(Nsn))
                av_intra_rot[j+1] += traj[i][j].base_pair2.intra_coord.rot/(float(Nsn))
                
    for i in range(0,Nj):
        unified_av[i,0] = av_intra_tr[i,0]
        unified_av[i,1] = av_intra_tr[i,1]
        unified_av[i,2] = av_intra_tr[i,2]
        unified_av[i,3] = av_intra_rot[i,0]*(math.pi/180.*5.)
        unified_av[i,4] = av_intra_rot[i,1]*(math.pi/180.*5.)
        unified_av[i,5] = av_intra_rot[i,2]*(math.pi/180.*5.)
        unified_av[i,6] = av_inter_tr[i,0]
        unified_av[i,7] = av_inter_tr[i,1]
        unified_av[i,8] = av_inter_tr[i,2]
        unified_av[i,9] = av_inter_rot[i,0]*(math.pi/180.*5.)
        unified_av[i,10] = av_inter_rot[i,1]*(math.pi/180.*5.)
        unified_av[i,11] = av_inter_rot[i,2]*(math.pi/180.*5.)
        
    print(unified_av)
                
    for i in range(0,Nsn) :      
        print("cov-"+str(i))
        for j in range(1,Nj-1) :
            unified_coord1[0] = traj[i][j].base_pair1.intra_coord.tran[0]
            unified_coord1[1] = traj[i][j].base_pair1.intra_coord.tran[1]
            unified_coord1[2] = traj[i][j].base_pair1.intra_coord.tran[2]
            unified_coord1[3] = traj[i][j].base_pair1.intra_coord.rot[0]*(math.pi/180.*5.)
            unified_coord1[4] = traj[i][j].base_pair1.intra_coord.rot[1]*(math.pi/180.*5.)
            unified_coord1[5] = traj[i][j].base_pair1.intra_coord.rot[2]*(math.pi/180.*5.)
            unified_coord1[6] = traj[i][j].inter_coord.tran[0]
            unified_coord1[7] = traj[i][j].inter_coord.tran[1]
            unified_coord1[8] = traj[i][j].inter_coord.tran[2]
            unified_coord1[9] = traj[i][j].inter_coord.rot[0]*(math.pi/180.*5.)
            unified_coord1[10] = traj[i][j].inter_coord.rot[1]*(math.pi/180.*5.)
            unified_coord1[11] = traj[i][j].inter_coord.rot[2]*(math.pi/180.*5.)
            
            for z in range(j,Nj-1) :
                unified_coord2[0] = traj[i][z].base_pair1.intra_coord.tran[0]
                unified_coord2[1] = traj[i][z].base_pair1.intra_coord.tran[1]
                unified_coord2[2] = traj[i][z].base_pair1.intra_coord.tran[2]
                unified_coord2[3] = traj[i][z].base_pair1.intra_coord.rot[0]*(math.pi/180.*5.)
                unified_coord2[4] = traj[i][z].base_pair1.intra_coord.rot[1]*(math.pi/180.*5.)
                unified_coord2[5] = traj[i][z].base_pair1.intra_coord.rot[2]*(math.pi/180.*5.)
                unified_coord2[6] = traj[i][z].inter_coord.tran[0]
                unified_coord2[7] = traj[i][z].inter_coord.tran[1]
                unified_coord2[8] = traj[i][z].inter_coord.tran[2]
                unified_coord2[9] = traj[i][z].inter_coord.rot[0]*(math.pi/180.*5.)
                unified_coord2[10] = traj[i][z].inter_coord.rot[1]*(math.pi/180.*5.)
                unified_coord2[11] = traj[i][z].inter_coord.rot[2]*(math.pi/180.*5.)
                
                for m in range(12):
                    for n in range(12):
                        cov[j*12+m,z*12+n] += (unified_coord1[m]-unified_av[j,m])*(unified_coord2[n]-unified_av[z,n])/(1.*Nsn)
                        
    for i in range(len(cov)) :
        for j in range(i+1,len(cov)):
            cov[j,i] = cov[i,j]
    
    return(topo, av_intra_tr, av_intra_rot, av_inter_tr, av_inter_rot, cov)



#read trajectory and topology (argv[1], argv[2])

"""
if len(sys.argv) != 5 :
    print("Usage: " + sys.argv[0]+" trajectory_file topology_file nreps snap0")
    sys.exit()


iname = sys.argv[1]
tname = sys.argv[2]
nreps = int(sys.argv[3])
snap0 = int(sys.argv[4])
"""

iname = 'trajectory.dat'
tname = 'generated.top'
snap0 = 100

    
ifile = open(iname,'r')
tfile = open(tname,'r')

tr_data = read_oxdna_trajectory_standard_order(ifile, tfile)

Nj = len(tr_data[0])

ifile.close()
tfile.close()

    
to, av_intra_tr, av_intra_rot, av_inter_tr, av_inter_rot, cov = average_and_cov_internal_coord_over_trajectory(tr_data)

oname = 'av_int_coord_mapv4.txt'
ofile = open(oname, 'w')

print('#btype shear stretch stagger buckle propeller opening shift slide rise tilt roll twist', file=ofile)
for i in range (0,len(to)) :
    line = to[i] + " " + str(av_intra_tr[i][0]) + " " + str(av_intra_tr[i][1]) + " " + str(av_intra_tr[i][2])
    line = line + " " + str(av_intra_rot[i][0]) + " " + str(av_intra_rot[i][1]) + " " + str(av_intra_rot[i][2])
    if i < len(to)-1 :
        line = line + " " + str(av_inter_tr[i][0]) + " " + str(av_inter_tr[i][1]) + " " + str(av_inter_tr[i][2])
        line = line + " " + str(av_inter_rot[i][0]) + " " + str(av_inter_rot[i][1]) + " " + str(av_inter_rot[i][2])
    print(line,file=ofile)
    
ofile.close()

oname = 'cov_mapv4.txt'
ofile = open(oname, 'w')

print('#shear stretch stagger buckle propeller opening shift slide rise tilt roll twist ...', file=ofile)
for i in range(len(cov)) :
    line = str(cov[i,0])
    for j in range(1,len(cov)) :
        line = line + " " + str(cov[i,j])        
    print(line,file=ofile)
    
ofile.close()



reduced_cov = np.zeros((Nj*5,Nj*5),dtype=float)

for i in range(10) :
    for j in range(len(reduced_cov)) :
        reduced_cov[i,j] = 0
        reduced_cov[j,i] = 0
        
for i in range(24) :
    for j in range(len(cov)) :
        cov[i,j] = 0
        cov[j,i] = 0

nrow = -1
for i in range(len(cov)) :
    if (i%12 == 4 or i%12 == 8 or i%12 == 9 or i%12 == 10 or i%12 == 11):
        nrow += 1
        ncol = -1
        for j in range(1,len(cov)) :
                if (j%12 == 4 or j%12 == 8 or j%12 == 9 or j%12 == 10 or j%12 == 11):
                    ncol+=1
                    reduced_cov[nrow,ncol] = cov[i,j]     



reduced_cov_nozeros = np.zeros(((Nj-4)*5,(Nj-4)*5),dtype=float)

for i in range(10,len(reduced_cov)-10) :
    for j in range(10,len(reduced_cov)-10) :
        reduced_cov_nozeros[i-10,j-10] = reduced_cov[i,j]
        
oname = 'reduced_cov_mapv4.txt'
ofile = open(oname, 'w')

print('#propeller rise tilt roll twist ...', file=ofile)
for i in range(len(reduced_cov_nozeros)) :
        line = str(reduced_cov_nozeros[i,0])
        for j in range(1,len(reduced_cov_nozeros)) :
                    line = line + " " + str(reduced_cov_nozeros[i,j])        
        print(line,file=ofile)
    
ofile.close()


inv_red_cov = np.linalg.inv(reduced_cov_nozeros)