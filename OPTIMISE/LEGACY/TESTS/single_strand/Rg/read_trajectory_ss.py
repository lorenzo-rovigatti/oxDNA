import numpy as np
import math
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

        y = -oxc.base_v
        z = -oxc.normal
        x = -oxc.base_norv
 
        ori = np.column_stack((x,y,z))
        
        p = np.zeros(3,dtype=float)       

        #hb_pos = oxc.int_centers(ty)[1]       

        #mapping to center of mass. Gives same result
       
        if ty == 'A' or ty == 'G' :
            p = oxc.center - np.dot(ori,TZU_C_PU)
        elif ty == 'C' or ty == 'T' :
            p = oxc.center - np.dot(ori,TZU_C_PI)
        
        self.frame = eframe(p,ori)
        
        
#junction class. Stores two base_pairs, a junction frame, and the inter coordinates (inter_coord)
class junction :
    def __init__(self,bp1,bp2) :
        self.base1 = bp1 #bp n
        self.base2 = bp2 #bp n+1

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



#returns list of list of base objects
def read_oxdna_trajectory_standard_order(tr_file, topo_file):
    
    Nb = 0
    Ns = 0
    nid = 0
    
    topology = []
    counts = 0
    for line in topo_file.readlines() :
        counts+=1
        vals = line.split()
        if len(vals) == 1 :
            Nb = int(vals[0])
            Ns = int(vals[1])
            if Ns != 1 :
                print("Number of strands in not 1.")
                exit(1)
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
                    for k in range(0,Nb-1) : #here bp1 and bp2 of junction will be base1 and base2
                        jun = junction(config[k],config[k+1])
                        juns.append(jun)
                    trajectory.append(juns)
            
    return trajectory
