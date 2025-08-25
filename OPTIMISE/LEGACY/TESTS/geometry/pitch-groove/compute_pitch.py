import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import CubicSpline
from scipy.signal import argrelextrema
#from oxdna_to_pdb import read_oxdna_trajectory as mapping_read_trj
#from oxdna_to_pdb import oxdna_nucleo_to_pdb

#from oxdna_to_internal_triadII import read_oxdna_trajectory_standard_order #Carlon style (see 2017 J. chem. phys)
from oxdna_to_internal_cgna import read_oxdna_trajectory_standard_order #Maddocks style (see theses)


if len(sys.argv) != 4:
    print("Unknown argument format")
    print("Usage: "+sys.argv[0] + " topo_file trajectory insnap")
    exit(1)

topo_file_name = sys.argv[1]
trj_file_name = sys.argv[2]
in_snap = int(sys.argv[3]) #ignore first in_snap snapshots (equilibration)


#things from pdb converter

P_proj_y = 9
P_y = P_proj_y*0.951056 #cos 18
P_x = P_proj_y*0.309016 #sin 18
offset_z = -1.4 #This is picked by hand!
hb_flag = False

Adenine = {
    'P': np.array([0.-P_x, 0.+P_y, 0.]),
    'Pideal': np.array([0., 0.+P_proj_y, offset_z]),
    "C1'": np.array([-2.479, 5.346, 0.]),
    'N9': np.array([-1.291, 4.498, 0.]),
    'C8': np.array([0.024, 4.897, 0.]),
    'N7': np.array([0.877, 3.902, 0.]),
    'C5': np.array([0.071, 2.771, 0.]),
    'C6': np.array([0.369, 1.398, 0.]),
    'N1': np.array([-0.668, 0.532, 0.]),
    'C2': np.array([-1.912, 1.023, 0.]),
    'N3': np.array([-2.320, 2.290, 0.]),
    'C4': np.array([-1.267, 3.124, 0.]),

}

if hb_flag :
    Adenine['N6'] = np.array([1.611, 0.909, 0.])
    
Guanine = {
    'P': np.array([0.-P_x, 0.+P_y, 0.]),
    'Pideal': np.array([0., 0.+P_proj_y, offset_z]),
    "C1'": np.array([-2.477, 5.399, 0.]),
    'N9': np.array([-1.289, 4.551, 0.]),
    'C8': np.array([0.023, 4.962, 0.]),
    'N7': np.array([0.870, 3.969, 0.]),
    'C5': np.array([0.071, 2.883, 0.]),
    'C6': np.array([0.424, 1.460, 0.]),
    'N1': np.array([-0.700, 0.641, 0.]),
    'C2': np.array([-1.999, 1.087, 0.]),
    'N3': np.array([-2.342, 2.364, 0.001]),
    'C4': np.array([-1.265, 3.177, 0.]),
}

if hb_flag :
    Guanine['O6'] = np.array([1.554, 0.955, 0.])
    Guanine['N2'] = np.array([-2.949, 0.139, -0.001])
    

Cytosine = {
    'P': np.array([0.-P_x, 0.+P_y, 0.]),
    'Pideal': np.array([0., 0.+P_proj_y, offset_z]),
    "C1'": np.array([-2.477, 5.402, 0.]),
    'N1': np.array([-1.285, 4.542, 0.]),
    'C2': np.array([-1.472, 3.158, 0.]),
    'N3': np.array([-0.391, 2.344, 0.]),
    'C4': np.array([0.837, 2.868, 0.]),
    'C5': np.array([1.056, 4.275, 0.]),
    'C6': np.array([-0.023, 5.068, 0.]),
   
}

if hb_flag :
    Cytosine['O2'] = np.array([-2.628, 2.709, 0.])
    Cytosine['N4'] = np.array([1.875, 2.027, 0.])

Thymine = {

    'P': np.array([0.-P_x, 0.+P_y, 0.]),
    'Pideal': np.array([0., 0.+P_proj_y, offset_z]),
    "C1'": np.array([-2.481, 5.354, 0.]),
    'N1': np.array([-1.284, 4.500, 0.]),
    'C2': np.array([-1.462, 3.135, 0.]),
    'N3': np.array([-0.298, 2.407, 0.]),
    'C4': np.array([0.994, 2.897, 0.]),
    'C5': np.array([1.106, 4.338, 0.]),
    #'C7': np.array([2.466, 4.961, 0.001]),
    'C6': np.array([-0.024, 5.057, 0.]),

}

if hb_flag :
    Thymine['O2'] = np.array([-2.562, 2.608, 0.])
    Thymine['O4'] = np.array([1.944, 2.119, 0.])

oxdna_cm_Tsukuba = {

    'A': np.array([0., 5.11, 0.]),
    'G': np.array([0., 5.11, 0.]),
    'C': np.array([0., 5.11, 0.]),
    'T': np.array([0., 5.11, 0.]),
}

OX_TO_ANG = 8.518

to_center = False
v_to_center = []

class oxdna_nucleotide:
    def __init__(self, id, c, bv, n, t):
        self.id = id
        self.center = c*OX_TO_ANG
        self.base_v = bv
        self.normal = n
        nrv = np.cross(n, bv)
        self.base_norv = nrv
        self.rot = np.column_stack((-nrv, -bv, -n))
        self.ty = t

# Covert oxdna_nucleotide (i.e. base vector, etc.), to pdb
def oxdna_nucleo_to_pdb(nucleo, cid):

    R = nucleo.rot

    # - nucleo.rot @ oxdna_cm_Tsukuba[nucleo.ty]
    translation = nucleo.center - R @ oxdna_cm_Tsukuba[nucleo.ty]
    pdb_coords = {}
    if nucleo.ty == 'A':
        for atom in Adenine:
            pdb_coords[atom] = R @ Adenine[atom] + translation

    elif nucleo.ty == 'C':
        for atom in Cytosine:
            pdb_coords[atom] = R @ Cytosine[atom] + translation

    elif nucleo.ty == 'G':
        for atom in Guanine:
            pdb_coords[atom] = R @ Guanine[atom] + translation

    elif nucleo.ty == 'T':
        for atom in Thymine:
            pdb_coords[atom] = R @ Thymine[atom] + translation
            
    if to_center:
        for atom in pdb_coords:
            pdb_coords[atom] -= v_to_center[cid]

    return pdb_coords

class topo:
    def __init__(self, nid, sid, bty, do, up):
        self.id = nid
        self.strand_id = sid
        self.base_type = bty
        self.down_id = do
        self.up_id = up

#read dna trajectory file (stand alone) and return oxdna_nucleotide trajectory and topology
def mapping_read_trj(tr_file, topo_file):

    global to_center
    global v_to_center

    Nb = 0
    Ns = 0
    nid = 0

    trajectory = []
    topology = {}
    counts = 0
    for line in topo_file.readlines():
        counts += 1
        vals = line.split()
        if counts == 1:
            Nb = int(vals[0])
            Ns = int(vals[1])
        else:
            to = topo(nid, int(vals[0]), vals[1], int(vals[2]), int(vals[3]))
            topology[nid] = to
            nid += 1
            
    counts = 0
    for line in tr_file.readlines():
        a = line.strip()[0]
        if a == 't':
            nid = 0
            # print(a)
            config = []
            counts += 1
            # print(counts)
            v_to_c = np.array([0.,0.,0.])
        else:
            if a != 'b' and a != 'E':
                vals = line.split()
                # print(vals[0])
                c = np.array([float(vals[0]), float(vals[1]), float(vals[2])])
                #print(c)
                bv = np.array([float(vals[3]), float(vals[4]), float(vals[5])])
                #print(bv)
                n = np.array([float(vals[6]), float(vals[7]), float(vals[8])])
                #print(n)
                b = oxdna_nucleotide(nid, c, bv, n, topology[nid].base_type)
                config.append(b)
                v_to_c += b.center
                if to_center == False :
                    if len("{:.{}f}".format(c[0], 2)) > 7 :
                        print("Warning! Coordinates are too large to be compatible with pdb format.")
                        print("Coordinates larger that 9.999999*10^4 lead to problems")
                        print("Triggering translation to center mode.")
                        print("Note that this might not be enough.")
                        print("You can rescale all coordinates, but you might lose precision.")
                        to_center = True
                nid += 1
                if len(config) == Nb:
                    trajectory.append(config)
                    v_to_center.append(v_to_c/Nb)

    return trajectory, topology

#Move phospate using backbone-backbone vector
def move_phosphate(pdb_c, pdb_c_p1_P) :
    pdb_c['P'] = pdb_c['P']/2 + pdb_c_p1_P/2



#read sequence
tofile = open(topo_file_name, 'r')

seq = ""

Njuns = 0
counts = 0
for line in tofile.readlines():
    vals = line.strip().split(" ")
    if counts == 0: Njuns = int(int(vals[0])/2)-1
    elif counts <= Njuns+1:
        seq+=vals[1]
    counts += 1

tofile.close()

print("Sequence: ", seq)

if Njuns-4 < 12:
    print("The sequence is too small. Can't compute pitch.")
    print("Remember that the first two and last two junctions are discarded")
    exit(1)


#ids of internal coordinates we use to compute pitch length. 8,11 = rise and twist
ids = [8,11]

in_j = 2 #ignore ends
fin_j = Njuns -2 #ignore ends

dimension = (Njuns-in_j-(Njuns-fin_j))*len(ids)

N_fit =  (Njuns-in_j-(Njuns-fin_j))

internal_coords = [] # internal coordinates

#average internal coordinates
mu_sampled = np.zeros(dimension, dtype = float) #average GS
mu_global_sampled = np.zeros(len(ids), dtype = float) #average int coordinates
cov_sampled = np.zeros((dimension,dimension), dtype = float)


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
    global cov_sampled
    
    Nsnaps = len(internal_coords)
    Ncoords = len(internal_coords[0])
    
    for i in range(Ncoords) :
        mu_sampled[i] = 0.
        for j in range(Ncoords) :
            cov_sampled[i][j] = 0.
    
    #sequence dependent
    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            mu_sampled[j] += internal_coords[i][j]/Nsnaps
            
    #everaged over sequence
    for i in range(Ncoords)        :
        mu_global_sampled[i%len(ids)]+=mu_sampled[i]/(Ncoords/len(ids))
    

    for i in range(Nsnaps) :
        for j in range(Ncoords) :
            for z in range(j,Ncoords) :
                cov_sampled[j][z] += (internal_coords[i][j] - mu_sampled[j])*(internal_coords[i][z] - mu_sampled[z])/Nsnaps
    
    for j in range(Ncoords) :
        for z in range(j+1,Ncoords) :
            cov_sampled[z][j] = cov_sampled[j][z]    
    
    return

def unscrumble_mu() :

    rise = []
    twist = []

    for i in range(len(mu_sampled)) :
        if i%2==0: rise.append(mu_sampled[i])
        else: twist.append(mu_sampled[i]) #we have to go from cgdna coords to degrees

    return rise, twist
    
    
def get_phs(trj,in_j,fin_j) :
    phs_W = []
    phs_C = []
    
    for i in range(len(trj)):
        phs_W1 = []
        phs_C1 = []   
        for j in range(1,int(len(trj[i])/2-1)):
            p1 = oxdna_nucleo_to_pdb(trj[i][j],0)['P']
            p2 = oxdna_nucleo_to_pdb(trj[i][j+1],0)['P']

            #phosphate is in beween backbones, more or less.
            p3 = p1/2+p2/2
            
            phs_W1.append(p3.tolist())
        for j in range(len(trj[i])-2,int(len(trj[i])/2),-1):
            p1 = oxdna_nucleo_to_pdb(trj[i][j],0)['P']
            p2 = oxdna_nucleo_to_pdb(trj[i][j+1],0)['P']

            #phosphate is in beween backbones, more or less.
            p3 = p1/2+p2/2
            
            phs_C1.append(p3)

        phs_W.append(phs_W1)
        #phs_C.append(np.flip(phs_C1))
        phs_C.append(phs_C1)
    return phs_W, phs_C

ifile = open(trj_file_name, 'r')
tfile = open(topo_file_name, 'r')

#Compute and store internal coordinates from sampled trajectories. Compute gs and cov

    
traj = read_oxdna_trajectory_standard_order(ifile, tfile)

ifile.close()
tfile.close()

ifile = open(trj_file_name, 'r')
tfile = open(topo_file_name, 'r')

tmp_trj , topo_mapping = mapping_read_trj(ifile,tfile)
#print(tmp_trj)
#print(len(tmp_trj))
trj_mapping = tmp_trj[in_snap+1:]
print(len(trj_mapping))

ifile.close()
tfile.close()

# True = overwrite, False = append
# append is for averaging over multiple trajectories

store_internal_coord(traj,ids,in_j,fin_j,in_snap,True)

#compute average coordinates
ave_stored()

#get phosphate positons, so we can get minor and major groove lengths
phs_W_tmp, phs_C_tmp = get_phs(trj_mapping,in_j,fin_j)

phs_W = np.array(phs_W_tmp)
phs_C = np.array(phs_C_tmp)

#print(phs_W[0])
#print(phs_C[0])

print(len(phs_W), len(phs_C))

dist_matrix = np.zeros((len(phs_C),len(phs_C)))

l1 = len(phs_W[0])
l2 = len(phs_C[0])
nbps_dmap = 10
grane = nbps_dmap*3 #20 points per base pair

pnt_W = phs_W[:,1:-1]

ofile = open("phs_W.dat", 'w')
for i in range(len(pnt_W[0])) :
    print(pnt_W[0,i,0], pnt_W[0,i,1], pnt_W[0,i,2], file = ofile)
ofile.close()

tng_W_tmp = []
for i in range(len(phs_W)):
    tng_W1 = []
    for j in range(1,l1-1):
        v = [0,0,0] 
        v[0] = phs_W[i][j+1][0]-phs_W[i][j-1][0]
        v[1] = phs_W[i][j+1][1]-phs_W[i][j-1][1]
        v[2] = phs_W[i][j+1][2]-phs_W[i][j-1][2]
        norm = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
        v[0] /= norm
        v[1] /= norm
        v[2] /= norm
        tng_W1.append(v) 

    tng_W_tmp.append(tng_W1)
    
tng_W = np.array(tng_W_tmp)
    
inter_p_W = []
for i in range(len(phs_W)):
    #interpolate phospate positions
    print(i, len(pnt_W[i]), len(tng_W[i]))
    t = np.arange(len(pnt_W[i]))
    spline_x = CubicSpline(t, pnt_W[i,:, 0], bc_type=((1, tng_W[i,0, 0]), (1, tng_W[i,-1, 0])))
    spline_y = CubicSpline(t, pnt_W[i,:, 1], bc_type=((1, tng_W[i,0, 1]), (1, tng_W[i,-1, 1])))
    spline_z = CubicSpline(t, pnt_W[i,:, 2], bc_type=((1, tng_W[i,0, 2]), (1, tng_W[i,-1, 2])))


    #compute distance map for 15 bps.
    #then, do running average
    
    ps1 = []

    for n in range(0,len(phs_W[0])-nbps_dmap) :
        t_fine = np.linspace(n, n+nbps_dmap, grane)
        x_interp = spline_x(t_fine)
        y_interp = spline_y(t_fine)
        z_interp = spline_z(t_fine)
        ps = np.column_stack((x_interp,y_interp,z_interp))
        ps1.append(ps)

    inter_p_W.append(ps1)
    
print("splined W")
    
pnt_C = phs_C[:,1:-1]

ofile = open("phs_C.dat", 'w')
for i in range(len(pnt_C[0])) :
    print(pnt_C[0,i,0], pnt_C[0,i,1], pnt_C[0,i,2], file = ofile)
ofile.close()

tng_C_tmp = []
for i in range(len(phs_C)):
    tng_C1 = []
    for j in range(1,l2-1): 
        v = [0,0,0] 
        v[0] = phs_C[i][j+1][0]-phs_C[i][j-1][0]
        v[1] = phs_C[i][j+1][1]-phs_C[i][j-1][1]
        v[2] = phs_C[i][j+1][2]-phs_C[i][j-1][2]
        norm = math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])
        v[0] /= norm
        v[1] /= norm
        v[2] /= norm
        tng_C1.append(v) 

    tng_C_tmp.append(tng_C1)
    
tng_C = np.array(tng_C_tmp)
inter_p_C = []
for i in range(len(phs_C)):
    print(i)
    t = np.arange(len(pnt_C[i]))
    spline_x = CubicSpline(t, pnt_C[i,:, 0], bc_type=((1, tng_C[i,0, 0]), (1, tng_C[i,-1, 0])))
    spline_y = CubicSpline(t, pnt_C[i,:, 1], bc_type=((1, tng_C[i,0, 1]), (1, tng_C[i,-1, 1])))
    spline_z = CubicSpline(t, pnt_C[i,:, 2], bc_type=((1, tng_C[i,0, 2]), (1, tng_C[i,-1, 2])))

    ps1 = []

    for n in range(0,len(phs_C[0])-nbps_dmap) :
        t_fine = np.linspace(n, n+nbps_dmap, grane)
        x_interp = spline_x(t_fine)
        y_interp = spline_y(t_fine)
        z_interp = spline_z(t_fine)
        ps = np.column_stack((x_interp,y_interp,z_interp))
        ps1.append(ps)

    inter_p_C.append(ps1)


def dist(a,b) :
    return math.sqrt( (b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]) + (b[2]-a[2])*(b[2]-a[2]))

dist_matrix = np.zeros((grane,grane))


print("splined C")

mins = []
majs = []
for n in range(len(inter_p_W[0])) :
    for i in range(grane) :
        for j in range(grane) :
             dist_matrix[i,j] = 0.
    for z in range(len(inter_p_W)) :
    #for n in range(5) :
        for i in range(grane) :
            for j in range(grane) :
                #print(z,n,i,j)
                dist_matrix[i,j] += dist(inter_p_W[z][n][i],inter_p_C[z][n][j])/len(inter_p_W)#/len(inter_p_W[z])
                
    AntiDiagonal = []

    for i in range(len(dist_matrix)):
        AntiDiagonal.append(dist_matrix[i][len(dist_matrix[i])-i-1])
    np_adiag = np.array(AntiDiagonal)
    local_minima = argrelextrema(np_adiag, np.less)
    print("Identified local minima [A] (subtract 5.8A to get groove widths (suggested by CURVES+))")
    grooves = np_adiag[local_minima]
    if len(grooves) ==2 :
        mins.append(grooves[0])
        majs.append(grooves[1])
    print(np_adiag[local_minima])

ofile = open("grooves.txt",'w')
for i in range(len(mins)): print(mins[i], majs[i], file=ofile)
ofile.close()


"""
print("Identified local minima [A] (subtract 5.8A to get groove widths (suggested by CURVES+))")
print(np_adiag[local_minima])

ofile = open("bb_distance_matrix.dat",'w')
for i in range(len(dist_matrix)):
    for j in range(len(dist_matrix[i])) :
        print(i,j,dist_matrix[i][j], file=ofile)

ofile.close()

print("printed dd matrix")


#searching for local minima of the anti diagonal. These are the groove widths
AntiDiagonal = []

for i in range(len(dist_matrix)):
    AntiDiagonal.append(dist_matrix[i][len(dist_matrix[i])-i-1])

np_adiag = np.array(AntiDiagonal)

local_minima = argrelextrema(np_adiag, np.less)

print("Identified local minima [A] (subtract 5.8A to get groove widths (suggested by CURVES+))")
print(np_adiag[local_minima])
"""

#pitch length

rise,twist = unscrumble_mu()

print(rise)
print(twist)

SD_pitch = []
SD_pitch_nbps = []

#print(len(mu_sampled))
#print(len(rise), len(twist))

#compute SD pitch
for i in range(len(rise)):
    pitch = 0
    angle = 0
    pitch_nbps = 0
    for j in range(i,len(twist)) :
        #print(i,j)
        #angle += twist[j]
        #pitch += rise[j]
        if angle+twist[j] >= 360 :
            pitch_nbps += (360-angle)/twist[j]
            pitch += rise[j]*((360-angle)/twist[j])
            SD_pitch.append(pitch)
            SD_pitch_nbps.append(pitch_nbps)
            break
        else:
            pitch_nbps += 1
            angle += twist[j]
            pitch += rise[j]
        
ave_pitch = 0
ave_pitch_nbps = 0
for i in range(len(SD_pitch)) : ave_pitch += SD_pitch[i]/len(SD_pitch)
for i in range(len(SD_pitch_nbps)) : ave_pitch_nbps += SD_pitch_nbps[i]/len(SD_pitch_nbps)

print("ave pitch [A]: " + str(ave_pitch))
print("ave pitch [nbps]: " + str(ave_pitch_nbps))


ofile = open("pitch.dat", 'w')

print("#id type pitch rise twist", file = ofile)
for i in range(len(SD_pitch)):
    print(str(i+2)+ " " + seq[i+2] + " " + str(SD_pitch_nbps[i])+" " + str(SD_pitch[i]) + " " + str(rise[i]) + " " + str(twist[i]),file=ofile)

ofile.close()

