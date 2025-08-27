#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 18:26:39 2025

@author: Andrea Bonato

This code converts any given oxdna trajectory (output from stand alone) to a
 simplified pdb trajectory. The conversion is based on the nucleotide to Tzukuba ideal nucleotides mapping.
 The atoms included in the pdb are the rings of the bases, the phosphates in the backbone (just P, not the whole group), and the atom of the sugars bonded to the bases.
 Optionally, the atoms involved in the hb bonds can be added, these are handy to distinguish between AT and CG.
 It works for both oxdna2 and oxdna3.

"""

import numpy as np
import sys
import argparse
import copy

print("USAGE")
print("Usage: python3 "+sys.argv[0]+ " trajectory_file topology_file")
print("*********************************")
print("OPTIONS:")
print("-hy: include atoms involved in hydrogen bonds")
print("-c: move all configurations to the center of the simulation box")
print("*********************************")
print("NOTES")
print("The output is called trajectory.pdb.")
print("For each nucleotide, the output will contain about 10 atoms, making it significantly larger in size.")
print("Due to the nature of the pdb format, there are some limitations: ")
print("Number of strands per configuration < 1156.")
print("Total number of nucleotides per configuration < 10^4.")
print("Coordinates < 10^5; Precision is 0.01.")
print("If the coordinates are too large, the configurations will be automatically moved to the center of the simulation box.")
print("*********************************")

if len(sys.argv) < 3 :
    print("Can't decode arguments.")
    exit(1)

tr_name = sys.argv[1]
topo_name = sys.argv[2]

hb_flag = False #True if we include hb atoms
to_center = False #True if we ant to move all configurations to the center of the box
v_to_center = []

if len(sys.argv) >= 3:
    for n in range(3,len(sys.argv)):
        if sys.argv[n] == "-hy" : hb_flag = True
        if sys.argv[n] == "-c" : to_center = True

#Note: if the nucleotide is bonded, we use P and the backbone-backbone vector
#to place the phospate, otherwise we use P ideal, and no backbone-backbone.
#We have a well defined mapping just for the bases.
#The position of the phospate is a reconstructed.
#The average projection of P on the ideal base plane in about [0,9,0]. To place P between
#bonded nucleotides, we rotate its projection by half the average twist (average shift and slide are about 0) and then add half of the backbone-backbone vector

P_proj_y = 9
P_y = P_proj_y*0.951056 #cos 18
P_x = P_proj_y*0.309016 #sin 18
offset_z = -1.4 #This is picked by hand!

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


class oxdna_nucleotide:
    def __init__(self, id, c, bv, n, t):
        self.id = id
        nrv = np.cross(n, bv)
        self.center = c*OX_TO_ANG
        self.backbone = (c+0.3408*nrv-0.34*bv)*OX_TO_ANG #oxview style: for oxview style we need to know where the backbone site is
        self.base_v = bv
        self.normal = n
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
            
    pdb_coords['P'] = nucleo.backbone  ##oxview style: place P on the backbone site.
            
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
def read_oxdna_trajectory(tr_file, topo_file):

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
            # if counts == 501:
            #   break
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
    
#Convert integer sting id to chain identifier.
#Note that there can be at most 1156 diffenet strands!
#In pdb the chain id goes in coulms 21-22, and we can use numbers and letters (34 charachters).
#I won't use special characters (don't know if they are allowed)
def chain_number(sid) :

    c1 = sid%34
    c2 = (sid//34)%34
    string = ""
    if c1 < 10:
        string += str(c1)
    else:
        string += chr(65+c1)
        
    if sid >= 34 :
        if c2 < 10:
            string += str(c2)
        else:
            string += chr(65+c2)
    
    return(string)


#Given oxdna_nucleotide trajectory, prints pdb. Ofile = output file
def print_pdb_file(trajectory, topology, ofile):

    P_list = []

    #Convert trajectory to pdb
    print("HEADER    DNA STRUCTURE", file=ofile)
    count_m = 0
    for conf in trajectory:
        count_m += 1
        print("MODEL    "+str(count_m), file=ofile)
        count_n = 0
        count_at = 0
        P_list = []
        pdb_c =  oxdna_nucleo_to_pdb(conf[0],count_m-1)
        pdb_c_p1 = None
        for n in range(1,len(conf)+1):
            nucleo = conf[n-1]
            count_n += 1
            strand = chain_number(topology[nucleo.id].strand_id-1)
            if n < len(conf) :
                pdb_c_p1 = oxdna_nucleo_to_pdb(conf[n],count_m-1)
                #move_phosphate(pdb_c, pdb_c_p1['P'])  #oxview style: place backbone where it is, without moving it
            
            # print(pdb_c)
            P_list_1 = ["","",""]
            P_list_1[0] = strand
            counts_innuc = 0
            for atom in pdb_c:
                
                #print(counts_innuc)
                
                if (topology[nucleo.id].up_id != -1) and atom == "Pideal":
                    counts_innuc += 1
                    continue
                if (topology[nucleo.id].up_id == -1) :
                    if atom == "P":
                        counts_innuc += 1
                        continue
                
                count_at += 1
                counts_innuc += 1
                if atom == 'P' or atom == 'Pideal':
                    P_list_1[1]=str(count_at)
                if (nucleo.ty == 'A' or nucleo.ty == 'G') and atom == "C1'" :
                    P_list_1[2] = str(count_at)
                if (nucleo.ty == 'C' or nucleo.ty == 'T') and atom == "C1'" :
                    P_list_1[2] = str(count_at)
                if counts_innuc == len(pdb_c):
                    P_list.append(P_list_1)
                string = "ATOM"
                nspaces = 7 - len(str(count_at))
                for i in range(nspaces):
                    string += " "
                nspaces = 0
                if atom != 'Pideal':
                    string += str(count_at) + "  " + atom
                    nspaces = 5-len(atom)
                else:
                    string += str(count_at) + "  " + 'P'
                    nspaces = 5-len('P')
                for i in range(nspaces):
                    string += " "
                string += "D"+nucleo.ty
                if len(strand) == 1:
                    string += " " + strand
                else:
                    string += strand
                nspaces = 4 - len(str(count_n))
                for i in range(nspaces):
                    string += " "
                string += str(count_n) + "    "
                nspaces_x = 8 - len("{:.{}f}".format(pdb_c[atom][0], 2))
                nspaces_y = 8 - len("{:.{}f}".format(pdb_c[atom][1], 2))
                nspaces_z = 8 - len("{:.{}f}".format(pdb_c[atom][2], 2))
                for i in range(nspaces_x) :
                    string += " "
                string += "{:.{}f}".format(pdb_c[atom][0], 2)
                for i in range(nspaces_y) :
                    string += " "
                string += "{:.{}f}".format(pdb_c[atom][1], 2)
                for i in range(nspaces_z) :
                    string += " "
                string += "{:.{}f}".format(pdb_c[atom][2], 2) + "  1.00  1.00           " + atom[0]
                print(string, file=ofile)
                
            pdb_c = pdb_c_p1
            
            
        #no triangles style: to remove the triangles we remove one bond per P atom.    
        for i in range(0, len(P_list)):
            string = "CONECT"
            nspaces = 5 - len(P_list[i][1])
            for z in range(nspaces):
                string += " "
            string += P_list[i][1]
            nspaces = 5 - len(P_list[i][2])
            for z in range(nspaces):
                string += " "
            string += P_list[i][2]
            print(string, file=ofile)
            
            if i < len(P_list) -1 and P_list[i][0] == P_list[i+1][0]:
                string = "CONECT"
                nspaces = 5 - len(P_list[i][1])
                for z in range(nspaces):
                    string += " "
                string += P_list[i][1]
                nspaces = 5 - len(P_list[i+1][1])
                for z in range(nspaces):
                    string += " "
                string += P_list[i+1][1]
                print(string, file=ofile)
            
        print("TER", file=ofile)
        print("ENDMDL", file=ofile)

    return

hb_pos_ox3 = {
	'A' : 0.43,
	'G' : 0.43,
	'C' : 0.37,
	'T' : 0.37,
}

#This function prints bary, backbone and hydr of oxdna nucleotides in lammps format.
#It is useful to compare nucleotide geometry to Tsukuba mapping
#Can do both oxdna2 and oxdna3 nucleotides
def print_xyz(trajectory, topology, ofile, ox_v = "ox3"):

    count_m = 0
    for conf in trajectory:
        count_m += 1
        print("ITEM: TIMESTEP", file = ofile)
        print(str(count_m), file=ofile)
        print( "ITEM: NUMBER OF ATOMS ", file = ofile )
        print(str(len(conf)*3),file=ofile)
        print("ITEM: BOX BOUNDS pp pp pp", file = ofile)
        print("-3.0000000000000000e+02 3.0000000000000000e+02", file = ofile)
        print("-3.0000000000000000e+02 3.0000000000000000e+02", file = ofile)
        print("-3.0000000000000000e+02 3.0000000000000000e+02", file = ofile)
        print("ITEM: ATOMS id mol type x y z ix iy iz", file = ofile)
        count_n = 0
        count_at = 0
        
        for n in range(0,len(conf)):
            count_n += 1
            c = copy.deepcopy(conf[n].center)
            if to_center : 
            	c -= v_to_center[count_m-1]
            string = str(count_n) + " 1 1 " + str(c[0]) + " " + str(c[1]) + " " + str(c[2]) + " 0 0 0"
            print(string, file = ofile)
            count_n += 1
            if ox_v == "ox2" : hb = conf[n].center + 0.4*OX_TO_ANG*conf[n].base_v
            else: hb = conf[n].center + hb_pos_ox3[conf[n].ty]*OX_TO_ANG*conf[n].base_v
            if to_center : hb -= v_to_center[count_m-1]
            string = str(count_n) + " 1 2 " + str(hb[0]) + " " + str(hb[1]) + " " + str(hb[2]) + " 0 0 0"
            print(string, file = ofile)
            count_n += 1
            bb = conf[n].center - 0.34*OX_TO_ANG*conf[n].base_v + 0.34*OX_TO_ANG*conf[n].base_norv
            if to_center : bb -= v_to_center[count_m-1] 
            string = str(count_n) + " 1 3 " + str(bb[0]) + " " + str(bb[1]) + " " + str(bb[2]) + " 0 0 0"
            print(string, file = ofile)
    return




cfile = open(tr_name, 'r')
tfile = open(topo_name, 'r')

traj, topo = read_oxdna_trajectory(cfile, tfile)

cfile.close()
tfile.close()

ofile = open("trajectory.pdb", 'w')

print_pdb_file(traj,topo, ofile)

ofile.close()

#ofile = open("trajectory.lammpstrj", 'w')

#print_xyz(traj,topo, ofile)

#ofile.close()
