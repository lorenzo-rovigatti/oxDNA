#!/usr/bin/env python

import sys
import itertools
import numpy as np
from math import sqrt
import copy

import base
import utils
from readers import LorenzoReader

BASE_SHIFT = 1.13
COM_SHIFT = 0.5

class Nucleotide(object):
    serial_residue = 1
    def __init__(self, name, idx):
        object.__init__(self)
        self.name = name
        self.base = name[1:]
        self.idx = idx
        self.base_atoms = []
        self.phosphate_atoms = []
        self.sugar_atoms = []
        self.named_atoms = {}
        self.ring_names = ["C2", "C4", "C5", "C6", "N1", "N3"]
        self.chain_id = None

    def get_atoms(self):
        return self.base_atoms + self.phosphate_atoms + self.sugar_atoms

    def add_atom(self, a):
        if 'P' in a.name or a.name == "HO5'": self.phosphate_atoms.append(a)
        elif "'" in a.name: self.sugar_atoms.append(a)
        else: self.base_atoms.append(a)

        self.named_atoms[a.name] = a
        if self.chain_id == None: self.chain_id = a.chain_id

    def get_com(self, atoms=None):
        if atoms == None: atoms = self.atoms
        com = np.array([0., 0., 0.])
        for a in atoms:
            com += a.pos

        return com/len(atoms)

    def compute_a3(self):
        base_com = self.get_com(self.base_atoms)
        parallel_to = self.get_com(self.phosphate_atoms) - base_com
        self.a3 = np.array([0., 0., 0.])
        
        for perm in itertools.permutations(self.ring_names, 3):
            p = self.named_atoms[perm[0]]
            q = self.named_atoms[perm[1]]
            r = self.named_atoms[perm[2]]
            v1 = p.pos - q.pos
            v2 = p.pos - r.pos
            v1 /= sqrt(np.dot(v1, v1))
            v2 /= sqrt(np.dot(v2, v2))
            if abs(np.dot(v1, v2)) > 0.01 or 1:
                a3 = np.cross(v1, v2)
                a3 /= sqrt(np.dot(a3, a3))
                if np.dot(a3, parallel_to) < 0.: a3 = -a3
                self.a3 += a3

        self.a3 /= sqrt(np.dot(self.a3, self.a3))

    def compute_a1(self):
        if "C" in self.name or "T" in self.name:
            pairs = [ ["N3", "C6"], ["C2", "N1"], ["C4", "C5"] ]
        else:
            pairs = [ ["N1", "C4"], ["C2", "N3"], ["C6", "C5"] ]

        self.a1 = np.array([0., 0., 0.])
        for pair in pairs:
            p = self.named_atoms[pair[0]]
            q = self.named_atoms[pair[1]]
            diff = p.pos - q.pos
            self.a1 += diff

        self.a1 /= sqrt(np.dot(self.a1, self.a1))

    def compute_as(self):
        self.compute_a1()
        self.compute_a3()
        self.a2 = np.cross(self.a3, self.a1)
        self.check = abs(np.dot(self.a1, self.a3))

    def to_pdb(self, chain_identifier):
        res = []
        for a in self.atoms:
            res.append(a.to_pdb(chain_identifier, Nucleotide.serial_residue))

        Nucleotide.serial_residue += 1
        return "\n".join(res)

    def to_mgl(self):
        res = []
        for a in self.atoms:
            res.append(a.to_mgl())

        return "\n".join(res)

    def rotate(self, R):
        com = self.get_com()
        for a in self.atoms:
            a.pos = np.dot(R, a.pos)

        self.compute_as()

    def set_com(self, new_com):
        com = self.get_com()
        for a in self.atoms:
            a.pos += new_com - com - COM_SHIFT*self.a1

    def set_base(self, new_base_com):
        com = self.get_com()
        atoms = [v for k, v in self.named_atoms.iteritems() if k in self.ring_names]
        ring_com = self.get_com(atoms)
        for a in self.atoms:
            a.pos += new_base_com - ring_com - BASE_SHIFT*self.a1

        self.compute_as()

    atoms = property(get_atoms)


class Atom(object):
    serial_atom = 1
    def __init__(self, pdb_line):
        object.__init__(self)
        # http://cupnet.net/pdb-format/
        self.name = pdb_line[12:16].strip()
        self.chain_id = pdb_line[21:22].strip()
        self.residue = pdb_line[17:20].strip()
        self.residue_idx = int(pdb_line[22:26])
        # convert to oxDNA's length unit
        self.pos = np.array([float(pdb_line[31:38]), float(pdb_line[38:46]), float(pdb_line[46:54])])# / 8.518

    def is_hydrogen(self):
        return "H" in self.name

    def to_pdb(self, chain_identifier, serial_residue):
        #res = "%-6s%5d %4s%1s%3s %1s%4d%1s   %11.3f%11.3f%11.3f%6.2f%6.2f          %2s%2s" % ("ATOM", Atom.serial_atom, self.name, " ", self.residue, chain_identifier, serial_residue, " ", self.pos[0], self.pos[1], self.pos[2], 1.00, 0.00, " ", " ")
        res = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%-4s%-2s%-2s" % ("ATOM", Atom.serial_atom, self.name, " ", self.residue, chain_identifier, serial_residue, " ", self.pos[0], self.pos[1], self.pos[2], 1.00, 0.00, " ", " ", " ")
        Atom.serial_atom += 1
        return res

    def to_mgl(self):
        colors = {"C" : "0,1,1", "P" : "1,1,0", "O" : "1,0,0", "H" : "0.5,0.5,0.5", "N" : "0,0,1"}
        for c in colors:
            if c in self.name: color = colors[c]
        r = 0.5
        return "%s %s %s @ %f C[%s]" % (self.pos[0], self.pos[1], self.pos[2], r, color)


if len(sys.argv) < 4:
    print >> sys.stderr, "Usage is %s dd_file configuration topology" % sys.argv[0]
    sys.exit(1)

with open(sys.argv[1]) as f:
    nucleotides = []
    old_residue = ""
    for line in f.readlines():
        if len(line) > 77:
            na = Atom(line)
            if na.residue_idx != old_residue:
                nn = Nucleotide(na.residue, na.residue_idx)
                nucleotides.append(nn)
                old_residue = na.residue_idx
            nn.add_atom(na)
       
bases = {}
for n in nucleotides:
    n.compute_as()
    if n.base in bases:
        if n.check < bases[n.base].check: bases[n.base] = copy.deepcopy(n)
    else: bases[n.base] = n

for n in nucleotides:
    n.a1, n.a2, n.a3 = utils.get_orthonormalized_base(n.a1, n.a2, n.a3)

#for b in ("A", "T", "C", "G"):
#    print bases[b].to_pdb("A")

lr = LorenzoReader(sys.argv[2], sys.argv[3])
s = lr.get_system()

def align(full_base, ox_base):
    #print full_base.base, np.dot(full_base.a1, ox_base._a1), np.dot(full_base.a3, ox_base._a3)

    theta = utils.get_angle(full_base.a3, ox_base._a3)
    axis = np.cross(full_base.a3, ox_base._a3)
    axis /= sqrt(np.dot(axis, axis))
    #print np.dot(axis, ox_base._a1)
    R = utils.get_rotation_matrix(axis, theta)
    full_base.rotate(R)

    theta = utils.get_angle(full_base.a1, ox_base._a1)
    axis = np.cross(full_base.a1, ox_base._a1)
    axis /= sqrt(np.dot(axis, axis))
    R = utils.get_rotation_matrix(axis, theta)
    full_base.rotate(R)
    old_a1 = np.array(full_base.a1)

    #print " ", np.dot(full_base.a1, ox_base._a1), np.dot(full_base.a3, ox_base._a3)

#print ".Box:100,100,100"
ox_nucleotides = []
s.map_nucleotides_to_strands()

print "HEADER    t=1.12"
print "MODEL     1"
print "REMARK ## 0,0"
first_id = True
old_chain_id = -1
for nucleotide in s._nucleotides:
    nb = base.number_to_base[nucleotide._base]
    my_base = copy.deepcopy(bases[nb])
    my_base.chain_id = s._nucleotide_to_strand[nucleotide.index]
    if first_id:
        old_chain_id = my_base.chain_id
        first_id = False
    if my_base.chain_id != old_chain_id:
        print >> sys.stderr, "TERRING: ", old_chain_id, my_base.chain_id
        old_chain_id = my_base.chain_id
        print "TER"
    my_base.idx = (nucleotide.index % 12) + 1
    align(my_base, nucleotide)
    my_base.set_base(nucleotide.pos_base*8.518)
    ox_nucleotides.append(my_base)
    #my_base.set_com(nucleotide.cm_pos*8.518)
    print my_base.to_pdb("A")
    #print my_base.to_mgl()
print "REMARK ## 0,0"
print "TER"
print "ENDMDL"


# remove the next line to print the distance between P atoms
sys.exit(1)
def print_P_distance_dd(nucls, divide_by=1):
    if len(nucls) != 24:
        print >> sys.stderr, "%s can be used for dodecamers only" % print_P_distance_dd.__name__
        sys.exit(1)
    # check the distances between P's
    for i, n1 in enumerate(nucls):
        for n2 in nucls[i+1: ]:
            # only for base pairs who have phospate groups
            if n1.idx + n2.idx == 13 and n1.chain_id != n2.chain_id and "P" in n1.named_atoms and "P" in n2.named_atoms:
                diff = n1.named_atoms["P"].pos - n2.named_atoms["P"].pos
                print n1.name, n2.name, sqrt(np.dot(diff, diff)) / divide_by


print_P_distance_dd(nucleotides)
print ""
print_P_distance_dd(ox_nucleotides)
