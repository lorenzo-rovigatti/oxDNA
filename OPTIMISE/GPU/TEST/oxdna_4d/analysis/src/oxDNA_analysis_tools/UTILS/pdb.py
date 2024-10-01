import numpy as np
import itertools
import copy
from typing import List, Dict

BASE_SHIFT = 3.4
COM_SHIFT = 0.5
FROM_OXDNA_TO_ANGSTROM = 8.518
FROM_ANGSTROM_TO_OXDNA = 1. / FROM_OXDNA_TO_ANGSTROM

NAME_TO_BASE = {
        "ADE" : "A",
        "CYT" : "C",
        "GUA" : "G",
        "THY" : "T",
        "URA" : "U",
    }
NAME_TO_AA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
aa_to_number = {'A':-1, 'R':-2, 'N':-3, 'D':-4, 'C':-5, 'E':-6, 'Q':-7, 'G':-8, 'H':-9, 'I':-10, 'L':-11, 'K':-12, 'M':-13, 'F':-14, 'P':-15, 'S':-16, 'T':-17, 'W':-18, 'Y':-19, 'V':-20, 'Z':-21, 'X':0}
BASES = ["A", "T", "G", "C", "U"]
AAS = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

class Atom(object):
    serial_atom = 1

    def __init__(self, pdb_line:str):
        object.__init__(self)
        # http://cupnet.net/pdb-format/
        self.name = pdb_line[12:16].strip()
        # some PDB files have * in place of '
        if "*" in self.name:
            self.name = self.name.replace("*", "'")
        
        self.alternate = pdb_line[16]
        self.residue = pdb_line[17:20].strip()
        self.chain_id = pdb_line[21:22].strip()
        self.residue_idx = int(pdb_line[22:26])
        self.pos = np.array([float(pdb_line[30:38]), float(pdb_line[38:46]), float(pdb_line[46:54])])
        
    def is_hydrogen(self) -> bool:
        return "H" in self.name

    def shift(self, diff:np.ndarray):
        self.pos += diff

    def to_pdb(self, residue_suffix:str, bfactor:float) -> Dict:
        return {
            "name" : self.name,
            "residue_name" : self.residue + residue_suffix,
            "pos" : self.pos,
            "bfactor" : bfactor
        }

    def to_mgl(self):
        colors = {"C" : "0,1,1", "P" : "1,1,0", "O" : "1,0,0", "H" : "0.5,0.5,0.5", "N" : "0,0,1"}
        for c in colors:
            if c in self.name: color = colors[c]
        r = 0.5
        return "%s %s %s @ %f C[%s]" % (self.pos[0], self.pos[1], self.pos[2], r, color)

class PDB_Nucleotide(object):
    
    def __init__(self, name, idx):
        self.name = name.strip()
        if self.name in NAME_TO_BASE.keys():
            self.base = NAME_TO_BASE[self.name]
        elif self.name in BASES:
            self.base = self.name
        else:
            self.base = name.strip('D53')
        self.idx = idx
        self.base_atoms = []
        self.phosphate_atoms = []
        self.sugar_atoms = []
        self.named_atoms = {}
        self.ring_names = ["C2", "C4", "C5", "C6", "N1", "N3"]
        self.chain_id = None

    def get_atoms(self):
        return self.base_atoms + self.phosphate_atoms + self.sugar_atoms

    def add_atom(self, a:Atom):
        if 'P' in a.name or a.name == "HO5'": 
            self.phosphate_atoms.append(a)
        elif "'" in a.name: 
            self.sugar_atoms.append(a)
        else: 
            self.base_atoms.append(a)
        
        self.named_atoms[a.name] = a
        if self.chain_id == None: 
            self.chain_id = a.chain_id

    def get_com(self, atoms:List[Atom]=None):
        if atoms == None: 
            atoms = self.atoms
        com = np.array([0., 0., 0.])
        for a in atoms:
            com += a.pos

        return com / len(atoms)

    def compute_a3(self):
        base_com = self.get_com(self.base_atoms)
        # the O4' oxygen is always (at least for non pathological configurations, as far as I know) oriented 3' -> 5' with respect to the base's centre of mass
        parallel_to = self.named_atoms["O4'"].pos - base_com
        self.a3 = np.array([0., 0., 0.])
        
        for perm in itertools.permutations(self.ring_names, 3):
            p = self.named_atoms[perm[0]]
            q = self.named_atoms[perm[1]]
            r = self.named_atoms[perm[2]]
            v1 = p.pos - q.pos
            v2 = p.pos - r.pos
            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)
            
            if abs(np.dot(v1, v2)) > 0.01 or abs(np.dot(v1, v2)) < 1:
                a3 = np.cross(v1, v2)
                a3 /= np.linalg.norm(a3)
                if np.dot(a3, parallel_to) < 0.: 
                    a3 = -a3
                self.a3 += a3

        self.a3 /= np.linalg.norm(self.a3)

    def compute_a1(self):
        if "C" in self.name or "T" in self.name or "U" in self.name:
            pairs = [ ["N3", "C6"], ["C2", "N1"], ["C4", "C5"] ]
        else:
            pairs = [ ["N1", "C4"], ["C2", "N3"], ["C6", "C5"] ]

        self.a1 = np.array([0., 0., 0.])
        for pair in pairs:
            p = self.named_atoms[pair[0]]
            q = self.named_atoms[pair[1]]
            diff = p.pos - q.pos
            self.a1 += diff

        self.a1 /= np.linalg.norm(self.a1)

    def compute_as(self):
        self.compute_a1()
        self.compute_a3()
        self.a2 = np.cross(self.a3, self.a1)
        self.check = abs(np.dot(self.a1, self.a3))
        
    def correct_for_large_boxes(self, box):
        map(lambda x: x.shift(-np.rint(x.pos / box ) * box), self.atoms)

    def to_pdb(self, print_H:bool, residue_type:str, bfactor:float) -> List[Dict]:
        res:List[Dict] = []
        for a in self.atoms:
            if not print_H and 'H' in a.name:
                continue
            if residue_type == "5": 
                if 'P' in a.name:
                    if a.name == 'P':
                        phosphorus = a
                    continue
                elif a.name == "O5'":
                    O5prime = a
            elif residue_type == "3":
                if a.name == "O3'":
                    O3prime = a
            res.append(a.to_pdb(residue_type, bfactor))
            
        # if the residue is a 3' or 5' end, it requires one more hydrogen linked to the O3' or O5', respectively
        if residue_type == "5":
            new_hydrogen = copy.deepcopy(phosphorus)
            new_hydrogen.name = "HO5'"
            
            # we put the new hydrogen at a distance 1 Angstrom from the O5' oxygen along the direction that, in a regular nucleotide, connects O5' and P
            dist_P_O = phosphorus.pos - O5prime.pos
            dist_P_O /= np.linalg.norm(dist_P_O)
            new_hydrogen.pos = O5prime.pos + dist_P_O
            res.append(new_hydrogen.to_pdb(residue_type, bfactor))
        elif residue_type == "3":
            new_hydrogen = copy.deepcopy(O3prime)
            new_hydrogen.name = "HO3'"
            
            # we put the new hydrogen at a distance 1 Angstrom from the O3' oxygen along a direction which is a linear combination of the three 
            # orientations that approximately reproduce the crystallographic position
            new_distance = 0.2 * self.a2 - 0.2 * self.a1 - self.a3
            new_distance /= np.linalg.norm(new_distance)
            new_hydrogen.pos = O3prime.pos + new_distance
            res.append(new_hydrogen.to_pdb(residue_type, bfactor))

        return res

    def to_mgl(self):
        res = []
        for a in self.atoms:
            res.append(a.to_mgl())

        return "\n".join(res)

    def rotate(self, R):
        com = self.get_com()
        for a in self.atoms:
            a.pos = np.dot(R, a.pos - com) + com

        self.compute_as()

    def set_com(self, new_com):
        com = self.get_com()
        for a in self.atoms:
            a.pos += new_com - com - COM_SHIFT * self.a1

    def set_base(self, new_base_com):
        atoms = [v for k, v in self.named_atoms.items() if k in self.ring_names]
        ring_com = self.get_com(atoms)
        for a in self.atoms:
            a.pos += new_base_com - ring_com - BASE_SHIFT * self.a1

        self.compute_as()

    atoms = property(get_atoms)

class PDB_AminoAcid(object):
    
    def __init__(self, name, idx):
        object.__init__(self)
        self.name = name.strip()
        if self.name in NAME_TO_AA.keys():
            self.base = NAME_TO_AA[self.name]
        elif self.name in AAS:
            self.base = self.name
        else:
            self.base = name[1:]
        self.idx = idx
        self.backbone_atoms = []
        self.sidechain_atoms = []
        self.ca_atom = None
        self.chain_id = None

    def get_atoms(self):
        return self.backbone_atoms + self.sidechain_atoms + [self.ca_atom]

    def add_atom(self, a):
        if a.name in ['N', 'C', 'O']:
            self.backbone_atoms.append(a)
        elif a.name == 'CA': 
            self.ca_atom = a
        else: 
            self.sidechain_atoms.append(a)
        
        if self.chain_id == None: 
            self.chain_id = a.chain_id

    def get_ca_pos(self):
        if self.ca_atom:
            com = np.array(self.ca_atom.pos)
        else:
            raise RuntimeError("No CA atom found in amino acid residue name:{}, id:{}, chain:{}".format(self.name, self.idx, self.chain_id))
        
        return com

    def get_com(self, atoms=None):
        if atoms == None:
            atoms = self.get_atoms()
        com = np.array([0., 0., 0.])
        for a in atoms:
            com += a.pos

        return com / len(atoms)
        
    def correct_for_large_boxes(self, box):
        map(lambda x: x.shift(-np.rint(x.pos / box ) * box), self.atoms)

    def to_pdb(self, print_H:bool, bfactor:float) -> List[Dict]:
        res = []
        for a in self.get_atoms():
            if not print_H and 'H' in a.name:
                continue
            res.append(a.to_pdb('', bfactor))

        return res

    def translate(self, t):
        for a in self.get_atoms():
            a.pos += t

    def rotate(self, R):
        com = self.get_com()
        for a in self.get_atoms():
            a.pos = np.dot(R, a.pos - com) + com

    def set_com(self, new_com):
        com = self.get_com()
        for a in self.get_atoms():
            a.pos += new_com - com

    def set_ca_pos(self, new_ca_pos):
        ca_pos = self.get_ca_pos()
        for a in self.get_atoms():
            a.pos += new_ca_pos - ca_pos