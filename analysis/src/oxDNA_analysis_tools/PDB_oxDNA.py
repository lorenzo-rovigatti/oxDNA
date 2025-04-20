import time
start_time = time.time()
import argparse
import numpy as np
from os import path
from typing import Tuple, List, Dict
from dataclasses import dataclass
from itertools import permutations
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.RyeReader import write_conf, write_top
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, System, Strand, Monomer

BASE_SHIFT = 1.13
COM_SHIFT = 0.5
OXDNA_TO_ANGSTROM = 8.518
ANGSTROM_TO_OXDNA = 1. / OXDNA_TO_ANGSTROM

@dataclass
class Atom():
    """A dataclass defining an atom"""
    id : int
    type : str
    alt : str
    resn : str
    chain : str
    resi : int
    pos : np.ndarray

class Residue():
    """A class defining a residue"""
    resn : str
    resi : int
    atoms : List[Atom]
    sugar_atoms : List[Atom]
    base_atoms : List[Atom]
    phosphate_atoms : List[Atom]
    atom_lookup : Dict[str, Atom]

    def __init__(self, resn, resi):
        self.resn = resn.strip('D35*')
        self.resi = resi
        self.atoms = []
        self.sugar_atoms = []
        self.phosphate_atoms = []
        self.base_atoms = []
        self.atom_lookup = {}
        if self.resn not in {'A', 'G', 'C', 'T', 'U'}:
            raise RuntimeError(f'Residue name {resn} at position {resi} could not be mapped to any standard nucleic acid name.')

    def parse_nucleic_atoms(self):
        for a in self.atoms:
            if 'P' in a.type or a.type == "HO5'":
                self.phosphate_atoms.append(a)
            elif "'" in a.type:
                self.sugar_atoms.append(a)
            else:
                self.base_atoms.append(a)

            self.atom_lookup[a.type] = a

    def calc_ox_properties(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        self.parse_nucleic_atoms() # should make this work for protein ANMs at some point
        ring_names = ["C2", "C4", "C5", "C6", "N1", "N3"] # These are found in all rings

        # Take atoms and figure out the pos, a1 and a3
        pos = np.mean([a.pos for a in self.atoms], axis=0) * ANGSTROM_TO_OXDNA # Center of mass of the atoms
        #sugar_com = np.mean([a.pos for a in self.sugar_atoms], axis=0)

        # Figure out the a3 vector by finding the average normal vector to the base
        a3 = np.zeros(3)
        # In the Taco script they use O4' and the center of the base as the reference, however this is a bad assumption for RNA
        # Instead I'll use C5' and C3' because that seems much more likely to point up the helix.
        ref_vec = self.atom_lookup["C5'"].pos - self.atom_lookup["C3'"].pos
        for perm in permutations(ring_names, 3):
            p = self.atom_lookup[perm[0]]
            q = self.atom_lookup[perm[1]]
            r = self.atom_lookup[perm[2]]
            v1 = p.pos - q.pos
            v2 = p.pos - r.pos
            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)
            a3_i = np.cross(v1, v2)
            a3_i /= np.linalg.norm(a3_i)
            if np.dot(a3_i, ref_vec) < 0:
                a3_i *= -1
            a3 += a3_i
        a3 /= np.linalg.norm(a3)

        # Figure out the a1 vector based on the ring atom positions
        if "C" in self.resn or "T" in self.resn or "U" in self.resn:
            pairs = [ ["N3", "C6"], ["C2", "N1"], ["C4", "C5"] ]
        else:
            pairs = [ ["N1", "C4"], ["C2", "N3"], ["C6", "C5"] ]
        a1 = np.zeros(3)
        for pair in pairs:
            p = self.atom_lookup[pair[0]]
            q = self.atom_lookup[pair[1]]
            v1 = p.pos - q.pos
            v1 /= np.linalg.norm(v1)
            a1 += v1
        a1 /= np.linalg.norm(a1)

        return (pos, a1, a3)
    
    def to_monomer(self, strand:Strand) -> Monomer:
        p, a1, a3 = self.calc_ox_properties()
        return Monomer(self.resi, self.resn, strand, None, None, None, p, a1, a3)

def parse_atom(l:str):
    return Atom(
        int(l[6:11].strip()),                                          # atom index
        l[12:16].strip(),                                              # atom name
        l[16].strip(),                                                 # alternate location
        l[17:20].strip(),                                              # resname
        l[21],                                                         # chain
        int(l[22:26].strip()),                                         # resid
        np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])  # xyz coords
    )

def PDB_oxDNA(pdb_str:str, old_top:bool=False) -> Tuple[List[Configuration], List[System]]:
    """
    Convert a PDB string to an oxDNA Configuration/System

    Parameters:
        pdb_str (str) : The contents of a PDB file as a string
        old_top (bool) : (optional) If True, create an old-style topology. default=False

    Returns:
        (tuple(System, Configuration)) : The system and configuration in oxDNA ready for write-out
    """
    
    # Cleanup function to run whenever we end a strand
    def end_strand():
        nonlocal a, r, strand, sys, prev_chain, prev_resi
        monomer = r.to_monomer(strand)
        if "O2'" in r.atom_lookup.keys(): # Check for the 2' hydroxyl 
            strand.type='RNA'
        strand.append(monomer)
        sys.append(strand)
        strand = Strand(strand.id+1)
        prev_chain = a.chain
        prev_resi = -1

    # Cleanup function to run whenever we end a model
    def end_system():
        nonlocal sys, old_top, systems, configs
        if old_top:
            for s in sys.strands:
                s.monomers = s.monomers[::-1]
                s.set_old(True)

        positions = np.array([m.pos for s in sys.strands for m in s])
        box = 1.5*(np.max(positions) - np.min(positions))
            
        conf = Configuration(0, np.array([box, box, box]), np.array([0, 0, 0]), positions, np.array([m.a1 for s in sys.strands for m in s]), np.array([m.a3 for s in sys.strands for m in s]))
        systems.append(sys)
        configs.append(conf)
        sys = System('', [])

    systems = []
    configs = []
    sys = System('')
    strand = Strand(1)
    prev_resi:int = -1
    prev_chain:str = ''

    lines = pdb_str.split('\n')
    for l in lines:
        if l.startswith('ATOM'): # We have found an atom
            a = parse_atom(l)

            # Catch case where there was no TER on the previous strand
            if prev_chain != '':
                if a.chain != prev_chain and len(strand) != 0:
                    end_strand()
            
            # Use the first available alternate position
            if a.alt:
                if not (a.alt == 'A' or a.alt == '1'):
                    continue 
                log(f"Alternate location for atom {a.id} of residue {a.resi} encountered, using location {a.alt}")
            
            # We're in a new residue
            if a.resi != prev_resi:
                if prev_resi != -1:
                    monomer = r.to_monomer(strand)
                    strand.append(monomer)
                r = Residue(a.resn, a.resi)

            r.atoms.append(a)
            prev_resi = a.resi
            prev_chain = a.chain
            continue

        # End of a chain
        elif l.startswith('TER'):
            end_strand()
            continue

        # End the system, start a new one (if there is one)
        # Also catch the last strand if there was no TER identifier
        elif l.startswith('END'):
            if len(strand) > 0:
                end_strand()
            if len(sys.strands) > 0:
                end_system()
            continue
    
    # Catch the case where there was no TER identifier
    if len(strand) > 0: 
        end_strand()
        
    # Catch the case where there was no END identifier
    if len(sys) > 0:
        end_system()
    
    return configs, systems

# This is what gets picked up by the cli documentation builder
def cli_parser(prog="program_name"):
    parser = argparse.ArgumentParser(prog = prog, description="Convert a PDB file to oxDNA")
    parser.add_argument('pdb_file', type=str, help='The pdb file you wish to convert')
    parser.add_argument('-o', '--output', metavar='output_file', help='The filename to save the output oxDNA files to')
    parser.add_argument('-b', '--backward', action='store_true', default=False, help="Use the old topology format? (strands will be listed 3'-5')")
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser

def main():
    # Get arguments from the CLI
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    #run system checks
    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy"])

    # Parse CLI input
    pdb_file = args.pdb_file

    old_top = args.backward
 
    # -o names the output file
    if args.output:
        outbase = args.output
    else:
        outbase = pdb_file.replace('.pdb', '')
        log(f"No outfile name provided, defaulting to \"{outbase}\"")

    # Get the string from the PDB
    with open(pdb_file, 'r') as f:
        pdb_str = f.read()

    # Get the oxDNA-style data structures
    configs, systems = PDB_oxDNA(pdb_str, old_top=old_top)

    for i, (conf, sys) in enumerate(zip(configs, systems)):
        # Figure out the names
        if len(configs) == 1:
            sysn = ''
        else:
            sysn = '_'+str(i)
        outtop = outbase+sysn+'.top'
        outdat = outbase+sysn+'.dat'
        # Write out the files
        write_top(outtop, sys, old_format=old_top)
        write_conf(outdat, conf)
        log(f"Wrote outfiles {outtop}, {outdat}")

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()