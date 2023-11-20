#!/usr/bin/env python

import sys
import os
import numpy as np
import copy
import argparse
from collections import defaultdict
from typing import List, Dict, Tuple
from io import TextIOWrapper

from oxDNA_analysis_tools.UTILS.pdb import Atom, PDB_Nucleotide, PDB_AminoAcid, FROM_OXDNA_TO_ANGSTROM
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, strand_describe, inbox
from oxDNA_analysis_tools.UTILS.data_structures import Strand, Configuration
import oxDNA_analysis_tools.UTILS.utils as utils

DD12_PDB_PATH = "./UTILS/dd12_na.pdb"
RNA_PDB_PATH = "./UTILS/2jxq.pdb"

number_to_DNAbase = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T'}
number_to_RNAbase = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'U'}


base_to_number = {'A' : 0, 'a' : 0, 'G' : 1, 'g' : 1,
                  'C' : 2, 'c' : 2, 'T' : 3, 't' : 3,
                  'U' : 3, 'u' : 3, 'D' : 4}
                  
aa_to_number = {'A':-1, 'R':-2, 'N':-3, 'D':-4, 'C':-5, 
                'E':-6, 'Q':-7, 'G':-8, 'H':-9, 'I':-10, 
                'L':-11, 'K':-12, 'M':-13, 'F':-14, 
                'P':-15, 'S':-16, 'T':-17, 'W':-18, 
                'Y':-19, 'V':-20, 'Z':-21, 'X':0}

number_to_aa = {-1:'A', -2:'R', -3:'N', -4:'D', -5:'C', 
                -6:'E', -7:'Q', -8:'G', -9:'H', -10:'I', 
                -11:'L', -12:'K', -13:'M', -14:'F', 
                -15:'P', -16:'S', -17:'T', -18:'W', 
                -19:'Y', -20:'V', -21:'Z', 0:'X'}

na_pdb_names = ['DA', 'DT', 'DG', 'DC', 'DI', 
                'A', 'U', 'G', 'C', 'I',
                'DA5', 'DT5', 'DG5', 'DC5', 'DI5',
                'DA3', 'DT3', 'DG3', 'DC3', 'DI3',
                'A5', 'U5', 'G5', 'C5', 'I5',
                'A3', 'U3', 'G3', 'C3', 'I3'
                ]

def align(full_base, ox_base):
        theta = utils.get_angle(full_base.a3, ox_base['a3'])
        axis = np.cross(full_base.a3, ox_base['a3'])
        axis /= np.linalg.norm(axis)
        R = utils.get_rotation_matrix(axis, theta)
        full_base.rotate(R)
    
        theta = utils.get_angle(full_base.a1, ox_base['a1'])
        axis = np.cross(full_base.a1, ox_base['a1'])
        axis /= np.linalg.norm(axis)
        R = utils.get_rotation_matrix(axis, theta)
        full_base.rotate(R)

def get_nucs_from_PDB(file:str) -> List[PDB_Nucleotide]:
    """
        Extract nucleotides from a PDB file

        Parameters:
            file (str) : The path to the PDB file to read

        Returns:
            List[PDB_Nucleotide] : A list of nucleotide objects from the PDB file
    """
    with open(file) as f:
        nucleotides = []
        old_residue = ""
        for line in f.readlines():
            if len(line) > 77 and line[17:20].strip() in na_pdb_names:
                na = Atom(line)
                if na.residue_idx != old_residue:
                    nn = PDB_Nucleotide(na.residue, na.residue_idx)
                    nucleotides.append(nn)
                    old_residue = na.residue_idx
                nn.add_atom(na)

    return nucleotides

def choose_reference_nucleotides(nucleotides:List[PDB_Nucleotide]) -> Dict[str, PDB_Nucleotide]:
    """
        Find nucleotides that most look like an oxDNA nucleotide (orthogonal a1 and a3 vectors).

        Parameters:
            nucleotides (List[PDB_Nucleotide]) : List of nucleotides to compare.
        
        Returns:
            Dict[str, PDB_Nucleotide] : The best nucleotide for each type in the format `{'C' : PDB_Nucleotide}`.
    """
    bases = {}
    for n in nucleotides:
        n.compute_as()
        if n.base in bases:
            if n.check < bases[n.base].check: # Find the most orthogonal a1/a3 in the reference
                bases[n.base] = copy.deepcopy(n)
                bases[n.base].a1, bases[n.base].a2, bases[n.base].a3 = utils.get_orthonormalized_base(n.a1, n.a2, n.a3)
        else:
            bases[n.base] = copy.deepcopy(n)
            bases[n.base].a1, bases[n.base].a2, bases[n.base].a3 = utils.get_orthonormalized_base(n.a1, n.a2, n.a3)

    return bases

def get_AAs_from_PDB(pdbfile:str, start_res:int=0, n_res:int=-1) -> Tuple[int, List[PDB_AminoAcid]]:
    """
        Get amino acid descriptions from a PDB file.

        Parameters:
            pdbfile (str) : File path to PDB file
            start_res (int) : Line number to start parsing from, default 0.
            n_res (int) : Get only this many residues, default until end of file.

        Returns
            (Tuple[int, List[PDB_AminoAcid]]) : The line number in the PDB file the read left off at and a list of amino acids extracted from the file.
    """
    with open(pdbfile) as pdbf:
        amino_acids = []
        old_residue = ''
        lines = pdbf.readlines()
        for l in lines:
            if l.startswith('ATOM'):
                # We assume that this file just contains proteins.  Die if you find a nucleic acid
                if l[17:20].strip() in na_pdb_names:
                    raise RuntimeError("Invalid residue name {} in {}. The reference PDB file must only contain protein residues.".format(l[17:20].strip(), pdbfile))
                
                a = Atom(l)
                if a.residue_idx != old_residue:
                    aa = PDB_AminoAcid(a.residue, a.residue_idx)
                    amino_acids.append(aa)
                    old_residue = a.residue_idx
                aa.add_atom(a)
    
    if n_res == -1:
        end = len(amino_acids)
    else:
        end = start_res + n_res

    if end == len(amino_acids):
        next_pos = -1
    else:
        next_pos = end

    return next_pos, amino_acids[start_res:end]

def peptide_to_pdb(strand:Strand, conf:Configuration, pdbfile:str, reading_position:int) -> Tuple[int, List[PDB_AminoAcid]]:
    """
        Convert a Strand object to an all-atom representation based on a reference pdb file.
        
        Parameters:
            strand (Strand) : The strand to convert
            conf (Configuraiton) : The entire oxDNA configuration
            pdbfile (str) : Path to a pdb file to get r-group orientations from
            reading_position (int) : Starting residue in the PDB file for this strand

        Returns:
            (Tuple[int, List[PDB_AminoAcid]]) : The next reading position in the PDB file (or -1 if the file was finished) and the list of PDB-ready amino acid objects
    """
    coord = np.array([conf.positions[m.id] for m in strand.monomers])  # amino acids only go from nterm to cterm (pdb format does as well)
    coord = coord * FROM_OXDNA_TO_ANGSTROM

    reading_position, amino_acids = get_AAs_from_PDB(pdbfile, reading_position, len(coord))
    
    # Translate CA COM of PDB structure to COM of oxDNA structure
    ca_poses = np.array([a.get_ca_pos() for a in amino_acids])
    pdb_com = np.mean(ca_poses, axis=0)
    ox_com = np.mean(coord, axis=0)
    centered_ca_poses = np.array([a.get_ca_pos() for a  in amino_acids]) - pdb_com
    centered_ox_poses = coord - ox_com
    centered_atom_poses = np.array([a.pos - pdb_com for aa in amino_acids for a in aa.get_atoms()])

    # Get rotation matrix
    M = np.dot(centered_ca_poses.T, centered_ox_poses)
    u,s,v = np.linalg.svd(M, compute_uv=True)
    R = np.dot(v.T, u.T)

    # reposition and rotate each amino acid
    for i, aa in enumerate(amino_acids):
        aa.set_ca_pos(coord[i])
        aa.rotate(R)

    return(reading_position, amino_acids)

def write_strand_to_PDB(strand_pdb:List[Dict], chain_id:str, atom_counter:int, out:TextIOWrapper) -> int:
    """
        Write a list of nucleotide property dictionaries as a new chain to an open PDB file

        Parameters:
            strand_pdb (List[Dict]) : A list of dicts with nucleotide properties which define a single strand.
            chain_id (str) : The current chainID
            atom_counter (int) : The starting atom ID for this chain
            out (io.TextIOWrapper) : An open file handle to write to

        Returns:
            (int) : The next atom ID
    """
    #re-index and create PDB string
    for nid, n in enumerate(strand_pdb, 1):
        for a in n:
            print("{:6s}{:5d} {:^4s}{:1s}{:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
                .format(
                    "ATOM",            #record
                    atom_counter,      #atom_id
                    a['name'],         #atom_name
                    " ",               #alt_loc
                    a['residue_name'], #res_name
                    chain_id,          #chain_id
                    nid,               #res_id
                    " ",               #ins_code
                    a['pos'][0],       #coord_x
                    a['pos'][1],       #coord_y
                    a['pos'][2],       #coord_z
                    1.00,              #residency
                    a['bfactor'],      #b-factor
                    " ", " "           #element,charge
                ),
                file=out
            )
            atom_counter = (atom_counter+1) % 9999
    print("TER", file=out)

    return(atom_counter)

def cli_parser(prog="oxDNA_PDB.py"):
    parser = argparse.ArgumentParser(prog=prog, description="Convert oxDNA files to PDB.  This converter can handle oxDNANM protein simulation files.")
    parser.add_argument('topology', type=str,
                        help='the oxDNA topology file for the structure')
    parser.add_argument('configuration', type=str,
                        help='the configuration file you wish to convert')
    parser.add_argument('direction', type=str,
                        help='the direction of strands in the oxDNA files, either 35 or 53.  Most oxDNA files are 3-5.')
    parser.add_argument('pdbfiles', type=str, nargs='*',
                        help='PDB files for the proteins present in your structure.  The strands in the PDB file(s) must be in the same order as your oxDNA file. If there are multiple of the same protein, you must provide that PDB file that many times.')
    parser.add_argument('-o', '--output', type=str, 
                        help='The name of the output pdb file.  Defaults to name of the configuration+.pdb')
    parser.add_argument('-d', '--output_direction', type=str,
                        help='Direction to save nucleic acid strands in.  Should be "53" or "35".  Defaults to same as input direction.')
    parser.add_argument('-H', '--hydrogen', action='store_true', default=True,
                        help='if you want to include hydrogen atoms in the output PDB file')
    parser.add_argument('-u', '--uniform-residue-names', action='store_true', default=False,
                        help='if you want to use uniform residue names in the output PDB file')
    parser.add_argument('-1', '--one_file_per_strand', action='store_true',
                        default=False, help='if you want to have one PDB file per strand')
    parser.add_argument('-r', '--rmsf-file', dest='rmsf_bfactor', type=str, nargs=1, 
                        help='A RMSF file from deviations.  Will be used to fill the b-factors field in the output PDB (only for D(R)NA)')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    # Parse positional arguments
    top_file = args.topology
    conf_file = args.configuration
    direction = args.direction
    if direction not in ["35", "53"]:
        raise RuntimeError("Error: Direction must be either 35 or 53")
    if args.pdbfiles:
        protein_pdb_files = args.pdbfiles
    else:
        protein_pdb_files = None

    # Parse optional arguments
    if args.output:
        out_basename = args.output.rstrip('.pdb')
    else:
        out_basename = conf_file
    reverse = False
    if args.output_direction:
        if args.output_direction not in ["35", "53"]:
            raise RuntimeError("Error: Output direction must be either 35 or 53")
        if args.output_direction != direction:
            reverse = True

    hydrogen = args.hydrogen
    uniform_residue_names = args.uniform_residue_names
    one_file_per_strand = args.one_file_per_strand
    rmsf_file = args.rmsf_bfactor

    # Open PDB File of nice lookin duplexes to get base structures from
    DNAnucleotides = get_nucs_from_PDB(os.path.join(os.path.dirname(__file__), DD12_PDB_PATH))
    RNAnucleotides = get_nucs_from_PDB(os.path.join(os.path.dirname(__file__), RNA_PDB_PATH))

    DNAbases = choose_reference_nucleotides(DNAnucleotides)
    RNAbases = choose_reference_nucleotides(RNAnucleotides)

    # Read oxDNA configuration
    system, _ = strand_describe(top_file)
    ti, di = describe(top_file, conf_file)
    conf = get_confs(ti, di, 0, 1)[0]
    conf = inbox(conf, center=True)
    box_angstrom = conf.box * FROM_OXDNA_TO_ANGSTROM

    # Handle RMSF -> bFactor conversion
    if rmsf_file:
        with open(rmsf_file) as f:
            try:
                # .json format from oat deviations
                substrings = f.read().split("[")[1].split("]")[0].split(",")
            except Exception as e:
                raise RuntimeError("Parsing error in RMSF file. Invalid Format: %s" % e)
            try:
                rmsf_per_nucleotide = {i: float(s) for i, s in enumerate(substrings)}
            except Exception as e:
                raise RuntimeError("Parsing error in RMSF file. Conversion to float failed : %s" % e)
    else:
        rmsf_per_nucleotide = defaultdict(lambda: 1.00)

    # Process optional conditionals
    correct_for_large_boxes = False
    if np.any(box_angstrom[box_angstrom > 999]):
        print("INFO: At least one of the box sizes is larger than 999: all the atoms which are outside of the box will be brought back through periodic boundary conditions", file=sys.stderr)
        correct_for_large_boxes = True
    
    if one_file_per_strand:
        out_name = out_basename+"_{}.pdb".format(system.strands[0].id)
    else:
        out_name = out_basename+".pdb"

    # Start writing the output file
    with open(out_name, 'w+') as out:
        reading_position = 0 
        chain_id = 'A'
        atom_counter = 1

        # Iterate over strands in the oxDNA file
        for strand in system.strands:
            strand_pdb = []
            nucleotides_in_strand = strand.monomers
            sequence = [n.btype for n in nucleotides_in_strand]
            isDNA = True #This should be in the strand parser instead.
            if 'U' in sequence or 'u' in sequence: #Turns out, this is a bad assumption but its all we got.
                isDNA = False

            print("\rINFO: Converting strand {}".format(strand.id), file=sys.stderr)

            # Handle protein
            if strand.id < 0 and protein_pdb_files:
                # Map oxDNA configuration onto R-group orientations from pdb file
                s_pdbfile = iter(protein_pdb_files)
                pdbfile = next(s_pdbfile)
                reading_position, amino_acids = peptide_to_pdb(strand, conf, pdbfile, reading_position)
                if reading_position == -1:
                    try:
                        pdbfile = next(s_pdbfile)
                        reading_position = 0
                    except StopIteration:
                        protein_pdb_files = [] #had better be nucleic acids next or we're going to the error in the elif.

                # Convert AminoAcid objects to write-ready dicts    
                for aa in amino_acids:
                    amino_acid_pdb = aa.to_pdb(
                        hydrogen,
                        bfactor=rmsf_per_nucleotide[aa.idx],
                    )
                    strand_pdb.append(amino_acid_pdb)
                
                # Write residue to file
                atom_counter = write_strand_to_PDB(strand_pdb, chain_id, atom_counter, out)
                
            elif strand.id < 0 and not protein_pdb_files:
                raise RuntimeError("You must provide PDB files containing just the protein for each protein in the scene.")

            # Nucleic Acids
            elif strand.id >= 0:
                for nucleotide in nucleotides_in_strand:
                    # Get paragon DNA or RNA nucleotide
                    if type(nucleotide.btype) != str:
                        if isDNA:
                            nb = number_to_DNAbase[nucleotide.btype]
                        else:
                            nb = number_to_RNAbase[nucleotide.btype]
                    else: 
                        nb = nucleotide.btype
                    
                    if isDNA:
                        my_base = copy.deepcopy(DNAbases[nb])
                    else:
                        my_base = copy.deepcopy(RNAbases[nb])

                    # end residue identifiers
                    residue_type = ""
                    if not uniform_residue_names:
                        if nucleotide == strand.monomers[0] and not strand.is_circular():
                            residue_type = "3"
                        elif nucleotide == strand.monomers[-1]:
                            residue_type = "5"

                    nuc_data = {
                        'pos' : conf.positions[nucleotide.id],
                        'a1' : conf.a1s[nucleotide.id],
                        'a3' : conf.a3s[nucleotide.id] 
                    }

                    # Align paragon nucleotide to the oxDNA nucleotide
                    my_base.set_com(nuc_data['pos'] * FROM_OXDNA_TO_ANGSTROM)
                    align(my_base, nuc_data)

                    if correct_for_large_boxes:
                        my_base.correct_for_large_boxes(box_angstrom)

                    # Turn nucleotide object into a dict for output
                    nucleotide_pdb = my_base.to_pdb(
                        hydrogen,
                        residue_type,
                        bfactor=rmsf_per_nucleotide[nucleotide.id],
                    )
                    strand_pdb.append(nucleotide_pdb)

                # Reverse the strand if the nucleotides should be flipped
                if reverse:
                    strand_pdb = strand_pdb[::-1]

                # Write the current strand to the pdb file.
                atom_counter = write_strand_to_PDB(strand_pdb, chain_id, atom_counter, out)
                
            # Either open a new file or increment chain ID
            # Chain ID can be any alphanumeric character.  Convention is A-Z, a-z, 0-9
            if one_file_per_strand:
                out.close()
                print("INFO: Wrote strand {}'s data to {}".format (strand.id, out_name))
                chain_id = 'A'
                if strand != system.strands[-1]:
                    out_name = out_basename + "_{}.pdb".format(strand.id, )
                    out = open(out_name, "w")
            else:
                chain_id = chr(ord(chain_id)+1)
                if chain_id == chr(ord('Z')+1):
                    chain_id = 'a'
                elif chain_id == chr(ord('z')+1):
                    chain_id = '1'
                elif chain_id == chr(ord('0')+1):
                    print("WARNING: More than 62 chains identified, looping chain identifier...", file=sys.stderr)
                    chain_id = 'A'

    print("INFO: Wrote data to '{}'".format(out_name), file=sys.stderr)
        
    print("\nINFO: DONE", file=sys.stderr)

if __name__ == '__main__':
    main()
