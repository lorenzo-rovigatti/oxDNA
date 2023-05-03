#!/usr/bin/env python

import sys
import os
import numpy as np
import copy
import argparse
import string
from collections import defaultdict
from typing import List, Dict

from oxDNA_analysis_tools.UTILS.pdb import Atom, Nucleotide, AminoAcid, FROM_OXDNA_TO_ANGSTROM
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, strand_describe, inbox
import oxDNA_analysis_tools.UTILS.protein_to_pdb as pro
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

def get_nucs_from_PDB(file:str) -> List[Nucleotide]:
    """
        Extract nucleotides from a PDB file

        Parameters:
            file (str) : The path to the PDB file to read

        Returns:
            List[Nucleotide] : A list of nucleotide objects from the PDB file
    """
    with open(file) as f:
        nucleotides = []
        old_residue = ""
        for line in f.readlines():
            if len(line) > 77 and line[17:20].strip() in na_pdb_names:
                na = Atom(line)
                if na.residue_idx != old_residue:
                    nn = Nucleotide(na.residue, na.residue_idx)
                    nucleotides.append(nn)
                    old_residue = na.residue_idx
                nn.add_atom(na)

    return nucleotides

def choose_reference_nucleotides(nucleotides:List[Nucleotide]) -> Dict[str, Nucleotide]:
    """
        Find nucleotides that most look like an oxDNA nucleotide (orthogonal a1 and a3 vectors).

        Parameters:
            nucleotides (List[Nucleotide]) : List of nucleotides to compare.
        
        Returns:
            Dict[str, Nucleotide] : The best nucleotide for each type in the format `{'C' : Nucleotide}`.
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

def cli_parser(prog="oxDNA_PDB.py"):
    parser = argparse.ArgumentParser(prog=prog, description="Convert oxDNA files to PDB.  This converter can handle oxDNANM protein simulation files.")
    parser.add_argument('topology', type=str,
                        help='the oxDNA topology file for the structure')
    parser.add_argument('configuration', type=str,
                        help='the configuration file you wish to convert')
    parser.add_argument('direction', type=str,
                        help='the direction of strands in the oxDNA files, either 35 or 53.  Most oxDNA files are 3-5.')
    parser.add_argument('pdbfiles', type=str, nargs='?',
                        help='PDB files for the proteins present in your structure.  If you have multiple proteins, you must specify multiple PDB files.')
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
        protein_pdb_files = args.pdbfiles.split(' ')
    else:
        protein_pdb_files = None

    # Parse optional arguments
    if args.output:
        out_basename = args.output.strip('.pdb')
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

    # get protein reference files
    if protein_pdb_files:
        s_pdbfile = iter(protein_pdb_files)
        pdbfile = next(s_pdbfile)

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

        # Iterate over strands in the oxDNA file
        for strand in system.strands:
            strand_pdb = []
            nucleotides_in_strand = strand.monomers
            sequence = [n.type for n in nucleotides_in_strand]
            isDNA = True
            if 'U' in sequence or 'u' in sequence:
                isDNA = False

            print("\rINFO: Converting strand {}".format(strand.id), file=sys.stderr)

            # Handle protein
            if strand.id < 0 and protein_pdb_files:
                coord = np.array([conf.positions[m.id] for m in strand.monomers])  # amino acids only go from nterm to cterm (pdb format does as well)
                next_reading_position = pro.oxdna_to_pdb(out, coord, pdbfile, np.array([0, 0, 0]), reading_position)
                if next_reading_position == -1:
                    try:
                        pdbfile = next(s_pdbfile)
                        reading_position = 0
                    except StopIteration:
                        continue
                else:
                     reading_position = next_reading_position
            elif strand.id < 0 and not protein_pdb_files:
                raise RuntimeError("You must provide PDB files containing just the protein for each protein in the scene.")

            # Nucleic Acids
            elif strand.id >= 0:
                for nucleotide in nucleotides_in_strand:
                    if type(nucleotide.type) != str:
                        if isDNA:
                            nb = number_to_DNAbase[nucleotide.type]
                        else:
                            nb = number_to_RNAbase[nucleotide.type]
                    else: 
                        nb = nucleotide.type
                    
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

                    my_base.set_com(nuc_data['pos'] * FROM_OXDNA_TO_ANGSTROM)
                    align(my_base, nuc_data)

                    if correct_for_large_boxes:
                        my_base.correct_for_large_boxes(box_angstrom)

                    # Make nucleotide line from pdb.py
                    nucleotide_pdb = my_base.to_pdb(
                        hydrogen,
                        residue_type,
                        bfactor=rmsf_per_nucleotide[nucleotide.id],
                    )
                    # Append to strand_pdb
                    strand_pdb.append(nucleotide_pdb)

                if reverse:
                    strand_pdb = strand_pdb[::-1]

                #re-index and create PDB string
                chain_id = 'A'
                atom_counter = 1
                for nid, n in enumerate(strand_pdb,1):
                    for a in n:
                        print("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:1s}{:4d}{:1s}  {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
                            .format("ATOM", 
                                atom_counter, 
                                a['name'], 
                                " ", 
                                a['residue_name'], 
                                chain_id,
                                " ", 
                                nid, 
                                " ", 
                                a['pos'][0], 
                                a['pos'][1], 
                                a['pos'][2], 
                                0.00, 
                                a['bfactor'], 
                                " ", " ", " "
                            ), 
                            file=out
                        )
                        atom_counter = (atom_counter+1) % 9999

                print("TER", file=out)
                
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

        if protein_pdb_files:  
            # #Must Now renumber and restrand, (sorry)
            out.seek(0)
            outw = open('.'.join(top_file.split('.')[:-1]) + "pdb", 'w+')
            resid, atmid, chainid = -1, 1, -1
            #check against next atom entry
            pres, pchain = 0, 0
            alpha = list(string.ascii_uppercase)
            a2 = copy.deepcopy(alpha)
            for i in range(26):
                alpha += [a2[i]+x for x in a2]
            for line in out:
                write_chain_end = False
                if line.startswith('ATOM'):
                    data = line.split()
                    cres_id = data[5]
                    curr_chainid = data[4]
                    if curr_chainid != pchain:
                        if chainid != -1:
                            write_chain_end = True
                        chainid += 1
                        pchain = curr_chainid 
                    if cres_id != pres:
                        resid += 1
                        pres = cres_id
                    data[1] = str(atmid)
                    data[5] = str(resid + 1)
                    data[4] = alpha[chainid]

                    coord = line[29:53]
                    xd, yd, zd = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    for x in [xd, yd, zd]:
                        spcs = 7-len(str(x))
                        x = ''.join([' ' for _ in range(spcs)]) + str(x)
                    #print(xd, yd, zd)
                    if len(list(data[-1])) == 1:
                        del data[-1]

                    if len(data) == 8:
                        try:
                            bfactor = float(data[7])
                            occ = float(data[6])
                        except ValueError:
                            occ = data[7][:4]
                            bfactor = data[7][4:]
                    elif len(data) == 9:
                        try:
                            bfactor = float(data[8])
                            occ = float(data[7])
                        except ValueError:
                            occ = data[8][:4]
                            bfactor = data[8][4:]
                    elif len(data) == 10:
                        try:
                            bfactor = float(data[9])
                            occ = float(data[8])
                        except ValueError:
                            occ = data[9][:4]
                            bfactor = data[9][4:]
                    elif len(data) == 11:
                        try:
                            bfactor=float(data[10])
                            occ = float(data[9])
                        except ValueError:
                            occ = data[10][:4]
                            bfactor = data[10][4:]
                    else: 
                        bfactor=0
                        occ = 1

                    for x in [occ, bfactor]:
                        spcs = 4-len(str(x))
                        x = ''.join([' ' for i in range(spcs)]) + str(x)

                    # data indice -> PDB FIELD LENGTH
                    dls = {1: 6, 2:4, 3:3, 4:2, 5:6}
                    for i in range(1,6):
                        x = len(data[i])
                        if x != dls[i]:
                            diff = dls[i] - x
                            empty = [' ' for j in range(diff)]
                            data[i] = ''.join(empty)+data[i]
                            #print(data[i])


                    print('ATOM', data[1], data[2], data[3], data[4], data[5], coord, occ, bfactor, file=outw)
                    if write_chain_end:
                        print('TER ', data[1], '    ', data[4], data[5], file=outw)
                    atmid += 1

            outw.close()
    print("INFO: Wrote data to '{}'".format(out_name), file=sys.stderr)
        
    print("\nINFO: DONE", file=sys.stderr)

if __name__ == '__main__':
    main()
