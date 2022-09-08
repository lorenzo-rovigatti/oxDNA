#!/usr/bin/env python

import sys
import os
import numpy as np
import copy
import string
import argparse

import oxpy

from oxDNA_analysis_tools.UTILS.pdb import Atom, Nucleotide, AminoAcid, FROM_OXDNA_TO_ANGSTROM
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, strand_describe, get_input_parameter
import oxDNA_analysis_tools.UTILS.protein_to_pdb as pro
import oxDNA_analysis_tools.UTILS.utils as utils

DD12_PDB_PATH = "./UTILS/dd12_na.pdb"

def main():
    print("This script isn't working yet.  Check back later.")
    sys.exit()

    parser = argparse.ArgumentParser(prog = os.path.basename(__file__), description="Convert oxDNA files to PDB.  This converter can handle oxDNANM protein simulation files.")
    parser.add_argument('input', type=str, nargs=1, help='the input file used to run the oxDNA simulation')
    parser.add_argument('configuration', type=str, nargs=1, help='the configuration file you wish to convert')
    parser.add_argument('direction', type=str, nargs=1, help='the direction of strands in the oxDNA files, either 35 or 53.  Most oxDNA files are 3-5.')
    parser.add_argument('pdbfiles', type=str, nargs='?', help='PDB files for the proteins present in your structure.  If you have multiple proteins, you must specify multiple PDB files.')
    parser.add_argument('-H', '--hydrogen', action='store_true', default=True, help='if you want to include hydrogen atoms in the output PDB file')
    parser.add_argument('-u', '--uniform-residue-names', action='store_true', default=False, help='if you want to use uniform residue names in the output PDB file')
    parser.add_argument('-o', '--one_file_per_strand', action='store_true', default=False, help='if you want to have one PDB file per strand')
    parser.add_argument('-s', '--same-pdb-all-protein-strands', action='store_true', default=False, help='if you want to have the same PDB file for all protein strands')
    args = parser.parse_args()

    inputfile = args.input[0]
    conf_file = args.configuration[0]
    top_file = get_input_parameter(inputfile, "topology")
    direction = args.direction[0]
    if direction not in ["35", "53"]:
        print("Error: direction must be either 35 or 53")
        sys.exit(1)
    
    if args.pdbfiles:
        protein_pdb_files = args.pdbfiles
    else:
        protein_pdb_files = None

    oxDNA_direction = 1 if direction == "35" else 0
    hydrogen = args.hydrogen
    uniform_residue_names = args.uniform_residue_names
    one_file_per_strand = args.one_file_per_strand
    same_pdb_all_protein_strands = args.same_pdb_all_protein_strands

    top_info, traj_info = describe(top_file, conf_file)
    system, monomers = strand_describe(top_file)

    # Open PDB File of nice lookin duplex
    with open(os.path.join(os.path.dirname(__file__), DD12_PDB_PATH)) as f:
        nucleotides = []
        aminoacids = []
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
        else: 
            bases[n.base] = n
    
    for n in nucleotides:
        n.a1, n.a2, n.a3 = utils.get_orthonormalized_base(n.a1, n.a2, n.a3)


    with oxpy.Context():
        inp = oxpy.InputFile()
        inp.init_from_filename(inputfile)
        inp["trajectory_file"] = traj_info.path
        inp["confs_to_analyse"] = str(1)

        backend = oxpy.analysis.AnalysisBackend(inp)
        backend.read_next_configuration()

        flat_conf = backend.flattened_conf
        conf = backend.config_info()

    com = np.mean(flat_conf.positions, axis = 0)

    # This is awful, but I can't get a box from oxpy at the moment.
    conf_with_box = get_confs(top_info, traj_info, 0, 1)
    box_angstrom = conf_with_box.box * FROM_OXDNA_TO_ANGSTROM

    correct_for_large_boxes = False
    if np.any(box_angstrom[box_angstrom > 999]):
        print("At least one of the box sizes is larger than 999: all the atoms which are outside of the box will be brought back through periodic boundary conditions", file=sys.stderr)
        correct_for_large_boxes = True
    
    if one_file_per_strand:
        out_name = conf_file + "_1.pdb"
    else:
        out_name = conf_file + ".pdb"

    with open(out_name, "w") as f:
        current_base_identifier = 'A'   
        reading_position = 0 
        if protein_pdb_files:
            s_pdbfile = iter(protein_pdb_files)
            pdbfile = next(s_pdbfile)
            
        for sid, strand in enumerate(conf.molecules):
            strand_pdb = []
            monomers = strand.particles
            if not oxDNA_direction:
                monomers = reversed(monomers)

            # protein part
            print(f"INFO: Converting Strand {sid}", end='\r', file=sys.stderr) # would be great if I could get the strand ID as it is in the topology file

            # And that's where we're going to have to stop because the branch of the code with oxpy doesn't support proteins yet.
            # I could do it by making some awful franken-reader that combines the system object which does have good strand IDs with the current oxpy, but it would be gnarly conditionals.

if __name__ == '__main__':
    main()