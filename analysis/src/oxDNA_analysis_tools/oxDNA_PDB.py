#!/usr/bin/env python

import sys
import os
import numpy as np
import copy
import string
from math import sqrt, sin
# import Bio.PDB as bio

from oxDNA_analysis_tools.UTILS.pdb import Atom, Nucleotide, AminoAcid, FROM_OXDNA_TO_ANGSTROM

import oxDNA_analysis_tools.UTILS.protein_to_pdb as pro
import oxDNA_analysis_tools.UTILS.base as base
import oxDNA_analysis_tools.UTILS.utils as utils

from oxDNA_analysis_tools.UTILS.readers import LorenzoReader2

DD12_PDB_PATH = "./UTILS/dd12_na.pdb"
protein_pdb_files = []
def print_usage():
        print("USAGE:", file=sys.stderr)
        print("\t%s topology configuration direction pdbfiles(optional, needs flag -p) " % sys.argv[0], file =sys.stderr)
        print("\t[-H\--hydrogens=True] [-u\--uniform-residue-names] [-p\--contains-protein]", file=sys.stderr)
        exit(1)

def parse_options():
    shortArgs = 'H:uosp'
    longArgs = ['hydrogens=','uniform-residue-names', 'one-file-per-strand', 'contains-protein']
    
    opts = {
        "configuration" : "",
        "topology" : "",
        "oxDNA_direction" : True,
        "print_hydrogens" : True,
        "uniform_residue_names" : False,
        "one_file_per_strand" : False,
        "contains_protein" : False,
    }
    
    try:
        import getopt
        args, positional_args = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-H' or k[0] == '--hydrogens':
                k_arg = k[1].lower()
                if k_arg.lower() == "true":
                    opts["print_hydrogens"] = True
                elif k_arg == "false":
                    print >> sys.stderr, "## Hydrogen atoms will *not* be printed"
                    opts["print_hydrogens"] = False
                else:
                    print >> sys.stderr, "The argument of '%s' should be either 'true' or 'false' (got '%s' instead)" % (k[0], k[1])
                    exit(1)
            elif k[0] == '-u' or k[0] == '--uniform-residue-names':
                    opts["uniform_residue_names"] = True
            elif k[0] == '-o' or k[0] == '--one-file-per-strand':
                    opts["one_file_per_strand"] = True
            elif k[0] == '-s' or k[0] == '--same-pdb-all-protein-strands':
                    opts["same_pdb_all_protein_strands"] = True
            elif k[0] == '-p' or k[0] == '--contains-protein':
                    opts["contains_protein"] = True
            
        opts['topology'] = positional_args[0]
        opts['configuration'] = positional_args[1]
        direction = positional_args[2]
        if opts["contains_protein"]:
            protein_pdb_files = positional_args[3:]
        else:
            protein_pdb_files=[]
        
        if direction == "35":
            opts["oxDNA_direction"] = True
        elif direction == "53":
            opts["oxDNA_direction"] = False
        else:
            print >> sys.stderr, "The 'direction' argument should be either 35 or 53"
            exit(1)
    except Exception:
        print_usage()
        
    return opts, protein_pdb_files

def align(full_base, ox_base):
        theta = utils.get_angle(full_base.a3, ox_base._a3)
        # if the two bases are already essentially aligned then we do nothing
        if sin(theta) > 1e-3:
            axis = np.cross(full_base.a3, ox_base._a3)
            axis /= sqrt(np.dot(axis, axis))
            R = utils.get_rotation_matrix(axis, theta)
            full_base.rotate(R)
    
        theta = utils.get_angle(full_base.a1, ox_base._a1)
        if sin(theta) > 1e-3:
            axis = np.cross(full_base.a1, ox_base._a1)
            axis /= sqrt(np.dot(axis, axis))
            R = utils.get_rotation_matrix(axis, theta)
            full_base.rotate(R)

def main():
    opts, protein_pdb_files = parse_options()
        
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
    
    #Lorenzo Reader2 Perhaps?
    try:
        lr = LorenzoReader2(opts['configuration'], opts['topology'])
        s = lr._get_system()
    except Exception as e:
        print("Parser error: %s" % e, file=sys.stderr)
        exit(1)
    
    ox_nucleotides = []
    s.map_nucleotides_to_strands()
    com = np.array([0., 0., 0.])
    for my_strand in s._strands:
        com += my_strand.cm_pos
    com /= s.N_strands
    
    box_angstrom = s._box * FROM_OXDNA_TO_ANGSTROM
    correct_for_large_boxes = False
    if np.any(box_angstrom[box_angstrom > 999]):
        print("At least one of the box sizes is larger than 999: all the atoms which are outside of the box will be brought back through periodic boundary conditions", file=sys.stderr)
        correct_for_large_boxes = True
    
    if opts['one_file_per_strand']:
        out_name = opts['configuration'] + "_1.pdb"
    else:
        out_name = opts['configuration'] + ".pdb"
        
    out = open(out_name, "w+")
    
    current_base_identifier = 'A'   
    reading_position = 0 
    s_pdbfile = iter(protein_pdb_files)
    pdbfile = next(s_pdbfile)
    for s_id, strand in enumerate(s._strands):
        strand_pdb = []
        nucleotides_in_strand = strand._nucleotides
        if not opts['oxDNA_direction']:
            nucleotides_in_strand = reversed(nucleotides_in_strand)
        #Protein
        sys.stdout.write("\rINFO: Converting Strand %i" % strand.index)
        sys.stdout.flush()
        if strand.index < 0:
            #s_pdbfile = protein_pdb_files[abs(strand.index)-1]
            coord = [n.cm_pos for n in strand._nucleotides]  # amino acids only go from nterm to cterm (pdb format does as well)
            #print('oxpositions ', coord)
            #print('reading_pos', reading_position)
            next_reading_position = pro.oxdna_to_pdb(out, coord, pdbfile, com, reading_position)
            if next_reading_position == -1:
                try:
                    pdbfile = next(s_pdbfile)
                    reading_position = 0
                except StopIteration:
                    continue
            #if pdbfile == "end":
                #break
            #else:
                #reading_position = 0
            else:
                 reading_position = next_reading_position

        # DNA
        elif strand.index >= 0:
            for n_idx, nucleotide in enumerate(nucleotides_in_strand, 1):
                nb = base.number_to_base[nucleotide._base]
                my_base = copy.deepcopy(bases[nb])
                my_base.chain_id = s._nucleotide_to_strand[nucleotide.index]
                residue_type = ""
                
                # 3' end
                if nucleotide == strand._nucleotides[0] and not strand._circular:
                    residue_type = "3"
                # 5' end
                elif nucleotide == strand._nucleotides[-1]:
                    residue_type = "5" 
                    
                if opts["uniform_residue_names"] == True:
                    residue_suffix = ""
                else:
                    residue_suffix = residue_type

                align(my_base, nucleotide)
                my_base.set_base((nucleotide.pos_base - com) * FROM_OXDNA_TO_ANGSTROM)

                if correct_for_large_boxes:
                    my_base.correct_for_large_boxes(box_angstrom)
                
                residue_serial = n_idx % 9999
                base_identifier = current_base_identifier
                # Make nucleotide line from pdb.py
                nucleotide_pdb = my_base.to_pdb(base_identifier, opts['print_hydrogens'], residue_serial, residue_suffix, residue_type)
                # Append to strand_pdb
                strand_pdb.append(nucleotide_pdb)
            
            print("\n".join(x for x in strand_pdb), file=out)
            print("TER", file=out)
        
        if opts['one_file_per_strand']:
            out.close()
            print >> sys.stderr, "## Wrote strand %d's data to '%s'" % (s_id + 1, out_name)
            # open a new file if needed
            if strand != s._strands[-1]:
                out_name = opts['configuration'] + "_%d.pdb" % (s_id + 2, )
                out = open(out_name, "w")
        else:
            # we update the base identifier only if a single file is printed
            if current_base_identifier == 'Z':
                current_base_identifier = 'A'
            else:
                current_base_identifier = chr(ord(current_base_identifier) + 1)
    
    # #Must Now renumber and restrand, (sorry)
    out.seek(0)
    outw = open('TM2.pdb', 'w')
    resid, atmid, chainid = -1, 1, -1
    #check against next atom entry
    pres, pchain = 0, 0
    alpha = list(string.ascii_uppercase)
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
                x = ''.join([' ' for i in range(spcs)]) + str(x)
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
            
            
            # nS.init_atom(data[2], [xd, yd, zd], bfactor, occ, alpha[chainid], data[2], serial_number=atmid)
    out.close()
    # io = bio.PDBIO()
    # io.set_structure(nS.get_structure())
    # io.save(out_name)     

    if not opts['one_file_per_strand']:
        outw.close()
        print("\n## Wrote data to '%s'" % out_name, file=sys.stderr)
        
    print("## DONE", file=sys.stderr)

if __name__ == '__main__':
    main()

