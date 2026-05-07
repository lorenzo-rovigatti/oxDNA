#!/usr/bin/env python
import time
start_time = time.time()
import os
import re
import inspect
import itertools
import numpy as np
import copy
import argparse
from collections import defaultdict
from typing import List, Dict, Tuple, Union, Iterator
from io import TextIOWrapper
from pathlib import Path

from oxDNA_analysis_tools.UTILS.pdb import Atom, PDB_Nucleotide, PDB_AminoAcid, FROM_OXDNA_TO_ANGSTROM, NAME_TO_AA
from oxDNA_analysis_tools.UTILS.RyeReader import get_confs, describe, strand_describe, inbox
from oxDNA_analysis_tools.UTILS.data_structures import Strand, Configuration, System
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.mmcif import (
    write_data_block_header, write_atom_site_block,
    write_entity_block, write_entity_poly_block,
    parse_atom_site, is_mmcif,
)
import oxDNA_analysis_tools.UTILS.utils as utils

# This is terrible code, but it *does* get the path of this file whether its run as a script or imported as a module or whatever Sphinx does to build docs.
PDB_PATH = Path(os.path.dirname(os.path.realpath(inspect.getfile(inspect.currentframe())))+"/UTILS/pdb_templates/")

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
                'A3', 'U3', 'G3', 'C3', 'I3',
                'RA5', 'RU5', 'RG5', 'RC5', 'RI5',
                'RA3', 'RU3', 'RG3', 'RC3', 'RI3',
                'PSU', '5HC', 'S6G'
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
            if 'ATOM' in line[0:7] or 'HETATM' in line[0:7] and line[17:20].strip() in na_pdb_names:
                na = Atom(line)
                if na.residue_idx != old_residue:
                    nn = PDB_Nucleotide(na.residue, na.residue_idx)
                    nucleotides.append(nn)
                    old_residue = na.residue_idx
                nn.add_atom(na)

    return nucleotides

# Helper functions for getting specific atoms relative to COM
get_base_center = lambda nuc: np.mean([a.pos - nuc.get_com() for a in nuc.base_atoms], axis=0)
get_base_C2 = lambda nuc: nuc["C2"].pos - nuc.get_com()
get_base_edge = lambda nuc: nuc['N1'].pos - nuc.get_com() if any((c in nuc.name for c in ['A', 'G'])) else nuc['N3'].pos - nuc.get_com()
get_base_hex_bottom = lambda nuc: nuc['N3'].pos - nuc.get_com() if any((c in nuc.name for c in ['A', 'G'])) else nuc['C2'].pos - nuc.get_com()
get_O3s = lambda nuc: nuc["O3'"].pos - nuc.get_com()
get_C5s = lambda nuc: nuc["C5'"].pos - nuc.get_com()
get_C4s = lambda nuc: nuc["C4'"].pos - nuc.get_com()
get_phosphate = lambda nuc: nuc['P'].pos - nuc.get_com() if not '5' in nuc.name \
                else nuc["HO5'"].pos - nuc.get_com() if "HO5'" in nuc.named_atoms \
                else nuc["O5'"].pos - nuc.get_com() if "O5'" in nuc.named_atoms \
                else np.array([np.nan, np.nan, np.nan])

# Empirically computed best reference points for DNA/RNA
# Voodoo: it works best if you use O3 for generating DNA fragments, but C4 for actual mapping
#DNA_funcs = {"back" : get_O3s, "base" : get_base_center}
DNA_funcs = {"back" : get_C4s, "base" : get_base_center}
RNA_funcs = {"back" : get_C5s, "base" : get_base_edge}

# Don't delete this function!
def choose_reference_nucleotides(nucleotides:List[PDB_Nucleotide]) -> Dict[str, PDB_Nucleotide]:
    """
        Find nucleotides that most look like an oxDNA nucleotide geometry.

        This function is never used in any production code, but it is used for building the reference library by `development code <https://github.com/ErikPoppleton/oxDNA_backmapping>`_. 

        Parameters:
            nucleotides (List[PDB_Nucleotide]) : List of nucleotides to compare.
        
        Returns:
            Dict[str, PDB_Nucleotide] : The best nucleotide for each type in the format `{'C' : PDB_Nucleotide}`.
    """
    ref_a1 = np.array([1, 0, 0])
    ref_a3 = np.array([0, 0, 1])
    ref_a2 = np.array([0, 1, 0])
    bs = utils.get_pos_base(np.zeros(3), ref_a1, ref_a3) * FROM_OXDNA_TO_ANGSTROM

    log("Scoring bases... Lower scores are better.")
    bases = {}
    for n in nucleotides:
        if 'D' in n.name:
            funcs = DNA_funcs
            bbs = utils.get_pos_back(np.zeros(3), ref_a1, ref_a3, type='DNA') * FROM_OXDNA_TO_ANGSTROM
        else:
            # Things that aren't DNA are treated as RNA.
            funcs = RNA_funcs
            bbs = utils.get_pos_back(np.zeros(3), ref_a1, ref_a3, type='RNA') * FROM_OXDNA_TO_ANGSTROM

        # Get the all-atom proxies for the oxDNA sites
        n.compute_as()
        proxies = np.array([
            n.a1,
            n.a3,
            n.a2,
            funcs["back"](n),
            funcs["base"](n)
        ])
        ref = np.array([
            ref_a1,
            ref_a3,
            ref_a2,
            bbs,
            bs
        ])

        # Align the all-atom representation to the oxDNA bead
        utils.kabsch_align(proxies, ref, center=False, inplace=True)

        # Score the alignment and keep the best
        diff = np.mean(np.linalg.norm(proxies - ref, axis=1))
        
        if n.base in bases:
            if diff < bases[n.base].pdb_score: # Find the most oxDNA-like nucleotide for each base type
                bases[n.base] = copy.deepcopy(n)
                bases[n.base].pdb_score = diff
        else:
            bases[n.base] = copy.deepcopy(n)
            bases[n.base].pdb_score = diff

    for k, v in bases.items():
        log(f"Base {k} : best score {v.pdb_score:.3f}")

    return bases

def _atom_from_cif_row(row: Dict) -> Atom:
    """Construct a UTILS/pdb.py Atom from an _atom_site mmCIF row dict."""
    raw = row.get('label_seq_id') or row.get('auth_seq_id', '0')
    return Atom.from_dict(
        name     = (row.get('label_atom_id') or row.get('auth_atom_id', '')).strip(),
        residue  = (row.get('label_comp_id') or row.get('auth_comp_id', '')).strip(),
        chain_id = (row.get('label_asym_id') or row.get('auth_asym_id', '')).strip(),
        residue_idx = int(raw) if raw not in ('.', '?') else 0,
        pos      = np.array([float(row['Cartn_x']), float(row['Cartn_y']), float(row['Cartn_z'])]),
    )


def get_AAs_from_PDB(pdbfile:str, start_res:int=0, n_res:int=-1) -> Tuple[int, List[PDB_AminoAcid]]:
    """
        Get amino acid descriptions from a PDB or mmCIF file.

        Parameters:
            pdbfile (str) : File path to PDB or mmCIF file (format auto-detected).
            start_res (int) : Residue index to start from, default 0.
            n_res (int) : Get only this many residues, default until end of file.

        Returns
            (Tuple[int, List[PDB_AminoAcid]]) : The next residue position (or -1 if
            the file was exhausted) and the list of amino acids extracted.
    """
    with open(pdbfile) as f:
        content = f.read()

    if is_mmcif(content):
        # ── mmCIF path ────────────────────────────────────────────────────────
        rows = [r for r in parse_atom_site(content)
                if r.get('group_PDB', 'ATOM') == 'ATOM']
        amino_acids: List[PDB_AminoAcid] = []
        old_key = None
        aa = None
        for row in rows:
            resn = (row.get('label_comp_id') or row.get('auth_comp_id', '')).strip()
            if resn not in NAME_TO_AA:
                continue   # skip nucleic acids, ligands, water, etc.
            chain = (row.get('label_asym_id') or row.get('auth_asym_id', '')).strip()
            raw   = row.get('label_seq_id') or row.get('auth_seq_id', '0')
            resi  = int(raw) if raw not in ('.', '?') else 0
            key   = (chain, resi)
            if key != old_key:
                aa = PDB_AminoAcid(resn, resi)
                amino_acids.append(aa)
                old_key = key
            aa.add_atom(_atom_from_cif_row(row))
    else:
        # ── PDB path ──────────────────────────────────────────────────────────
        amino_acids = []
        old_residue = ''
        for l in content.splitlines():
            if l.startswith('ATOM'):
                if l[17:20].strip() not in NAME_TO_AA:
                    continue   # skip nucleic acids, ligands, water, etc.
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

    next_pos = -1 if end >= len(amino_acids) else end
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
    
    # Align PDB CA COM to oxDNA bead COM
    ca_poses = np.array([a.get_ca_pos() for a in amino_acids])
    pdb_com  = np.mean(ca_poses, axis=0)
    ox_com   = np.mean(coord,    axis=0)
    centered_ca_poses = ca_poses - pdb_com
    centered_ox_poses = coord    - ox_com
    
    _, R_row = utils.kabsch_align(centered_ca_poses, centered_ox_poses,
                                   center=False, return_rot=True)
    R = R_row.T

    # reposition and rotate each amino acid
    for i, aa in enumerate(amino_acids):
        aa.set_ca_pos(coord[i])
        aa.rotate(R)

    return(reading_position, amino_acids)

_HY36_DIGITS_UPPER = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
_HY36_DIGITS_LOWER = '0123456789abcdefghijklmnopqrstuvwxyz'

def _hy36_encode_pure(digits: str, value: int) -> str:
    n = len(digits)
    result = []
    while value:
        result.append(digits[value % n])
        value //= n
    result.reverse()
    return ''.join(result) if result else digits[0]

def _format_res_seq(n: int) -> str:
    """Format a PDB residue sequence number using the hybrid-36 convention (Grosse-Kunstleve 2007).

    1–9999: decimal.  10000–1223055: uppercase base-36 (A000–ZZZZ).
    1223056–2436111: lowercase base-36 (a000–zzzz).
    """
    if n < 10000:
        return f"{n:4d}"
    i = n - 10000
    if i < 26 * 36**3:
        return _hy36_encode_pure(_HY36_DIGITS_UPPER, i + 10 * 36**3)
    i -= 26 * 36**3
    if i < 26 * 36**3:
        return _hy36_encode_pure(_HY36_DIGITS_LOWER, i + 10 * 36**3)
    raise ValueError(f"residue sequence number {n} out of hybrid-36 range (max 2436111)")

def _format_atom_serial(n: int) -> str:
    """Format a PDB atom serial using the hybrid-36 convention (Grosse-Kunstleve 2007).

    1–99999: decimal.  100000–43770015: uppercase base-36 (A0000–ZZZZZ).
    43770016–87440031: lowercase base-36 (a0000–zzzzz).
    """
    if n < 100000:
        return f"{n:5d}"
    i = n - 100000
    if i < 26 * 36**4:
        return _hy36_encode_pure(_HY36_DIGITS_UPPER, i + 10 * 36**4)
    i -= 26 * 36**4
    if i < 26 * 36**4:
        return _hy36_encode_pure(_HY36_DIGITS_LOWER, i + 10 * 36**4)
    raise ValueError(f"atom serial {n} out of hybrid-36 range (max 87440031)")

def _load_templates() -> tuple:
    """Load DNA and RNA nucleotide fragment templates from pdb_templates/."""
    DNAbases: Dict = {}
    RNAbases: Dict = {}
    for f in PDB_PATH.iterdir():
        base = get_nucs_from_PDB(str(f))[0]
        base.compute_as()
        base.a1, base.a2, base.a3 = utils.get_orthonormalized_base(base.a1, base.a2, base.a3)
        if 'D' in f.stem:
            DNAbases[f.stem] = base
        else:
            RNAbases[f.stem] = base
    return DNAbases, RNAbases


def _preprocess_templates(DNAbases: Dict, RNAbases: Dict, hydrogen: bool) -> Dict:
    """
    Extract numpy arrays from PDB_Nucleotide template objects.

    Called once per oxDNA_PDB invocation.  The returned dict maps each
    template key (e.g. 'DA3', 'RG5') to a sub-dict containing:
      'proxies'       : (5, 3) reference proxy points for Kabsch alignment
      'atoms_centered': (n_atoms, 3) atom positions centred at their COM
      'n_atoms'       : int
      'a1'            : (3,) template a1 vector (used for the post-rotation shift)
      'atom_names'    : list[str]
      'residue_name'  : str
    """
    result: Dict = {}
    for name, base in DNAbases.items():
        result[name] = _extract_template_arrays(base, DNA_funcs, hydrogen)
    for name, base in RNAbases.items():
        result[name] = _extract_template_arrays(base, RNA_funcs, hydrogen)
    return result


def _extract_template_arrays(base, funcs: Dict, hydrogen: bool) -> Dict:
    # The original per-nucleotide loop deep-copied the template then called
    # compute_as(), which overwrites the orthonormalized a1/a2/a3 stored by
    # _load_templates() with the raw ring-atom values.  The proxies must use
    # those same raw values to produce an identical rotation matrix R.
    saved_a1, saved_a2, saved_a3 = base.a1.copy(), base.a2.copy(), base.a3.copy()
    base.compute_as()
    a1 = base.a1.copy()
    a3 = base.a3.copy()
    a2 = base.a2.copy()
    proxies = np.array([a1, a3, a2, funcs["back"](base), funcs["base"](base)], dtype=float)
    base.a1, base.a2, base.a3 = saved_a1, saved_a2, saved_a3   # restore

    # Original code centred on the COM of ALL atoms (including H) before
    # rotating, even when hydrogen=False.  Replicate that here.
    all_atoms = base.get_atoms()
    all_positions = np.array([a.pos for a in all_atoms], dtype=float)
    com = all_positions.mean(axis=0)

    out_atoms = [a for a in all_atoms if hydrogen or 'H' not in a.name]
    out_positions = np.array([a.pos for a in out_atoms], dtype=float)
    atoms_centered = out_positions - com

    return {
        'proxies':        proxies,
        'atoms_centered': atoms_centered,
        'n_atoms':        len(out_atoms),
        'a1':             a1,           # raw (non-orthonormalized), matches loop
        'atom_names':     [a.name for a in out_atoms],
        'residue_name':   base.name,
    }


def _nucleotide_template_key(nucleotide, strand, isDNA: bool, uniform_residue_names: bool) -> str:
    """Return the template dict key for a single nucleotide."""
    if isinstance(nucleotide.btype, int):
        nb = number_to_DNAbase[nucleotide.btype % 4] if isDNA else number_to_RNAbase[nucleotide.btype % 4]
    elif isinstance(nucleotide.btype, str):
        nb = nucleotide.btype
    else:
        raise RuntimeError(f"Bad base type: {nucleotide.btype} on nucleotide id {nucleotide.id}")

    if uniform_residue_names:
        end_type = ""
    elif (nucleotide is strand.monomers[0] or nucleotide is strand.monomers[-1]) and not strand.is_circular():
        if strand.is_old():
            end_type = "3" if nucleotide is strand.monomers[0] else "5"
        else:
            end_type = "5" if nucleotide is strand.monomers[0] else "3"
    else:
        end_type = ""

    sugar_type = 'D' if strand.type == 'DNA' else ('R' if strand.type == 'RNA' and end_type else '')
    return sugar_type + nb + end_type


def _batch_kabsch(mobile: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """
    Vectorised Kabsch rotation matrices for N pairs of point sets.

    Matches ``utils.kabsch_align`` exactly:
        cov = mobile.T @ ref
        u, _, vt = svd(cov)
        rot = (vt.T @ u.T).T  =  u @ vt
        if det(rot) < 0: vt[2] = -vt[2]; rot = u @ vt  (reflection correction)

    Parameters
    ----------
    mobile : (N, K, 3)  — the "from" point clouds
    ref    : (N, K, 3)  — the "to"   point clouds

    Returns
    -------
    R : (N, 3, 3) such that ``mobile @ R ≈ ref`` for each sample.
    """
    M = mobile.swapaxes(-1, -2) @ ref   # (N, 3, 3): M[n] = mobile[n].T @ ref[n]
    U, _s, Vt = np.linalg.svd(M)       # each (N, 3, 3); numpy returns Vt
    # Reflection correction: flip last row of Vt when det(U @ Vt) < 0
    d = np.linalg.det(U @ Vt)          # (N,)
    D = np.zeros((len(M), 3, 3))
    D[:, 0, 0] = D[:, 1, 1] = 1.0
    D[:, 2, 2] = d
    return U @ D @ Vt                   # (N, 3, 3)


def _build_nucleic_strand_atoms(
    strand,
    conf: Configuration,
    templates: Dict,
    rmsf_per_nucleotide,
    box_angstrom: np.ndarray,
    reverse: bool,
    uniform_residue_names: bool,
) -> List[List[Dict]]:
    """
    Convert one oxDNA DNA/RNA strand into a list of per-residue atom dicts.

    Uses vectorised numpy operations: all Kabsch alignments and atom
    rotations for the strand are computed in a single batch, with no
    Python loop over nucleotides during the heavy linear-algebra work.
    """
    monomers = strand.monomers
    N = len(monomers)
    if N == 0:
        return []

    isDNA = strand.get_kwdata()['type'] == 'DNA'

    # ── 1. Resolve template key for every nucleotide ──────────────────────────
    nb_list = [_nucleotide_template_key(nuc, strand, isDNA, uniform_residue_names)
               for nuc in monomers]

    # ── 2. Allocate batch arrays ──────────────────────────────────────────────
    max_atoms = max(templates[nb]['n_atoms'] for nb in nb_list)

    proxies_batch  = np.empty((N, 5, 3))          # template reference frames
    ox_sites_batch = np.empty((N, 5, 3))          # target oxDNA reference frames
    atoms_padded   = np.zeros((N, max_atoms, 3))  # zero-padded centered positions
    a1_tpl         = np.empty((N, 3))             # template a1 for post-rotation shift
    pos_batch      = np.empty((N, 3))             # target positions in Å
    n_atoms_arr    = np.empty(N, dtype=int)

    for i, (nuc, nb) in enumerate(zip(monomers, nb_list)):
        tpl = templates[nb]
        n   = tpl['n_atoms']

        proxies_batch[i]     = tpl['proxies']
        atoms_padded[i, :n]  = tpl['atoms_centered']
        a1_tpl[i]            = tpl['a1']
        n_atoms_arr[i]       = n

        pos = conf.positions[nuc.id] * FROM_OXDNA_TO_ANGSTROM
        a1  = conf.a1s[nuc.id]
        a3  = conf.a3s[nuc.id]
        a2  = np.cross(a3, a1)
        bbs = (utils.get_pos_back(pos, a1, a3, type=strand.type) - pos) * FROM_OXDNA_TO_ANGSTROM
        bs  = (utils.get_pos_base( pos, a1, a3, type=strand.type) - pos) * FROM_OXDNA_TO_ANGSTROM
        ox_sites_batch[i] = (a1, a3, a2, bbs, bs)
        pos_batch[i]      = pos

    # ── 3. Batch Kabsch: one SVD call for the whole strand ────────────────────
    R_batch = _batch_kabsch(proxies_batch, ox_sites_batch)  # (N, 3, 3)

    # ── 4. Rotate and translate all atom clouds simultaneously ─────────────────
    # atoms_padded @ R_batch : (N, max_atoms, 3) @ (N, 3, 3) → (N, max_atoms, 3)
    final_atoms = atoms_padded @ R_batch + pos_batch[:, None, :]

    # ── 5. a1 shift (empirical correction that slightly improves RNA geometry) ─
    new_a1 = (a1_tpl[:, None, :] @ R_batch).squeeze(1)  # (N, 3)
    final_atoms -= 0.5 * new_a1[:, None, :]

    # ── 6. Reconstruct the List[List[Dict]] output format ────────────────────
    strand_pdb = []
    for i, (nuc, nb) in enumerate(zip(monomers, nb_list)):
        tpl  = templates[nb]
        n    = tpl['n_atoms']
        bfac = rmsf_per_nucleotide[nuc.id]
        strand_pdb.append([
            {'name': tpl['atom_names'][j], 'residue_name': tpl['residue_name'],
             'pos': final_atoms[i, j], 'bfactor': bfac}
            for j in range(n)
        ])

    if reverse:
        strand_pdb = strand_pdb[::-1]

    return strand_pdb


def _strand_sequence(strand) -> str:
    """Return the one-letter sequence of a DNA/RNA strand in storage order."""
    isDNA = strand.get_kwdata()['type'] == 'DNA'
    bases = []
    for m in strand.monomers:
        if type(m.btype) == int:
            bases.append(number_to_DNAbase[m.btype % 4] if isDNA else number_to_RNAbase[m.btype % 4])
        else:
            bases.append(str(m.btype))
    return ''.join(bases)


def _resolve_output_format(format: str, system: System, hydrogen: bool = True) -> str:
    """Return ``'pdb'`` or ``'mmcif'`` for the given format hint and system size.

    ``'auto'`` selects mmCIF when the estimated atom count exceeds 99 999 or
    the total residue count exceeds 9 999 — the plain-decimal PDB limits.
    """
    if format == 'mmcif':
        return 'mmcif'
    if format == 'pdb':
        return 'pdb'
    total_residues = sum(len(s.monomers) for s in system.strands)
    atoms_per_nuc = 33 if hydrogen else 22   # conservative upper bound
    est_atoms = total_residues * atoms_per_nuc
    if est_atoms > 99999 or total_residues > 9999:
        log(f"System has ~{est_atoms} estimated atoms / {total_residues} residues; "
            "automatically using mmCIF output")
        return 'mmcif'
    return 'pdb'


def _chain_id_generator(mmcif: bool) -> Iterator[str]:
    """Yield unique chain IDs for each strand.

    For PDB output the IDs are single characters (A–Z, a–z, 0–9); a warning
    is logged if the pool of 62 is exhausted and IDs start cycling.
    For mmCIF output multi-character IDs are used once the 26 uppercase
    single-letter pool is exhausted: A–Z, AA–AZ, BA–BZ, …, ZA–ZZ, AAA–AAZ, …
    """
    if mmcif:
        for length in itertools.count(1):
            for combo in itertools.product('ABCDEFGHIJKLMNOPQRSTUVWXYZ', repeat=length):
                yield ''.join(combo)
    else:
        chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
        i = 0
        while True:
            if i > 0 and i % len(chars) == 0:
                log("More than 62 chains identified, looping chain identifier...", level='warning')
            yield chars[i % len(chars)]
            i += 1


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
            print("{:6s}{:5s} {:^4s}{:1s}{:>3s} {:1s}{:4s}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
                .format(
                    "ATOM",                         #record
                    _format_atom_serial(atom_counter), #atom_id
                    a['name'],                      #atom_name
                    " ",                            #alt_loc
                    a['residue_name'],              #res_name
                    chain_id,                       #chain_id
                    _format_res_seq(nid),           #res_id
                    " ",                            #ins_code
                    a['pos'][0],                    #coord_x
                    a['pos'][1],                    #coord_y
                    a['pos'][2],                    #coord_z
                    1.00,                           #residency
                    a['bfactor'],                   #b-factor
                    " ", " "                        #element,charge
                ),
                file=out
            )
            atom_counter = (atom_counter + 1) % 87440032  # hybrid36 max: zzzzz
    print("TER", file=out)

    return(atom_counter)

def oxDNA_PDB(conf:Configuration, system:System, out_basename:str, protein_pdb_files:Union[List[str], None]=None, reverse:bool=False, hydrogen:bool=True, uniform_residue_names:bool=False, one_file_per_strand:bool=False, rmsf_file:str='', format:str='auto'):
    """
        Convert an oxDNA file to a PDB or mmCIF file using fragment assembly.

        Parameters:
            conf (Configuration) : The Configuration to convert
            system (System) : The system topology for the configuration to convert
            out_basename (str) : Filename base (without extension) to write out to
            protein_pdb_files (List[str]) : Filenames of pdb files corresponding to proteins present in the oxDNA files. (Default: [])
            reverse (bool) :  Reverse nucleic acid strands from oxDNA files (If you want 'correct' PDB files from backwards oxDNA files). (Default: False)
            hydrogen (bool) : Write hydrogens to output file (Default: True)
            uniform_residue_names (bool) : Don't add '5' and '3' to the names of the terminal residues (Default: False)
            one_file_per_strand (bool) : Split each strand into a separate file (appends strand id to out_basename). (Default: False)
            rmsf_file (str) : Write oxDNA (json-formatted from deviations) RMSFs into the b-factor field. (Default: '')
            format (str) : Output format: 'pdb', 'mmcif', or 'auto'. 'auto' selects mmCIF when the
                system exceeds the plain-decimal PDB limits (>99999 atoms or >9999 residues). (Default: 'auto')
    """
    use_mmcif = _resolve_output_format(format, system, hydrogen) == 'mmcif'
    ext = '.cif' if use_mmcif else '.pdb'

    DNAbases, RNAbases = _load_templates()
    templates = _preprocess_templates(DNAbases, RNAbases, hydrogen)
    box_angstrom = conf.box * FROM_OXDNA_TO_ANGSTROM

    # Handle RMSF -> bFactor conversion
    if rmsf_file:
        with open(rmsf_file) as f:
            try:
                substrings = f.read().split("[")[1].split("]")[0].split(",")
            except Exception as e:
                raise RuntimeError("Parsing error in RMSF file. Invalid Format: %s" % e)
            try:
                rmsf_per_nucleotide = {i: float(s) for i, s in enumerate(substrings)}
            except Exception as e:
                raise RuntimeError("Parsing error in RMSF file. Conversion to float failed : %s" % e)
    else:
        rmsf_per_nucleotide = defaultdict(lambda: 1.00)

    if one_file_per_strand:
        out_name = out_basename + "_{}".format(system.strands[0].id) + ext
    else:
        out_name = out_basename + ext

    # One PDB file per protein strand; create the iterator once so each strand
    # advances to the next file rather than re-reading from the first.
    protein_file_iter = iter(protein_pdb_files or [])

    with open(out_name, 'w+') as out:
        chain_iter = _chain_id_generator(use_mmcif)
        chain_id = next(chain_iter)
        atom_counter = 1
        entity_id = 1

        # mmCIF needs all strand data before writing entity blocks (combined file).
        all_chain_data = []
        entity_rows = []
        entity_poly_rows = []

        if use_mmcif and not one_file_per_strand:
            write_data_block_header(out, os.path.basename(out_basename))

        for i, strand in enumerate(system.strands):
            strand_pdb = []
            log("Converting strand {}".format(strand.id), end='\r')

            if strand.type == 'peptide':
                pdbfile = next(protein_file_iter, None)
                if pdbfile is None:
                    raise RuntimeError("You must provide a PDB file for each protein strand.")
                _, amino_acids = peptide_to_pdb(strand, conf, pdbfile, 0)
                for aa in amino_acids:
                    strand_pdb.append(aa.to_pdb(hydrogen, bfactor=rmsf_per_nucleotide[aa.idx]))

            elif strand.id < 0:
                raise RuntimeError("You must provide PDB files containing just the protein for each protein in the scene.")

            elif strand.type == 'DNA' or strand.type == 'RNA':
                strand_pdb = _build_nucleic_strand_atoms(
                    strand, conf, templates,
                    rmsf_per_nucleotide,
                    box_angstrom, reverse, uniform_residue_names,
                )

            else:
                log(f"Unknown strand type {strand.type} on strand {strand.id}. Skipping", level='warning')

            if strand_pdb:
                if use_mmcif:
                    if strand.type in ('DNA', 'RNA'):
                        poly_type = 'polydeoxyribonucleotide' if strand.type == 'DNA' else 'polyribonucleotide'
                        entity_rows.append((entity_id, f'{strand.type} strand'))
                        entity_poly_rows.append((entity_id, poly_type, _strand_sequence(strand)))
                    else:
                        entity_rows.append((entity_id, 'protein strand'))

                    if one_file_per_strand:
                        write_data_block_header(out, f'{os.path.basename(out_basename)}_{strand.id}')
                        write_entity_block(entity_rows[-1:], out)
                        if entity_poly_rows and entity_poly_rows[-1][0] == entity_id:
                            write_entity_poly_block(entity_poly_rows[-1:], out)
                        write_atom_site_block([(chain_id, strand_pdb, entity_id)], out)
                        entity_rows.clear()
                        entity_poly_rows.clear()
                    else:
                        all_chain_data.append((chain_id, strand_pdb, entity_id))
                else:
                    atom_counter = write_strand_to_PDB(strand_pdb, chain_id, atom_counter, out)

            if one_file_per_strand:
                if not use_mmcif:
                    print("END", file=out)
                out.close()
                log("Wrote strand {}'s data to {}".format(strand.id, out_name))
                chain_iter = _chain_id_generator(use_mmcif)
                chain_id = next(chain_iter)
                entity_id = 1
                if i < len(system) - 1:
                    next_strand = system.strands[i + 1]
                    out_name = out_basename + "_{}".format(next_strand.id) + ext
                    out = open(out_name, "w")
            else:
                chain_id = next(chain_iter)
                entity_id += 1

        if not one_file_per_strand:
            if use_mmcif and all_chain_data:
                write_entity_block(entity_rows, out)
                if entity_poly_rows:
                    write_entity_poly_block(entity_poly_rows, out)
                write_atom_site_block(all_chain_data, out)
            elif not use_mmcif:
                print("END", file=out)
        print()

    log("Wrote data to '{}'".format(out_name))
    log("DONE")


def cli_parser(prog="oxDNA_PDB.py"):
    parser = argparse.ArgumentParser(prog=prog, description="Convert oxDNA files to PDB or mmCIF.  This converter can handle oxDNANM protein simulation files.")
    parser.add_argument('topology', type=str,
                        help='the oxDNA topology file for the structure')
    parser.add_argument('configuration', type=str,
                        help='the configuration file you wish to convert')
    parser.add_argument('direction', type=str,
                        help='the direction of strands in the oxDNA files, either 35 or 53.  Old oxDNA files are 3-5. New oxDNA files are 5-3')
    parser.add_argument('pdbfiles', type=str, nargs='*',
                        help='PDB files for the proteins present in your structure.  The strands in the PDB file(s) must be in the same order as your oxDNA file. If there are multiple of the same protein, you must provide that PDB file that many times.')
    parser.add_argument('-o', '--output', type=str, 
                        help='The name of the output pdb file.  Defaults to name of the configuration+.pdb')
    parser.add_argument('-d', '--output_direction', type=str,
                        help='Direction to save nucleic acid strands in.  Should be "53" or "35".  Defaults to same as input direction.')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, 
                        help="Don't print 'INFO' messages to stderr")
    parser.add_argument('-H', '--hydrogen', action='store_true', default=True,
                        help='if you want to include hydrogen atoms in the output PDB file')
    parser.add_argument('-u', '--uniform-residue-names', action='store_true', default=False,
                        help='if you want to use uniform residue names in the output PDB file')
    parser.add_argument('-1', '--one_file_per_strand', action='store_true',
                        default=False, help='if you want to have one PDB file per strand')
    parser.add_argument('-r', '--rmsf-file', dest='rmsf_bfactor', type=str, nargs=1,
                        help='A RMSF file from deviations.  Will be used to fill the b-factors field in the output PDB (only for D(R)NA)')
    parser.add_argument('--format', choices=['auto', 'pdb', 'mmcif'], default='auto',
                        help="Output format: 'pdb' (always PDB), 'mmcif' (always mmCIF), "
                             "'auto' (mmCIF when system exceeds plain-decimal PDB limits). Default: auto")
    parser.add_argument('--no-inbox', dest='no_inbox', action='store_true', default=False,
                        help='Skip the inbox (periodic-boundary centering) step before conversion. '
                             'Use when the structure is already centred or spans a periodic box.')
    return parser

def main():
    parser = cli_parser(os.path.basename(__file__))
    args = parser.parse_args()

    # Parse positional arguments
    logger_settings.set_quiet(args.quiet)
    top_file = args.topology
    conf_file = args.configuration
    direction = args.direction
    if direction not in ["35", "53"]:
        raise RuntimeError("Error: Direction must be either 35 or 53")
    if args.pdbfiles:
        protein_pdb_files = args.pdbfiles
    else:
        protein_pdb_files = []

    # Parse optional arguments
    fmt = args.format
    if args.output:
        if args.output.endswith('.cif'):
            out_basename = args.output[:-4]
            if fmt == 'auto':
                fmt = 'mmcif'
        else:
            out_basename = re.sub(r"\.pdb$", "", args.output)
    else:
        out_basename = os.path.splitext(conf_file)[0]
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

    # Read oxDNA configuration
    system, _ = strand_describe(top_file)
    ti, di = describe(top_file, conf_file)
    conf = get_confs(ti, di, 0, 1)[0]
    if not args.no_inbox:
        conf = inbox(conf, center=True)

    oxDNA_PDB(conf, system, out_basename, protein_pdb_files, reverse, hydrogen, uniform_residue_names, one_file_per_strand, rmsf_file, format=fmt)

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()
