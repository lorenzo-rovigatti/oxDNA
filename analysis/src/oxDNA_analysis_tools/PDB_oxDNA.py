import time
start_time = time.time()
import argparse
import numpy as np
from os import path
from typing import Tuple, List, Dict
from dataclasses import dataclass
from itertools import permutations
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.mmcif import parse_atom_site, is_mmcif
from oxDNA_analysis_tools.UTILS.pdb import (
    NAME_TO_AA,
    FROM_OXDNA_TO_ANGSTROM as OXDNA_TO_ANGSTROM,
    FROM_ANGSTROM_TO_OXDNA as ANGSTROM_TO_OXDNA,
)
from oxDNA_analysis_tools.UTILS.RyeReader import write_conf, write_top
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, System, Strand, Monomer

# ── Vectorised geometry helpers ───────────────────────────────────────────────

# Pre-computed at import time — same for every call.
_RING_NAMES = ("C2", "C4", "C5", "C6", "N1", "N3")
_RING_IDX   = {name: i for i, name in enumerate(_RING_NAMES)}
# All 120 ordered triples drawn from the 6 ring atoms (same as permutations(6,3))
_PERM_ARR   = np.array(list(permutations(range(6), 3)), dtype=np.intp)  # (120, 3)
# Index pairs for a1: pyrimidines (C/T/U) and purines (A/G)
_PYR_PAIRS  = np.array([[_RING_IDX["N3"], _RING_IDX["C6"]],
                         [_RING_IDX["C2"], _RING_IDX["N1"]],
                         [_RING_IDX["C4"], _RING_IDX["C5"]]], dtype=np.intp)  # (3, 2)
_PUR_PAIRS  = np.array([[_RING_IDX["N1"], _RING_IDX["C4"]],
                         [_RING_IDX["C2"], _RING_IDX["N3"]],
                         [_RING_IDX["C6"], _RING_IDX["C5"]]], dtype=np.intp)  # (3, 2)


def _calc_ox_properties_batch(residues: list):
    """
    Vectorised equivalent of ``Residue.calc_ox_properties()`` for a list of
    residues belonging to the same strand.

    All 120 ring-atom permutations (for a3) and the 3 base-edge pairs (for a1)
    are evaluated across every residue simultaneously using numpy batch ops,
    replacing the per-residue Python loops that dominate runtime for large
    structures.

    Returns
    -------
    positions : (N, 3) — COM in oxDNA units
    a1s       : (N, 3) — base-pairing-edge unit vectors
    a3s       : (N, 3) — base-normal unit vectors
    """
    for r in residues:
        r.parse_nucleic_atoms()   # populate atom_lookup if not done

    # ── Ring-atom positions: (N, 6, 3) ────────────────────────────────────
    ring_atoms = np.array(
        [[r.atom_lookup[name].pos for name in _RING_NAMES] for r in residues],
        dtype=float,
    )

    # ── Reference vector for a3 sign: C5' − C3', (N, 3) ──────────────────
    ref_vecs = np.array(
        [r.atom_lookup["C5'"].pos - r.atom_lookup["C3'"].pos for r in residues],
        dtype=float,
    )

    # ── COM of all atoms, converted to oxDNA units: (N, 3) ────────────────
    positions = np.array(
        [np.mean([a.pos for a in r.atoms], axis=0) for r in residues],
        dtype=float,
    ) * ANGSTROM_TO_OXDNA

    # ── Batch a3 ──────────────────────────────────────────────────────────
    # Gather the three point clouds for all 120 permutations: each (N, 120, 3)
    p    = ring_atoms[:, _PERM_ARR[:, 0], :]
    q    = ring_atoms[:, _PERM_ARR[:, 1], :]
    r_pt = ring_atoms[:, _PERM_ARR[:, 2], :]

    v1 = p - q
    v2 = p - r_pt
    v1 /= np.linalg.norm(v1, axis=-1, keepdims=True)
    v2 /= np.linalg.norm(v2, axis=-1, keepdims=True)

    a3_all = np.cross(v1, v2)                                # (N, 120, 3)
    a3_all /= np.linalg.norm(a3_all, axis=-1, keepdims=True)

    # Flip contributions where the normal points away from the helix axis
    dots = np.einsum('nki,ni->nk', a3_all, ref_vecs)        # (N, 120)
    a3_all *= np.where(dots < 0, -1.0, 1.0)[:, :, None]

    a3 = a3_all.sum(axis=1)                                  # (N, 3)
    a3 /= np.linalg.norm(a3, axis=-1, keepdims=True)

    # ── Batch a1 ──────────────────────────────────────────────────────────
    is_pyr = np.array(
        [("C" in r.resn or "T" in r.resn or "U" in r.resn) for r in residues],
        dtype=bool,
    )

    def _pair_a1(pair_idx):
        p = ring_atoms[:, pair_idx[:, 0], :]  # (N, 3, 3)
        q = ring_atoms[:, pair_idx[:, 1], :]
        vecs = p - q
        vecs /= np.linalg.norm(vecs, axis=-1, keepdims=True)
        return vecs.sum(axis=1)               # (N, 3)

    pyr_a1 = _pair_a1(_PYR_PAIRS)
    pur_a1 = _pair_a1(_PUR_PAIRS)

    a1 = np.where(is_pyr[:, None], pyr_a1, pur_a1)          # (N, 3)
    a1 /= np.linalg.norm(a1, axis=-1, keepdims=True)

    return positions, a1, a3


def _flush_to_strand(residue_buffer: list, strand) -> None:
    """Batch-compute oxDNA properties for all buffered residues and append to strand."""
    if not residue_buffer:
        return
    if isinstance(residue_buffer[0], ProteinResidue):
        for r_b in residue_buffer:
            strand.append(r_b.to_monomer(strand))
    else:
        positions, a1s, a3s = _calc_ox_properties_batch(residue_buffer)
        for r_b, pos, a1, a3 in zip(residue_buffer, positions, a1s, a3s):
            strand.append(Monomer(r_b.resi, r_b.resn, strand, None, None, None, pos, a1, a3))
    residue_buffer.clear()


def _is_amino_acid(resn: str) -> bool:
    """Return True if *resn* is a standard amino acid 3-letter code."""
    return resn.strip() in NAME_TO_AA


class ProteinResidue:
    """Minimal all-atom amino acid residue for coarse-graining to an oxDNA ANM bead.

    The CG position is the Cα position in oxDNA units; a1/a3 are set to the
    identity frame because protein beads are isotropic.
    """

    def __init__(self, resn: str, resi: int):
        self.resn = resn.strip()
        self.resi = resi
        self.ca_pos: np.ndarray = None
        self.atoms: list = []

    def add_atom(self, atom) -> None:
        self.atoms.append(atom)
        if atom.type == 'CA':
            self.ca_pos = atom.pos.copy()

    def is_valid(self) -> bool:
        return self.ca_pos is not None

    def to_monomer(self, strand) -> Monomer:
        # Use the one-letter code so topology writers can join btypes as strings.
        btype = NAME_TO_AA.get(self.resn, 'X')
        pos   = self.ca_pos * ANGSTROM_TO_OXDNA
        return Monomer(self.resi, btype, strand, None, None, None,
                       pos, np.array([1., 0., 0.]), np.array([0., 0., 1.]))

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
        self.resn = resn.strip('DR35*')
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
    
    def has_backbone_atoms(self) -> bool:
        atom_names = {a.type for a in self.atoms}
        return "C3'" in atom_names and "C5'" in atom_names

    def to_monomer(self, strand:Strand) -> Monomer:
        p, a1, a3 = self.calc_ox_properties()
        return Monomer(self.resi, self.resn, strand, None, None, None, p, a1, a3)

_HY36_UPPER_VALUES = {c: i for i, c in enumerate('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ')}
_HY36_LOWER_VALUES = {c: i for i, c in enumerate('0123456789abcdefghijklmnopqrstuvwxyz')}

def _hy36_decode_pure(values: dict, s: str) -> int:
    result = 0
    for c in s:
        result = result * 36 + values[c]
    return result

def _parse_res_seq(s: str) -> int:
    """Parse a PDB residue sequence number using the hybrid-36 convention (Grosse-Kunstleve 2007).

    Decimal 1–9999, uppercase A000–ZZZZ (10000–1223055), lowercase
    a000–zzzz (1223056–2436111).
    """
    s = s.strip()
    if not s:
        return 0
    f = s[0]
    if f == '-' or f.isdigit():
        return int(s)
    if f.isupper():
        return _hy36_decode_pure(_HY36_UPPER_VALUES, s) - 10 * 36**3 + 10000
    if f.islower():
        return _hy36_decode_pure(_HY36_LOWER_VALUES, s) + 16 * 36**3 + 10000
    raise ValueError(f"invalid hybrid-36 residue sequence: {s!r}")

def _parse_atom_serial(s: str) -> int:
    """Parse a PDB atom serial using the hybrid-36 convention (Grosse-Kunstleve 2007).

    Decimal 1–99999, uppercase A0000–ZZZZZ (100000–43770015), lowercase
    a0000–zzzzz (43770016–87440031).
    """
    s = s.strip()
    if not s:
        return 0
    f = s[0]
    if f == '-' or f.isdigit():
        return int(s)
    if f.isupper():
        return _hy36_decode_pure(_HY36_UPPER_VALUES, s) - 10 * 36**4 + 100000
    if f.islower():
        return _hy36_decode_pure(_HY36_LOWER_VALUES, s) + 16 * 36**4 + 100000
    raise ValueError(f"invalid hybrid-36 atom serial: {s!r}")

def parse_atom(l:str):
    return Atom(
        _parse_atom_serial(l[6:11]),                                   # atom index
        l[12:16].strip().replace('*', "'"),                            # atom name (normalize * to ' for sugar atoms)
        l[16].strip(),                                                 # alternate location
        l[17:20].strip(),                                              # resname
        l[21],                                                         # chain
        _parse_res_seq(l[22:26]),                                      # resid
        np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])  # xyz coords
    )

_is_mmcif = is_mmcif   # local alias; public version lives in UTILS/mmcif.py


def _finalize_residue(r, residue_buffer: list, strand) -> None:
    """Validate and buffer a completed residue. No-op when *r* is None."""
    if r is None:
        return
    if isinstance(r, ProteinResidue):
        if r.is_valid():
            residue_buffer.append(r)
        else:
            log(f"Residue {r.resn}{r.resi} has no CA atom, skipping", level='warning')
    else:
        if any(atom.type == "O2'" for atom in r.atoms):
            strand.type = 'RNA'
        if not r.has_backbone_atoms():
            log(f"Residue {r.resn}{r.resi} is missing C3' or C5' atom, skipping", level='warning')
        else:
            residue_buffer.append(r)


def _finalize_system(sys: System, first_strand_type: str,
                     systems: list, configs: list) -> System:
    """Separate, relabel, reorder strands; build Configuration; return a fresh System."""
    protein_strands = [s for s in sys.strands if s.type == 'peptide']
    na_strands      = [s for s in sys.strands if s.type != 'peptide']
    for i, s in enumerate(protein_strands, 1):
        s.id = -i
    for i, s in enumerate(na_strands, 1):
        s.id = i
    sys.strands = (protein_strands + na_strands
                   if first_strand_type == 'peptide'
                   else na_strands + protein_strands)
    if not sys.strands:
        return System('', [])
    # Reassign monomer IDs to 0-based global indices so that
    # conf.positions[monomer.id] always resolves to the correct slot.
    # PDB/CIF-loaded monomers carry m.id = PDB residue number by default,
    # which would cause out-of-bounds or wrong-position access downstream.
    global_idx = 0
    for s in sys.strands:
        for m in s.monomers:
            m.id = global_idx
            global_idx += 1
    positions = np.array([m.pos for s in sys.strands for m in s])
    box = 1.5 * (np.max(positions) - np.min(positions))
    configs.append(Configuration(
        0, np.array([box, box, box]), np.array([0, 0, 0]),
        positions,
        np.array([m.a1 for s in sys.strands for m in s]),
        np.array([m.a3 for s in sys.strands for m in s]),
    ))
    systems.append(sys)
    return System('', [])


def _parse_mmcif(cif_str: str) -> Tuple[List[Configuration], List[System]]:
    """Convert an mmCIF string to oxDNA Configuration/System objects."""

    residue_buffer: list = []
    first_strand_type: str = None   # 'peptide' or 'na'

    def end_strand():
        nonlocal strand, sys, prev_chain, prev_resi
        _finalize_residue(r, residue_buffer, strand)
        _flush_to_strand(residue_buffer, strand)
        if len(strand) > 0:
            sys.append(strand)
        strand = Strand(strand.id + 1)
        prev_chain = a.chain
        prev_resi = -1

    def end_system():
        nonlocal sys
        sys = _finalize_system(sys, first_strand_type, systems, configs)

    rows = parse_atom_site(cif_str)

    systems: List[System] = []
    configs: List[Configuration] = []
    sys = System('')
    strand = Strand(1)
    prev_resi: int = -1
    prev_chain: str = ''
    prev_model: str = '1'
    r   = None
    a   = None

    for row in rows:
        # Match PDB behaviour: only process standard ATOM records, skip HETATM
        if row.get('group_PDB', 'ATOM') != 'ATOM':
            continue

        atom_name = (row.get('label_atom_id') or row.get('auth_atom_id', '')).replace('*', "'")
        alt_loc   = row.get('label_alt_id', '.')
        alt_loc   = '' if alt_loc in ('.', '?') else alt_loc
        resn      = row.get('label_comp_id') or row.get('auth_comp_id', '')
        chain     = row.get('label_asym_id') or row.get('auth_asym_id', '')
        raw_seq   = row.get('label_seq_id') or row.get('auth_seq_id', '0')
        resi      = int(raw_seq) if raw_seq not in ('.', '?') else 0
        model     = row.get('pdbx_PDB_model_num', '1')

        if model != prev_model and len(sys.strands) > 0:
            if prev_resi != -1:
                end_strand()
            end_system()
            strand = Strand(1)
            prev_resi = -1
            prev_chain = ''
        prev_model = model

        a = Atom(
            id=int(row.get('id', 0)),
            type=atom_name, alt=alt_loc, resn=resn, chain=chain, resi=resi,
            pos=np.array([float(row['Cartn_x']), float(row['Cartn_y']), float(row['Cartn_z'])]),
        )

        if a.alt and a.alt not in ('A', '1'):
            continue
        if a.alt:
            log(f"Alternate location for atom {a.id} of residue {a.resi} "
                f"encountered, using location {a.alt}")

        if prev_chain != '' and a.chain != prev_chain and prev_resi != -1:
            end_strand()

        if a.resi != prev_resi:
            if prev_resi != -1:
                _finalize_residue(r, residue_buffer, strand)
            if _is_amino_acid(a.resn):
                r = ProteinResidue(a.resn, a.resi)
                strand.type = 'peptide'
                if first_strand_type is None:
                    first_strand_type = 'peptide'
            else:
                try:
                    r = Residue(a.resn, a.resi)
                except RuntimeError as e:
                    log(str(e), level='warning')
                    r = None
                    prev_resi = a.resi
                    prev_chain = a.chain
                    continue
                if first_strand_type is None:
                    first_strand_type = 'na'

        if r is None:
            continue
        if isinstance(r, ProteinResidue):
            r.add_atom(a)
        else:
            r.atoms.append(a)
        prev_resi = a.resi
        prev_chain = a.chain

    if prev_resi != -1:
        _finalize_residue(r, residue_buffer, strand)
        _flush_to_strand(residue_buffer, strand)
        sys.append(strand)
    if len(sys.strands) > 0:
        end_system()

    return configs, systems


def PDB_oxDNA(pdb_str:str) -> Tuple[List[Configuration], List[System]]:
    """
    Convert a PDB or mmCIF string to an oxDNA Configuration/System.

    The format is detected automatically: strings that begin with a ``data_``
    block are treated as mmCIF; everything else is treated as PDB.

    Parameters:
        pdb_str (str) : The contents of a PDB or mmCIF file as a string

    Returns:
        (tuple(List[Configuration], List[System])) : The system and configuration in oxDNA ready for write-out
    """
    if _is_mmcif(pdb_str):
        return _parse_mmcif(pdb_str)

    residue_buffer: list = []
    first_strand_type: str = None   # 'peptide' or 'na'

    def end_strand():
        nonlocal strand, sys, prev_chain, prev_resi
        _finalize_residue(r, residue_buffer, strand)
        _flush_to_strand(residue_buffer, strand)
        if len(strand) > 0:
            sys.append(strand)
        strand = Strand(strand.id+1)
        prev_chain = a.chain
        prev_resi = -1

    def end_system():
        nonlocal sys
        sys = _finalize_system(sys, first_strand_type, systems, configs)

    systems = []
    configs = []
    sys = System('')
    strand = Strand(1)
    prev_resi: int = -1
    prev_chain: str = ''
    r = None

    lines = pdb_str.split('\n')
    for l in lines:
        if l.startswith('ATOM'):
            a = parse_atom(l)

            if prev_chain != '' and a.chain != prev_chain and prev_resi != -1:
                end_strand()

            if a.alt:
                if not (a.alt == 'A' or a.alt == '1'):
                    continue
                log(f"Alternate location for atom {a.id} of residue {a.resi} encountered, using location {a.alt}")

            if a.resi != prev_resi:
                if prev_resi != -1:
                    _finalize_residue(r, residue_buffer, strand)
                if _is_amino_acid(a.resn):
                    r = ProteinResidue(a.resn, a.resi)
                    strand.type = 'peptide'
                    if first_strand_type is None:
                        first_strand_type = 'peptide'
                else:
                    try:
                        r = Residue(a.resn, a.resi)
                    except RuntimeError as e:
                        log(str(e), level='warning')
                        r = None
                        prev_resi = a.resi
                        prev_chain = a.chain
                        continue
                    if first_strand_type is None:
                        first_strand_type = 'na'

            if r is None:
                continue   # unrecognised residue type; skip
            if isinstance(r, ProteinResidue):
                r.add_atom(a)
            else:
                r.atoms.append(a)
            prev_resi = a.resi
            prev_chain = a.chain
            continue

        elif l.startswith('TER'):
            end_strand()
            continue

        elif l.startswith('END'):
            if prev_resi != -1:
                end_strand()
            if len(sys.strands) > 0:
                end_system()
            continue

    if prev_resi != -1:
        end_strand()

    if len(sys) > 0:
        end_system()

    return configs, systems

# This is what gets picked up by the cli documentation builder
def cli_parser(prog="program_name"):
    parser = argparse.ArgumentParser(prog=prog, description="Convert a PDB or mmCIF file to oxDNA")
    parser.add_argument('pdb_file', type=str, help='The PDB or mmCIF file to convert (format is auto-detected)')
    parser.add_argument('-o', '--output', metavar='output_file', help='The filename to save the output oxDNA files to')
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

    # -o names the output file; strip known structure-file extensions from the basename
    if args.output:
        outbase = args.output
    else:
        outbase = path.splitext(pdb_file)[0]
        log(f"No outfile name provided, defaulting to \"{outbase}\"")

    # Get the string from the PDB
    with open(pdb_file, 'r') as f:
        pdb_str = f.read()

    # Get the oxDNA-style data structures
    configs, systems = PDB_oxDNA(pdb_str)

    for i, (conf, sys) in enumerate(zip(configs, systems)):
        sysn = '' if len(configs) == 1 else f'_{i}'
        outtop = outbase + sysn + '.top'
        outdat = outbase + sysn + '.dat'
        write_top(outtop, sys)
        write_conf(outdat, conf)
        log(f"Wrote outfiles {outtop}, {outdat}")

    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    main()