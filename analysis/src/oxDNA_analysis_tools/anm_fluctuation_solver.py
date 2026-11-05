"""
Build anisotropic network model (ANM / ANMT) parameters for oxDNA protein simulations.

A protein is coarse-grained to one bead per residue (the C-alpha atom). Every pair of
beads closer than a cutoff is joined by a harmonic spring with a single spring constant,
gamma, chosen so that the B-factors predicted by normal mode analysis match the
experimental B-factors from the structure file.  The resulting network is written in the
format read by the AC (ANM), DNANM, DNACT, RNANM and RNACT interactions.

For ANMT, the bending and torsion potential between consecutive residues is included in the
normal mode analysis (the frame rotations are integrated out with a Schur complement) and
the spring constant is fitted with the bending/torsion constants fixed.

Based on the ANM/ANMT parts of models.py from the ANM-oxDNA ANMUtils package
(Procyk et al., Coarse-grained nucleic acid-protein model for hybrid nanotechnology,
Soft Matter 2021), with these differences:
  - the Hessian is assembled with numpy and its pseudo-inverse diagonal is obtained from a
    Cholesky factorisation after removing the rigid-body modes analytically (no SVD);
  - models.py fitted ANMT with the spring-only fluctuations (its ANMT Hessian was never finished);
  - a1 points to the side-chain center of mass (models.py meant this, but an always-true test
    made it the all-atom center of mass) and a3 is normalised, so the cosines in the parameter
    file match what the engine computes from the configuration (models.py's did not);
  - kT = T / 3000 oxDNA energy units, as in the engine.
"""
import time
start_time = time.time()
import argparse
from os import path
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union
import numpy as np
from oxDNA_analysis_tools.UTILS.logger import log, logger_settings
from oxDNA_analysis_tools.UTILS.pdb import NAME_TO_AA, FROM_ANGSTROM_TO_OXDNA, FROM_OXDNA_TO_ANGSTROM
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, System, Strand, Monomer
from oxDNA_analysis_tools.UTILS.RyeReader import write_conf, write_top

#: Converts a mean square fluctuation (A^2) to a B-factor (A^2)
BFACTOR_CONVERSION = 8. * np.pi**2 / 3.
#: One oxDNA spring-constant unit in N/m (one energy unit / one length unit^2)
SIM_FORCE_CONSTANT = 0.05709
#: One oxDNA energy unit in (N/m) * A^2
SIM_ENERGY = SIM_FORCE_CONSTANT * FROM_OXDNA_TO_ANGSTROM**2
#: Boltzmann constant in (N/m) * A^2 / K, so that kB*T/gamma is in A^2 when gamma is in N/m.
#: Uses the oxDNA definition kT = T / 3000 energy units (models.py used 0.00138064852, 1e-4 smaller)
BOLTZMANN = SIM_ENERGY / 3000.

# Residue names that are not in NAME_TO_AA but are common variants of a standard amino acid
_AA_ALIASES = {'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'HSD': 'H', 'HSE': 'H', 'HSP': 'H',
               'CYX': 'C', 'CYM': 'C', 'ASH': 'D', 'GLH': 'E', 'LYN': 'K', 'MSE': 'M'}


# Backbone atoms excluded from the side-chain center of mass that defines a1 (CA is kept, as models.py intended)
_BACKBONE_ATOMS = {'N', 'C', 'O', 'OXT', 'H', 'H1', 'H2', 'H3', 'HN'}


def _one_letter(resname:str) -> Union[str, None]:
    resname = resname.strip().upper()
    return NAME_TO_AA.get(resname, _AA_ALIASES.get(resname))


####################################################################################
##########                        STRUCTURE READING                       ##########
####################################################################################

@dataclass
class ProteinStructure:
    """
        The coarse-grained (one bead per residue) representation of a protein.

        Parameters:
            chain_ids (List[str]) : The chain identifier of each chain, in file order
            sequences (List[str]) : One-letter sequence of each chain (N to C terminus)
            residue_ids (List[str]) : The residue number (plus insertion code) of each residue
            coords (np.ndarray) : (N, 3) C-alpha positions in Angstrom
            bfactors (np.ndarray) : (N,) C-alpha B-factors from the structure file
            com_vectors (np.ndarray) : (N, 3) vector from each C-alpha to the center of mass of the residue without its backbone N, C and O atoms (Angstrom), used to orient the beads for ANMT. Zero for glycine.
            n_vectors (np.ndarray) : (N, 3) vector from each C-alpha to the backbone N of the residue (Angstrom), NaN if the atom is missing
            c_vectors (np.ndarray) : (N, 3) vector from each C-alpha to the backbone C of the residue (Angstrom), NaN if the atom is missing
    """
    chain_ids : List[str]
    sequences : List[str]
    residue_ids : List[str]
    coords : np.ndarray
    bfactors : np.ndarray
    com_vectors : np.ndarray
    n_vectors : np.ndarray
    c_vectors : np.ndarray

    @property
    def n_residues(self) -> int:
        return len(self.coords)

    @property
    def chain_starts(self) -> np.ndarray:
        """Index of the first residue of each chain."""
        return np.cumsum([0] + [len(s) for s in self.sequences[:-1]])


class _ResidueBuilder:
    def __init__(self, chain, key, resname):
        self.chain = chain
        self.group = chain  # residues are split into chains whenever this changes
        self.key = key
        self.resname = resname
        self.atoms = {}  # atom name -> (occupancy, position, bfactor)

    def add(self, name, occupancy, pos, bfactor):
        # Keep the alternate location with the highest occupancy (first one on ties)
        if name not in self.atoms or occupancy > self.atoms[name][0]:
            self.atoms[name] = (occupancy, pos, bfactor)


def _select_models(residues:dict, order:list) -> List[_ResidueBuilder]:
    """
    Combine the models of a structure file, dropping repeated conformations of an ensemble.

    Some programs (e.g. Modeller/Chimera) write one chain per MODEL, while NMR-style ensembles repeat
    the same residues in every model. A model is dropped if its residues are identical to those of an
    earlier kept model and its C-alpha atoms overlap with that model (mean distance below 10 A; identical
    copies of a homo-oligomer are further apart). Chains of later models are labelled <model>:<chain> if
    their chain id is reused.
    """
    def ca_coords(keys):
        return np.array([residues[k].atoms['CA'][1] if 'CA' in residues[k].atoms else np.full(3, np.nan) for k in keys])

    model_residues = {}
    for key in order:
        model_residues.setdefault(key[0], []).append(key)
    kept, signatures, used_chains = [], [], set()
    for model, keys in model_residues.items():
        signature = [(k[1], k[2], residues[k].resname) for k in keys]
        ca = ca_coords(keys)
        if any(signature == sig and np.nanmean(np.linalg.norm(ca - ca_prev, axis=1)) < 10. for sig, ca_prev in signatures):
            log(f"Model {model} is another conformation of an earlier model and is ignored (only the first model of an ensemble is used)", level='warning')
            continue
        signatures.append((signature, ca))
        chains = {k[1] for k in keys}
        for k in keys:
            r = residues[k]
            r.group = (model, k[1])
            if chains & used_chains:
                r.chain = f'{model}:{k[1].strip()}' if k[1].strip() else model
            kept.append(r)
        used_chains |= chains
    return kept


def _residues_to_structure(residues:List[_ResidueBuilder], source:str) -> ProteinStructure:
    chain_ids, sequences, residue_ids, coords, bfactors, com_vectors, n_vectors, c_vectors = [], [], [], [], [], [], [], []
    group = None
    for r in residues:
        aa = _one_letter(r.resname)
        if aa is None:
            log(f"Residue {r.resname} {r.key} in chain {r.chain} is not a standard amino acid, skipping", level='warning')
            continue
        if 'CA' not in r.atoms:
            log(f"Residue {r.resname} {r.key} in chain {r.chain} has no CA atom, skipping", level='warning')
            continue
        if group != r.group:
            group = r.group
            chain_ids.append(r.chain)
            sequences.append('')
        sequences[-1] += aa
        residue_ids.append(r.key)
        _, ca, ca_b = r.atoms['CA']
        coords.append(ca)
        bfactors.append(ca_b)
        com_vectors.append(np.mean([a[1] for name, a in r.atoms.items() if name not in _BACKBONE_ATOMS], axis=0) - ca)
        n_vectors.append(r.atoms['N'][1] - ca if 'N' in r.atoms else np.full(3, np.nan))
        c_vectors.append(r.atoms['C'][1] - ca if 'C' in r.atoms else np.full(3, np.nan))

    if len(coords) == 0:
        raise RuntimeError(f"No amino acids with a CA atom were found in {source}")

    return ProteinStructure(chain_ids, sequences, residue_ids, np.array(coords, dtype=float),
                            np.array(bfactors, dtype=float), np.array(com_vectors, dtype=float),
                            np.array(n_vectors, dtype=float), np.array(c_vectors, dtype=float))


def _read_pdb(text:str, source:str) -> ProteinStructure:
    residues = {}
    order = []
    model = '1'
    for line in text.splitlines():
        if line.startswith('MODEL'):
            model = line[6:].strip() or str(int(model) + 1)
            continue
        if not line.startswith('ATOM'):
            continue
        name = line[12:16].strip()
        chain = line[21]
        key = (model, chain, line[22:26].strip() + line[26].strip())
        occupancy = float(line[54:60]) if line[54:60].strip() else 1.
        bfactor = float(line[60:66]) if line[60:66].strip() else 0.
        pos = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        if key not in residues:
            residues[key] = _ResidueBuilder(chain, key[2], line[17:20].strip())
            order.append(key)
        residues[key].add(name, occupancy, pos, bfactor)
    return _residues_to_structure(_select_models(residues, order), source)


def _read_mmcif(text:str, source:str) -> ProteinStructure:
    from oxDNA_analysis_tools.UTILS.mmcif import parse_atom_site

    def get(row, *keys, default=''):
        for k in keys:
            v = row.get(k, '')
            if v not in ('', '.', '?'):
                return v
        return default

    residues = {}
    order = []
    for row in parse_atom_site(text):
        if row.get('group_PDB', 'ATOM') != 'ATOM':
            continue
        model = get(row, 'pdbx_PDB_model_num', default='1')
        name = get(row, 'auth_atom_id', 'label_atom_id')
        chain = get(row, 'auth_asym_id', 'label_asym_id')
        key = (model, chain, get(row, 'auth_seq_id', 'label_seq_id') + get(row, 'pdbx_PDB_ins_code'))
        pos = np.array([float(row['Cartn_x']), float(row['Cartn_y']), float(row['Cartn_z'])])
        if key not in residues:
            residues[key] = _ResidueBuilder(chain, key[2], get(row, 'auth_comp_id', 'label_comp_id'))
            order.append(key)
        residues[key].add(name, float(get(row, 'occupancy', default='1')), pos, float(get(row, 'B_iso_or_equiv', default='0')))
    return _residues_to_structure(_select_models(residues, order), source)


def read_structure(structure_file:str) -> ProteinStructure:
    """
        Read the C-alpha coarse-grained representation of the protein in a PDB or mmCIF file.

        Only ATOM records of the first model are used.  Residues that are not standard amino
        acids or have no CA atom are skipped. Alternate locations are resolved by occupancy.

        Parameters:
            structure_file (str) : Path to a PDB or mmCIF file (format is auto-detected)

        Returns:
            ProteinStructure : The residues, coordinates and B-factors of the protein
    """
    from oxDNA_analysis_tools.UTILS.mmcif import is_mmcif
    with open(structure_file) as f:
        text = f.read()
    if is_mmcif(text):
        return _read_mmcif(text, structure_file)
    return _read_pdb(text, structure_file)


def read_bfactor_file(bfactor_file:str) -> np.ndarray:
    """
        Read B-factors from a text file containing one number per residue (whitespace separated, '#' starts a comment).

        Parameters:
            bfactor_file (str) : Path to the file

        Returns:
            np.ndarray : The B-factors in A^2
    """
    values = []
    with open(bfactor_file) as f:
        for line in f:
            values.extend(float(v) for v in line.split('#')[0].split())
    return np.array(values)


####################################################################################
##########                     NETWORK AND NORMAL MODES                   ##########
####################################################################################

def contact_pairs(coords:np.ndarray, cutoff:float) -> Tuple[np.ndarray, np.ndarray]:
    """
        Find every pair of beads closer than the cutoff.

        Parameters:
            coords (np.ndarray) : (N, 3) positions
            cutoff (float) : The cutoff distance, in the same units as coords

        Returns:
            (np.ndarray, np.ndarray) : The i < j indices of each pair, sorted by i then j
    """
    from scipy.spatial import cKDTree
    pairs = cKDTree(coords).query_pairs(cutoff, output_type='ndarray')
    if len(pairs) == 0:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    pairs = np.sort(pairs, axis=1)
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    return pairs[order, 0], pairs[order, 1]


def build_hessian(coords:np.ndarray, i:np.ndarray, j:np.ndarray, spring_constants:Union[float, np.ndarray]=1.) -> np.ndarray:
    """
        Assemble the (3N, 3N) ANM Hessian of a spring network.

        Parameters:
            coords (np.ndarray) : (N, 3) positions
            i, j (np.ndarray) : The indices of the beads joined by each spring
            spring_constants (float or np.ndarray) : The spring constant of all springs or of each spring

        Returns:
            np.ndarray : The Hessian
    """
    n = len(coords)
    d = coords[j] - coords[i]
    k = np.broadcast_to(np.asarray(spring_constants, dtype=float), i.shape)
    blocks = (k / np.einsum('pi,pi->p', d, d))[:, None, None] * d[:, :, None] * d[:, None, :]

    hess = np.zeros((n, 3, n, 3))
    hess[i, :, j, :] = -blocks
    hess[j, :, i, :] = -blocks
    diag = np.zeros((n, 3, 3))
    np.add.at(diag, i, blocks)
    np.add.at(diag, j, blocks)
    idx = np.arange(n)
    hess[idx, :, idx, :] = diag
    return hess.reshape(3 * n, 3 * n)


def angular_hessian(coords:np.ndarray, p:np.ndarray, a1s:np.ndarray, a3s:np.ndarray, kb:float, kt:float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
        Hessian of the ANMT bending/torsion potential at its minimum (DNANMInteraction::_protein_ang_pot).

        For each pair (p, q = p + 1) with r the unit vector from p to q, the energy is
        kb/2 [(r.a1p - a0)^2 + (-r.a1q - b0)^2] + kt/2 [(a1p.a1q - c0)^2 + (a3p.a3q - d0)^2].
        The degrees of freedom are the positions and an infinitesimal rotation theta of each particle's frame (da = theta x a).
        At the minimum each term contributes k g g^T, with g the gradient of the cosine.

        Parameters:
            coords (np.ndarray) : (N, 3) positions
            p (np.ndarray) : Index of the first particle of each angular pair (the second is p + 1)
            a1s, a3s (np.ndarray) : (N, 3) particle frames, as the engine sees them (orthonormal)
            kb, kt (float) : Bending and torsion constants, in energy units of (spring constant units) * (coords units)^2

        Returns:
            (np.ndarray, np.ndarray, np.ndarray) : The (3N, 3N) position-position, position-rotation and rotation-rotation blocks
    """
    n = len(coords)
    q = p + 1
    r = coords[q] - coords[p]
    length = np.linalg.norm(r, axis=1)
    u = r / length[:, None]
    a, b, A, B = a1s[p], a1s[q], a3s[p], a3s[q]

    def proj(v):  # (I - u u^T) v / |r|, the derivative of u.v with respect to the position of q
        return (v - u * np.einsum('ij,ij->i', u, v)[:, None]) / length[:, None]

    zero = np.zeros_like(u)
    # gradient of each cosine with respect to (x_p, x_q, theta_p, theta_q)
    terms = [(kb, -proj(a), proj(a), np.cross(a, u), zero),
             (kb, proj(b), -proj(b), zero, np.cross(u, b)),
             (kt, zero, zero, np.cross(a, b), np.cross(b, a)),
             (kt, zero, zero, np.cross(A, B), np.cross(B, A))]

    hxx = np.zeros((n, 3, n, 3))
    hxt = np.zeros((n, 3, n, 3))
    htt = np.zeros((n, 3, n, 3))
    for k, gxp, gxq, gtp, gtq in terms:
        gx = (gxp, gxq)
        gt = (gtp, gtq)
        idx = (p, q)
        for s in range(2):
            for t in range(2):
                # pairs are unique, but particle p appears in two pairs, so accumulate with add.at
                np.add.at(hxx, (idx[s], slice(None), idx[t]), k * np.einsum('pi,pj->pij', gx[s], gx[t]))
                np.add.at(hxt, (idx[s], slice(None), idx[t]), k * np.einsum('pi,pj->pij', gx[s], gt[t]))
                np.add.at(htt, (idx[s], slice(None), idx[t]), k * np.einsum('pi,pj->pij', gt[s], gt[t]))
    return hxx.reshape(3 * n, 3 * n), hxt.reshape(3 * n, 3 * n), htt.reshape(3 * n, 3 * n)


def eliminate_rotations(hxx:np.ndarray, hxt:np.ndarray, htt:np.ndarray) -> np.ndarray:
    """
        Effective position Hessian after integrating out the rotational degrees of freedom (Schur complement).

        H_eff = H_xx - H_xt H_tt^+ H_tx.  H_tt is banded (rotations only couple through consecutive particles),
        so a banded Cholesky factorisation is used.

        H_tt is singular: rotating every frame of a chain about its own a1 only changes the N-1 torsion (a3.a3)
        terms, which leaves one zero mode per chain.  These modes do not couple to the positions (the energy is
        invariant), so they are regularised with a diagonal shift.  With frames that follow the backbone, consecutive
        a3 vectors are nearly parallel and H_tt also has genuine soft modes (~1e-10 of its largest eigenvalue), which a
        large shift would stiffen, so the smallest shift that allows the factorisation is used (1e-14 * max(diag H_tt),
        then 1e-12 and 1e-10).  If the factorisation still fails, a thresholded pseudo-inverse is used.

        Parameters:
            hxx, hxt, htt (np.ndarray) : (3N, 3N) blocks of the full Hessian

        Returns:
            np.ndarray : (3N, 3N) effective Hessian of the positions
    """
    import scipy.linalg as sla
    active = np.where(np.diag(htt) > 0)[0]  # rotations of particles without angular terms are free and decoupled
    if len(active) == 0:
        return hxx.copy()
    htt = htt[np.ix_(active, active)]
    htx = hxt[:, active].T
    rows, cols = np.nonzero(htt)
    bw = int(np.max(np.abs(cols - rows)))
    try:
        if bw > 64:
            raise np.linalg.LinAlgError("not banded")
        banded = np.zeros((bw + 1, len(active)))
        for k in range(bw + 1):
            banded[bw - k, k:] = np.diagonal(htt, k)
        diag_max = banded[bw].max()
        chol = None
        for shift in (1e-14, 1e-12, 1e-10):
            shifted = banded.copy()
            shifted[bw] += shift * diag_max
            try:
                chol = sla.cholesky_banded(shifted, lower=False, check_finite=False)
                break
            except np.linalg.LinAlgError:
                continue
        if chol is None:
            raise np.linalg.LinAlgError("the rotational block is not positive definite")
        # H_tt = U^T U, so H_xt H_tt^-1 H_tx = Z^T Z with Z = U^-T H_tx
        z = _banded_upper_solve_transpose(chol, bw, htx)
    except np.linalg.LinAlgError:
        log("Using a pseudo-inverse for the rotational degrees of freedom", level='warning')
        w, v = np.linalg.eigh(htt)
        keep = w > 1e-15 * w[-1]
        z = (v[:, keep].T @ htx) / np.sqrt(w[keep])[:, None]
    return hxx - z.T @ z


def _banded_upper_solve_transpose(chol:np.ndarray, bw:int, rhs:np.ndarray) -> np.ndarray:
    """Solve U^T Z = rhs for Z, with U upper triangular in LAPACK banded storage."""
    from scipy.linalg import solve_banded
    n = chol.shape[1]
    # U^T is lower triangular with bandwidth bw; convert the storage of U to that of U^T
    lower = np.zeros((bw + 1, n))
    for k in range(bw + 1):
        lower[k, :n - k] = chol[bw - k, k:]
    return solve_banded((bw, 0), lower, rhs, check_finite=False)


def _rigid_body_modes(coords:np.ndarray, components:Optional[np.ndarray]=None) -> np.ndarray:
    """Orthonormal basis (3N, m) of the rigid translations and rotations of each connected component."""
    n = len(coords)
    if components is None:
        components = np.zeros(n, dtype=int)
    basis = []
    for c in np.unique(components):
        members = np.where(components == c)[0]
        x = coords[members] - coords[members].mean(axis=0)
        modes = np.zeros((n, 3, 6))
        for a in range(3):
            modes[members, a, a] = 1.
            e = np.zeros(3)
            e[a] = 1.
            modes[members, :, 3 + a] = np.cross(e, x)
        u, s, _ = np.linalg.svd(modes.reshape(3 * n, 6), full_matrices=False)
        basis.append(u[:, s > 1e-8 * s[0]])  # two-particle components have 5 modes, single particles 3
    return np.hstack(basis)


def mean_square_fluctuations(hessian:np.ndarray, coords:np.ndarray, temperature:float=300., method:str='auto',
                             components:Optional[np.ndarray]=None) -> np.ndarray:
    """
        Per-bead mean square fluctuations, kB*T * trace of the diagonal blocks of the Hessian pseudo-inverse.

        Parameters:
            hessian (np.ndarray) : (3N, 3N) Hessian in N/m
            coords (np.ndarray) : (N, 3) positions the Hessian was built from (used to remove the rigid-body modes)
            temperature (float) : Temperature in K
            method (str) : 'cholesky' removes the rigid-body modes analytically and inverts with a Cholesky factorisation; 'eigh' uses the full eigendecomposition; 'auto' (default) tries 'cholesky' and falls back to 'eigh' if the network has extra zero modes.
            components (np.ndarray) : Connected component of each bead in the network (default: all connected)

        Returns:
            np.ndarray : (N,) mean square fluctuations in A^2
    """
    import scipy.linalg as sla
    n = len(coords)
    if method not in ('auto', 'cholesky', 'eigh'):
        raise RuntimeError(f"Unknown method {method}, choose from 'auto', 'cholesky' or 'eigh'")
    q = _rigid_body_modes(coords, components)
    if method in ('auto', 'cholesky'):
        # H Q = 0 for the rigid-body modes Q, so (H + c Q Q^T)^-1 = H^+ + Q Q^T / c.
        # If the network has no other zero modes, H + c Q Q^T is positive definite.
        c = np.trace(hessian) / (3 * n)
        a = hessian + c * (q @ q.T)
        try:
            chol, info = sla.lapack.dpotrf(a, lower=1, overwrite_a=1, clean=1)
            if info != 0:
                raise np.linalg.LinAlgError("matrix is not positive definite")
            linv, info = sla.lapack.dtrtri(chol, lower=1, overwrite_c=1)
            if info != 0:
                raise np.linalg.LinAlgError("triangular inversion failed")
            # diag(A^-1) = squared column norms of L^-1
            diag = np.einsum('ij,ij->j', linv, linv) - np.einsum('ij,ij->i', q, q) / c
            return BOLTZMANN * temperature * diag.reshape(n, 3).sum(axis=1)
        except np.linalg.LinAlgError:
            if method == 'cholesky':
                raise
            log("The network has zero-energy modes besides the rigid-body motions (some residues are under-constrained), using the full eigendecomposition", level='warning')

    w, v = sla.eigh(hessian, driver='evd')
    nonzero = w > 1e-6 * w[-1]
    n_zero = np.sum(~nonzero)
    if n_zero != q.shape[1]:
        log(f"The Hessian has {n_zero} zero-energy modes instead of {q.shape[1]}, the fluctuations of under-constrained residues are not meaningful; consider increasing the cutoff", level='warning')
    diag = (v[:, nonzero]**2) @ (1 / w[nonzero])
    return BOLTZMANN * temperature * diag.reshape(n, 3).sum(axis=1)


def fit_spring_constant(msf:np.ndarray, exp_bfactors:np.ndarray, method:str='linear', max_deviations:float=2.,
                        bmax:float=10000., grid:Tuple[float, float, float]=(0.001, 5., 0.001)) -> Tuple[float, np.ndarray]:
    """
        Find the spring constant gamma for which the predicted B-factors best match the experimental ones.

        The predicted B-factors are BFACTOR_CONVERSION * msf / gamma, where msf are the fluctuations at gamma = 1 N/m.

        Parameters:
            msf (np.ndarray) : Mean square fluctuations with a spring constant of 1 N/m (A^2)
            exp_bfactors (np.ndarray) : Experimental B-factors (A^2)
            method (str) : 'linear' (default): least squares fit through the origin, ignoring residues whose predicted B-factor is more than max_deviations standard deviations from the mean or above bmax.  'grid': grid search minimising sum |B_pred^2 - B_exp^2| over all residues (the method used for the original ANM-oxDNA example files).
            max_deviations (float) : Outlier threshold for the linear fit, in standard deviations
            bmax (float) : Residues with a fitted B-factor above this value are excluded from the linear fit
            grid (tuple) : (start, stop, step) of the gamma grid in N/m for the grid search

        Returns:
            (float, np.ndarray) : gamma in N/m and a boolean mask of the residues used in the fit
    """
    flex = BFACTOR_CONVERSION * np.asarray(msf, dtype=float)
    exp_bfactors = np.asarray(exp_bfactors, dtype=float)

    if method == 'grid':
        start, stop, step = grid
        gammas = np.arange(start, stop, step)
        residual = np.abs((flex[None, :] / gammas[:, None])**2 - exp_bfactors[None, :]**2).sum(axis=1)
        return float(np.argmin(residual) * step + start), np.ones(len(flex), dtype=bool)

    if method != 'linear':
        raise RuntimeError(f"Unknown fitting method {method}, choose from 'linear' or 'grid'")

    keep = np.abs(flex - flex.mean()) < max_deviations * flex.std()
    slope = 1.
    used = keep.copy()
    for _ in range(len(flex) + 1):
        used = keep & (flex * slope < bmax)
        if not np.any(used):
            raise RuntimeError("No residues left to fit the spring constant to, check bmax and max_deviations")
        slope = np.dot(flex[used], exp_bfactors[used]) / np.dot(flex[used], flex[used])
        if np.all(flex[used] * slope < bmax):
            break
    if slope <= 0:
        raise RuntimeError("The fitted spring constant is negative, the experimental B-factors are not correlated with the network fluctuations")
    return float(1. / slope), used


def _perpendicular_unit(v:np.ndarray, ref:np.ndarray) -> np.ndarray:
    """Unit vector along the component of v perpendicular to the unit vector ref, or None if it is (nearly) parallel to it."""
    if v is None or not np.all(np.isfinite(v)):
        return None
    p = v - np.dot(v, ref) * ref
    norm = np.linalg.norm(p)
    if norm < 1e-3 * max(np.linalg.norm(v), 1e-12):
        return None
    return p / norm


def _any_perpendicular(v:np.ndarray) -> np.ndarray:
    """A unit vector perpendicular to v (last resort: it depends on the orientation of the input coordinates)."""
    axis = np.zeros(3)
    axis[np.argmin(np.abs(v))] = 1.
    p = np.cross(v, axis)
    return p / np.linalg.norm(p)


def orientation_vectors(structure:ProteinStructure) -> Tuple[np.ndarray, np.ndarray]:
    """
        Body-frame vectors of each bead, used by the ANMT bending and torsion potentials. Both are defined by the residue
        itself, so the parameters do not depend on how the structure is oriented in the input file.

        a1 points from the C-alpha to the center of mass of the side chain (all atoms except the backbone N, C and O, CA
        included). For glycine (and any residue without side-chain atoms) it points to the ideal position of the C-beta,
        computed from the backbone N, CA and C atoms, or, if they are missing, along the normal to the plane of the C-alpha
        and its two neighbours in the chain.

        a3 points along the backbone, from the N to the C atom of the residue (or from the previous to the next C-alpha if
        they are missing), made perpendicular to a1. As for nucleotides, it follows the direction of the chain.

        The frames are orthonormal, so the engine (which orthonormalises them when it reads the configuration) sees the same vectors.

        Parameters:
            structure (ProteinStructure) : The protein

        Returns:
            (np.ndarray, np.ndarray) : (N, 3) a1 and a3 vectors
    """
    n = structure.n_residues
    x = structure.coords
    chain = np.repeat(np.arange(len(structure.sequences)), [len(s) for s in structure.sequences])
    a1s = np.zeros((n, 3))
    a3s = np.zeros((n, 3))
    arbitrary = 0

    def neighbours(r):
        prev_r = r - 1 if r > 0 and chain[r - 1] == chain[r] else None
        next_r = r + 1 if r < n - 1 and chain[r + 1] == chain[r] else None
        return prev_r, next_r

    for r in range(n):
        prev_r, next_r = neighbours(r)
        nv, cv = structure.n_vectors[r], structure.c_vectors[r]
        has_nc = np.all(np.isfinite(nv)) and np.all(np.isfinite(cv))

        # a1: towards the side chain
        a1 = structure.com_vectors[r].copy()
        if np.linalg.norm(a1) < 1e-8:
            if has_nc:
                # ideal C-beta position (the coefficients of trRosetta/AlphaFold)
                b = -nv  # CA - N
                c = cv  # C - CA
                a1 = -0.58273431 * np.cross(b, c) + 0.56802827 * b - 0.54067466 * c
            elif prev_r is not None and next_r is not None:
                a1 = np.cross(x[next_r] - x[r], x[prev_r] - x[r])
        if np.linalg.norm(a1) < 1e-8:
            bond = x[next_r] - x[r] if next_r is not None else (x[r] - x[prev_r] if prev_r is not None else np.array([0., 0., 1.]))
            a1 = _any_perpendicular(bond)
            arbitrary += 1
        a1 = a1 / np.linalg.norm(a1)

        # a3: along the backbone, perpendicular to a1
        candidates = []
        if has_nc:
            candidates.append(cv - nv)
        if prev_r is not None and next_r is not None:
            candidates.append(x[next_r] - x[prev_r])
        if next_r is not None:
            candidates.append(x[next_r] - x[r])
        if prev_r is not None:
            candidates.append(x[r] - x[prev_r])
        a3 = None
        for v in candidates:
            a3 = _perpendicular_unit(v, a1)
            if a3 is not None:
                break
        if a3 is None:
            a3 = _any_perpendicular(a1)
            arbitrary += 1

        a1s[r] = a1
        a3s[r] = a3

    if arbitrary > 0:
        log(f"{arbitrary} bead orientation vectors could not be defined from the structure (missing backbone atoms and chain neighbours) and were chosen arbitrarily", level='warning')
    return a1s, a3s


@dataclass
class ANMNetwork:
    """
        A fitted elastic network, ready to be written as an oxDNA parameter file.

        Parameters:
            model (str) : 'anm' or 'anmt'
            cutoff (float) : Spring cutoff in Angstrom
            temperature (float) : Temperature used for the fit (K)
            i, j (np.ndarray) : 0-based indices of the residues joined by each spring (i < j)
            r0 (np.ndarray) : Equilibrium length of each spring in oxDNA length units
            k (np.ndarray) : Spring constant of each spring in oxDNA units
            gamma (float) : The fitted spring constant in N/m (multiply by 100 for pN/A)
            exp_bfactors (np.ndarray) : Experimental B-factors (A^2)
            pred_bfactors (np.ndarray) : B-factors predicted by the network (A^2)
            fitted (np.ndarray) : Boolean mask of the residues used in the fit
            a1s, a3s (np.ndarray) : Bead orientations, written to the configuration
            angles (np.ndarray or None) : For ANMT, the (a0, b0, c0, d0) cosines of each spring with j == i + 1, else None
            kb, kt (float or None) : ANMT bending and torsion constants (oxDNA units) assumed by the fit, None for ANM
            write_kbkt (bool) : Write kb and kt on every backbone line of the parameter file (otherwise set bending_k and torsion_k in the input file)
    """
    model : str
    cutoff : float
    temperature : float
    i : np.ndarray
    j : np.ndarray
    r0 : np.ndarray
    k : np.ndarray
    gamma : float
    exp_bfactors : np.ndarray
    pred_bfactors : np.ndarray
    fitted : np.ndarray
    a1s : np.ndarray
    a3s : np.ndarray
    angles : Optional[np.ndarray] = None
    kb : Optional[float] = None
    kt : Optional[float] = None
    write_kbkt : bool = False

    @property
    def rmsf(self) -> np.ndarray:
        """Predicted root mean square fluctuation of each residue (A)."""
        return np.sqrt(self.pred_bfactors / BFACTOR_CONVERSION)

    @property
    def correlation(self) -> float:
        """Pearson correlation between predicted and experimental B-factors over the fitted residues."""
        if np.sum(self.fitted) < 2 or np.std(self.exp_bfactors[self.fitted]) == 0:
            return float('nan')
        return float(np.corrcoef(self.pred_bfactors[self.fitted], self.exp_bfactors[self.fitted])[0, 1])


def _engine_frames(a1s:np.ndarray, a3s:np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Orthonormalise the frames the way the engine does when it reads a configuration."""
    a3 = a3s / np.linalg.norm(a3s, axis=1)[:, None]
    a1 = a1s / np.linalg.norm(a1s, axis=1)[:, None]
    a1 = a1 - a3 * np.einsum('ij,ij->i', a1, a3)[:, None]
    return a1 / np.linalg.norm(a1, axis=1)[:, None], a3


def anm_fluctuation_solver(structure:ProteinStructure, model:str='anm', cutoff:float=13., temperature:float=300.,
                           bfactors:Optional[np.ndarray]=None, spring_constant:Optional[float]=None, fit:str='linear',
                           max_deviations:float=2., bmax:float=10000., kb:float=50., kt:float=50., write_kbkt:bool=False,
                           method:str='auto') -> ANMNetwork:
    """
        Build an ANM or ANMT network for a protein and fit its spring constant to experimental B-factors.

        For ANMT, the bending and torsion terms between consecutive residues are included in the normal mode analysis
        (with the rotations of the particle frames integrated out) and only the spring constant is fitted, with kb and kt fixed.

        Parameters:
            structure (ProteinStructure) : The protein, from read_structure
            model (str) : 'anm' (springs only) or 'anmt' (springs plus bending/torsion between consecutive residues)
            cutoff (float) : Residues closer than this (Angstrom) are joined by a spring (default 13)
            temperature (float) : Temperature of the experiment in K (default 300)
            bfactors (np.ndarray) : Experimental B-factors, one per residue.  Defaults to the C-alpha B-factors in the structure.
            spring_constant (float) : Use this spring constant (oxDNA units) instead of fitting it to the B-factors
            fit (str) : How to fit the spring constant, 'linear' (default) or 'grid', see fit_spring_constant.  'grid' always uses the spring-only fluctuations (as the original ANM-oxDNA scripts did, also for ANMT).
            max_deviations (float) : Outlier threshold of the linear fit, in standard deviations (default 2)
            bmax (float) : Residues predicted above this B-factor are excluded from the linear fit (default 10000)
            kb, kt (float) : ANMT only. Bending and torsion constants in oxDNA units (default 50, 50). They must match the values used in the simulation.
            write_kbkt (bool) : ANMT only. Write kb and kt on every backbone line of the parameter file; otherwise set bending_k and torsion_k in the input file (default False)
            method (str) : Linear algebra used for the fluctuations, see mean_square_fluctuations

        Returns:
            ANMNetwork : The fitted network
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    model = model.lower()
    if model not in ('anm', 'anmt'):
        raise RuntimeError(f"Unknown model {model}, choose from 'anm' or 'anmt'")
    if kb < 0 or kt < 0:
        raise RuntimeError("kb and kt must be non-negative")
    if write_kbkt and model != 'anmt':
        raise RuntimeError("kb and kt can only be written for the anmt model")
    angular = model == 'anmt'

    n = structure.n_residues
    exp_b = structure.bfactors if bfactors is None else np.asarray(bfactors, dtype=float)
    if len(exp_b) != n:
        raise RuntimeError(f"Got {len(exp_b)} B-factors for {n} residues")

    xyz = structure.coords
    i, j = contact_pairs(xyz, cutoff)
    n_comp, components = connected_components(coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)), directed=False)
    log(f"{n} residues in {len(structure.sequences)} chains, {len(i)} springs with a {cutoff} A cutoff")
    if n_comp > 1:
        log(f"The network has {n_comp} disconnected parts, consider increasing the cutoff", level='warning')

    a1s, a3s = orientation_vectors(structure)
    spring_hessian = build_hessian(xyz, i, j)  # gamma = 1 N/m
    msf_springs = mean_square_fluctuations(spring_hessian, xyz, temperature, method, components)

    if angular:
        # the engine applies the angular potential to every spring with j == i + 1
        backbone = i[j == i + 1]
        a1e, a3e = _engine_frames(a1s, a3s)
        hxx, hxt, htt = angular_hessian(xyz, backbone, a1e, a3e, kb * SIM_ENERGY, kt * SIM_ENERGY)
        # the Schur complement only changes H_xx, so the angular part of H_eff does not depend on gamma
        angular_eff = eliminate_rotations(hxx, hxt, htt)
        del hxx, hxt, htt

        cache = {}

        def msf_at(gamma):
            if gamma not in cache:
                cache[gamma] = mean_square_fluctuations(spring_hessian * gamma + angular_eff, xyz, temperature, method, components)
            return cache[gamma]
    else:
        def msf_at(gamma):
            return msf_springs / gamma

    if spring_constant is not None:
        gamma = spring_constant * SIM_FORCE_CONSTANT
        fitted = np.ones(n, dtype=bool)
    else:
        if np.allclose(exp_b, exp_b[0]):
            raise RuntimeError("The B-factors are all the same, so the spring constant cannot be fitted; give one with spring_constant")
        gamma, fitted = fit_spring_constant(msf_springs, exp_b, fit, max_deviations, bmax)
        if not np.all(fitted):
            log(f"Residues ignored by the fit (0-based): {np.where(~fitted)[0].tolist()}")
        if angular and fit == 'grid':
            log("The grid fit uses the spring-only fluctuations, bending and torsion are not included in the fit")
        elif angular and (kb > 0 or kt > 0):
            gamma = _fit_angular(msf_at, exp_b, fitted, gamma)
    pred_b = BFACTOR_CONVERSION * msf_at(gamma)

    coords = xyz * FROM_ANGSTROM_TO_OXDNA
    d = coords[j] - coords[i]
    r0 = np.linalg.norm(d, axis=1)
    k = np.full(len(i), gamma / SIM_FORCE_CONSTANT)

    angles = None
    if angular:
        # cosines of the frames written to the configuration, as the engine will compute them
        bb = j == i + 1
        rhat = d[bb] / r0[bb, None]
        ii, jj = i[bb], j[bb]
        angles = np.stack([np.einsum('pi,pi->p', a1s[ii], rhat),
                           -np.einsum('pi,pi->p', a1s[jj], rhat),
                           np.einsum('pi,pi->p', a1s[ii], a1s[jj]),
                           np.einsum('pi,pi->p', a3s[ii], a3s[jj])], axis=1)

    return ANMNetwork(model, cutoff, temperature, i, j, r0, k, gamma, exp_b, pred_b, fitted, a1s, a3s, angles,
                      kb if angular else None, kt if angular else None, write_kbkt)


def _fit_angular(msf_at, exp_b:np.ndarray, fitted:np.ndarray, gamma0:float, tol:float=1e-9) -> float:
    """
    Find the spring constant for which the least-squares slope (through the origin) of the experimental on the
    predicted B-factors is 1, with the bending and torsion constants fixed.  Each evaluation is a normal mode
    analysis, so the root of log(slope) in log(gamma) is found with a safeguarded secant method.
    """
    t = exp_b[fitted]

    def log_slope(x):
        m = BFACTOR_CONVERSION * msf_at(np.exp(x))[fitted]
        return float(np.log(np.dot(t, m) / np.dot(m, m)))

    # log(slope) increases with log(gamma), with derivative 1 for springs only and < 1 with angular terms.
    # The spring-only fit is an upper bound because the angular terms only stiffen the network.
    x0 = np.log(gamma0)
    f0 = log_slope(x0)
    for _ in range(20):
        if f0 >= 0:
            break
        x0 += np.log(2)
        f0 = log_slope(x0)
    if abs(f0) < tol:
        return float(np.exp(x0))
    hi, f_hi = x0, f0
    lo, f_lo = None, None
    x1 = x0 - f0
    for _ in range(60):
        f1 = log_slope(x1)
        if abs(f1) < tol:
            return float(np.exp(x1))
        if f1 > 0:
            hi, f_hi = x1, f1
        else:
            lo, f_lo = x1, f1
        # secant step through the last two points, kept inside the bracket
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0) if f1 != f0 else x1 - f1
        x0, f0 = x1, f1
        if lo is None:
            # no point below the root yet: keep stepping down, at least as far as a unit-slope step
            x2 = min(x2, x1 - f1) if f1 > 0 else x2
        elif not (lo < x2 < hi):
            x2 = 0.5 * (lo + hi)
        if lo is None and x2 < np.log(gamma0) - 50:
            break
        x1 = x2
    if lo is None:
        raise RuntimeError("The bending and torsion terms alone make the protein stiffer than the experimental B-factors; lower kb and kt")
    return float(np.exp(x1))


####################################################################################
##########                          FILE WRITERS                          ##########
####################################################################################

def write_par(path:str, network:ANMNetwork, index_offset:int=0) -> None:
    """
        Write the oxDNA parameter file of the network.

        The first line is the number of amino acids, followed by one line per spring: ``i j r0 s k``.
        For ANMT, springs with j == i + 1 also carry the four equilibrium cosines ``a0 b0 c0 d0`` and, if set, ``kb kt``.

        Parameters:
            path (str) : The file to write
            network (ANMNetwork) : The network
            index_offset (int) : Added to every particle index, for proteins that are not the first molecules in the topology
    """
    angles = {}
    if network.angles is not None:
        bb = np.where(network.j == network.i + 1)[0]
        angles = {p: a for p, a in zip(bb, network.angles)}
    with open(path, 'w') as f:
        f.write(f'{len(network.a1s)}\n')
        for p, (i, j, r0, k) in enumerate(zip(network.i.tolist(), network.j.tolist(), network.r0.tolist(), network.k.tolist())):
            line = f'{i + index_offset} {j + index_offset} {r0} s {k}'
            if p in angles:
                line += ' ' + ' '.join(str(float(a)) for a in angles[p])
                if network.write_kbkt:
                    line += f' {float(network.kb)} {float(network.kt)}'
            f.write(line + '\n')


def protein_system(structure:ProteinStructure) -> System:
    """
        The oxDNA System (new-style topology) of the protein, one peptide strand per chain.

        Parameters:
            structure (ProteinStructure) : The protein

        Returns:
            System : The system
    """
    system = System()
    idx = 0
    for c, seq in enumerate(structure.sequences):
        strand = Strand(-(c + 1))
        for aa in seq:
            strand.append(Monomer(idx, aa, strand))
            idx += 1
        system.append(strand)
    return system


def get_old_top_string(structure:ProteinStructure, index_offset:int=0) -> str:
    """
        Old-style topology of the protein: ``<strand id < 0> <aa> <nside> <cside>`` per residue.

        The header ``N N_strands 0 N 0`` is accepted both by the protein-only (ANM/ANMT) and the hybrid (DNANM/RNANM) interactions.

        Parameters:
            structure (ProteinStructure) : The protein
            index_offset (int) : Added to the nside/cside indices

        Returns:
            str : The topology
    """
    n = structure.n_residues
    lines = [f'{n} {len(structure.sequences)} 0 {n} 0']
    idx = 0
    for c, seq in enumerate(structure.sequences):
        for r, aa in enumerate(seq):
            nside = idx - 1 + index_offset if r > 0 else -1
            cside = idx + 1 + index_offset if r < len(seq) - 1 else -1
            lines.append(f'{-(c + 1)} {aa} {nside} {cside}')
            idx += 1
    return '\n'.join(lines) + '\n'


def protein_configuration(structure:ProteinStructure, a1s:np.ndarray, a3s:np.ndarray) -> Configuration:
    """
        Starting configuration of the protein in oxDNA units, centred on the origin.

        Parameters:
            structure (ProteinStructure) : The protein
            a1s, a3s (np.ndarray) : Bead orientations

        Returns:
            Configuration : The configuration
    """
    pos = structure.coords * FROM_ANGSTROM_TO_OXDNA
    pos = pos - pos.mean(axis=0)
    box = 2.5 * np.ceil(np.max(pos.max(axis=0) - pos.min(axis=0)))
    return Configuration(0, np.array([box, box, box]), np.zeros(3), pos, a1s, a3s)


def plot_bfactors(path:str, network:ANMNetwork, structure:ProteinStructure) -> None:
    """
        Plot the experimental and predicted B-factors of each residue.

        Parameters:
            path (str) : Output image file
            network (ANMNetwork) : The fitted network
            structure (ProteinStructure) : The protein (for chain boundaries)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from oxDNA_analysis_tools.UTILS.constants import FIG_DPI

    x = np.arange(structure.n_residues)
    fig, ax = plt.subplots(figsize=(9, 4))
    for s in structure.chain_starts[1:]:
        ax.axvline(s - 0.5, color='#d0d0cc', linewidth=1, zorder=0)
    for s, cid in zip(structure.chain_starts, structure.chain_ids):
        ax.annotate(cid, (s, 1), xycoords=('data', 'axes fraction'), xytext=(3, -3), textcoords='offset points',
                    va='top', fontsize=9, color='#6b6b66')
    ax.plot(x, network.exp_bfactors, color='#2a78d6', linewidth=2, label='Experimental')
    ax.plot(x, network.pred_bfactors, color='#eb6834', linewidth=2,
            label=f'{network.model.upper()} ($\\gamma$ = {network.gamma * 100:.3g} pN/$\\AA$)')
    if not np.all(network.fitted):
        out = ~network.fitted
        ax.plot(x[out], network.pred_bfactors[out], linestyle='none', marker='o', markersize=4,
                markerfacecolor='none', markeredgecolor='#6b6b66', label='Excluded from fit')
    ax.set_xlabel('Residue')
    ax.set_ylabel(r'B-factor ($\AA^2$)')
    ax.set_xlim(-0.5, structure.n_residues - 0.5)
    ax.set_ylim(bottom=0)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, loc='lower right', bbox_to_anchor=(1, 1), ncol=3)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


def cli_parser(prog="anm_fluctuation_solver.py"):
    parser = argparse.ArgumentParser(prog=prog, description="Build ANM/ANMT network parameters for an oxDNA protein simulation, fitting the spring constant to experimental B-factors. Writes <output>.par, <output>.top and <output>.dat")
    parser.add_argument('structure', type=str, help='PDB or mmCIF file of the protein (only C-alpha atoms of standard amino acids in ATOM records of the first model are used)')
    parser.add_argument('-o', '--output', metavar='output_base', help='Base name of the output files (default: the structure file name without extension)')
    parser.add_argument('-m', '--model', choices=['anm', 'anmt'], default='anm', help='anm: springs only; anmt: springs plus bending/torsion between consecutive residues (default anm)')
    parser.add_argument('-c', '--cutoff', type=float, default=13., help='Spring cutoff distance in Angstrom (default 13)')
    parser.add_argument('-t', '--temperature', type=float, default=300., help='Temperature of the experiment in K (default 300)')
    parser.add_argument('-b', '--bfactors', metavar='bfactor_file', help='Text file with one experimental B-factor per residue (default: the C-alpha B-factors of the structure)')
    parser.add_argument('-k', '--spring-constant', type=float, dest='spring_constant', help='Use this spring constant (oxDNA units) instead of fitting one to the B-factors')
    parser.add_argument('-f', '--fit', choices=['linear', 'grid'], default='linear', help='linear: least-squares fit ignoring outliers (default); grid: grid search over all residues, as used for the original ANM-oxDNA examples')
    parser.add_argument('--max-deviations', type=float, default=2., dest='max_deviations', help='Residues whose predicted B-factor is further than this many standard deviations from the mean are ignored by the linear fit (default 2)')
    parser.add_argument('--kb', type=float, default=50., help='ANMT bending constant in oxDNA units, included in the fit (default 50). Must match bending_k in the input file unless --write-kbkt is used')
    parser.add_argument('--kt', type=float, default=50., help='ANMT torsion constant in oxDNA units, included in the fit (default 50). Must match torsion_k in the input file unless --write-kbkt is used')
    parser.add_argument('--write-kbkt', action='store_true', dest='write_kbkt', help='ANMT only. Write kb and kt on every backbone line of the parameter file (then do not set bending_k/torsion_k in the input file)')
    parser.add_argument('-i', '--index-offset', type=int, default=0, dest='index_offset', help='Added to all particle indices in the parameter file and old-style topology, for proteins that are not first in a hybrid system (default 0)')
    parser.add_argument('--old-top', action='store_true', dest='old_top', help='Also write an old-style topology to <output>_old.top')
    parser.add_argument('--plot', metavar='plot_file', help='Plot the experimental and predicted B-factors to this file')
    parser.add_argument('-q', '--quiet', metavar='quiet', dest='quiet', action='store_const', const=True, default=False, help="Don't print 'INFO' messages to stderr")
    return parser


def main():
    parser = cli_parser(path.basename(__file__))
    args = parser.parse_args()

    logger_settings.set_quiet(args.quiet)
    from oxDNA_analysis_tools.config import check
    check(["python", "numpy", "matplotlib", "sklearn"]) # scipy comes with scikit-learn

    if args.write_kbkt and args.model != 'anmt':
        parser.error("--write-kbkt is only valid with --model anmt")

    outbase = args.output if args.output else path.splitext(args.structure)[0]
    if not args.output:
        log(f"No output name provided, defaulting to \"{outbase}\"")

    structure = read_structure(args.structure)
    bfactors = read_bfactor_file(args.bfactors) if args.bfactors else None

    network = anm_fluctuation_solver(structure, model=args.model, cutoff=args.cutoff, temperature=args.temperature,
                                     bfactors=bfactors, spring_constant=args.spring_constant, fit=args.fit,
                                     max_deviations=args.max_deviations, kb=args.kb, kt=args.kt,
                                     write_kbkt=args.write_kbkt)
    log(f"Spring constant: {network.gamma * 100:.4g} pN/A = {network.k[0]:.6g} oxDNA units")
    if args.spring_constant is None:
        log(f"Correlation between predicted and experimental B-factors: {network.correlation:.3f}")
    if args.model == 'anmt' and not args.write_kbkt:
        log(f"Set bending_k = {args.kb} and torsion_k = {args.kt} in the oxDNA input file")

    write_par(outbase + '.par', network, args.index_offset)
    write_top(outbase + '.top', protein_system(structure))
    write_conf(outbase + '.dat', protein_configuration(structure, network.a1s, network.a3s))
    written = [outbase + ext for ext in ('.par', '.top', '.dat')]
    if args.old_top:
        with open(outbase + '_old.top', 'w') as f:
            f.write(get_old_top_string(structure, args.index_offset))
        written.append(outbase + '_old.top')
    if args.plot:
        plot_bfactors(args.plot, network, structure)
        written.append(args.plot)
    log(f"Wrote {', '.join(written)}")

    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    main()
