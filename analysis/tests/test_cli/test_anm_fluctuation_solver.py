"""
Tests for anm_fluctuation_solver.

The reference parameter files in tests/resources/anm were generated with models.py from the
ANM-oxDNA ANMUtils package ("KDPG Protein Example" notebook: cutoff 13 A, T = 300 K, grid-search
fit of the spring constant). kdpg.pdb is that example's structure (PDB 1WA3, one trimer) reduced
to its ATOM records.
"""
import gzip
import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.anm_fluctuation_solver import (
    anm_fluctuation_solver, read_structure, read_bfactor_file, build_hessian, contact_pairs,
    mean_square_fluctuations, fit_spring_constant, write_par, get_old_top_string, main, cli_parser,
    angular_hessian, eliminate_rotations, BOLTZMANN, BFACTOR_CONVERSION, SIM_FORCE_CONSTANT, SIM_ENERGY,
)
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, strand_describe

RESOURCES = Path(__file__).parent.parent / "resources" / "anm"
KDPG_PDB = RESOURCES / "kdpg.pdb"
N_KDPG = 605
KDPG_CHAIN_LENGTHS = [202, 202, 201]


def read_par(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        lines = [l.split() for l in f.read().splitlines() if l.strip()]
    return int(lines[0][0]), lines[1:]


@pytest.fixture(scope="module")
def kdpg():
    return read_structure(str(KDPG_PDB))


def _pdb_to_mmcif(pdb_path):
    """Minimal mmCIF _atom_site block with the ATOM records of a PDB file."""
    cols = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id',
            'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv',
            'auth_seq_id', 'auth_asym_id', 'pdbx_PDB_model_num']
    out = ['data_test', '#', 'loop_'] + [f'_atom_site.{c}' for c in cols]
    for l in open(pdb_path):
        if not l.startswith('ATOM'):
            continue
        resi = l[22:26].strip()
        out.append(' '.join(['ATOM', l[6:11].strip(), l[76:78].strip() or '?', l[12:16].strip(), l[16].strip() or '.',
                             l[17:20].strip(), l[21], resi, l[26].strip() or '?', l[30:38].strip(), l[38:46].strip(),
                             l[46:54].strip(), l[54:60].strip(), l[60:66].strip(), resi, l[21], '1']))
    return '\n'.join(out) + '\n#\n'


# =============================================================================
# Structure reading
# =============================================================================

class TestReadStructure:

    def test_pdb(self, kdpg):
        assert kdpg.n_residues == N_KDPG
        assert kdpg.chain_ids == ['A', 'B', 'C']
        assert [len(s) for s in kdpg.sequences] == KDPG_CHAIN_LENGTHS
        assert kdpg.sequences[0].startswith('KMEELF')
        assert np.array_equal(kdpg.chain_starts, [0, 202, 404])
        # first CA of the file: ATOM 2 CA LYS A 1 -50.508 22.494 10.403 B=55.93
        ca = [l for l in open(KDPG_PDB) if l.startswith('ATOM') and l[12:16] == ' CA '][0]
        assert np.allclose(kdpg.coords[0], [float(ca[30:38]), float(ca[38:46]), float(ca[46:54])])
        assert kdpg.bfactors[0] == float(ca[60:66])
        # consecutive CAs are 3.8 A apart within chains (2.9-3.0 A for the cis peptide 131-132), > 50 A between chains
        d = np.linalg.norm(np.diff(kdpg.coords, axis=0), axis=1)
        cis = np.array([129, 331, 532])
        assert kdpg.residue_ids[129] == '131' and np.all(np.abs(d[cis] - 2.97) < 0.05)
        within = np.ones(N_KDPG - 1, dtype=bool)
        within[[201, 403]] = False
        within[cis] = False
        assert np.all(np.abs(d[within] - 3.8) < 0.1)
        assert np.all(d[[201, 403]] > 50)

    def test_mmcif_matches_pdb(self, kdpg, tmp_path):
        cif = tmp_path / "kdpg.cif"
        cif.write_text(_pdb_to_mmcif(KDPG_PDB))
        s = read_structure(str(cif))
        assert s.sequences == kdpg.sequences
        assert s.chain_ids == kdpg.chain_ids
        assert np.allclose(s.coords, kdpg.coords)
        assert np.allclose(s.bfactors, kdpg.bfactors)
        assert np.allclose(s.com_vectors, kdpg.com_vectors)
        assert np.allclose(s.n_vectors, kdpg.n_vectors, equal_nan=True) and np.allclose(s.c_vectors, kdpg.c_vectors, equal_nan=True)

    def test_models(self, kdpg, tmp_path):
        atoms = [l for l in open(KDPG_PDB) if l.startswith('ATOM')]
        chain_a = [l for l in atoms if l[21] == 'A']
        # an ensemble repeats the same residues in each model: only the first model is used
        ensemble = tmp_path / "ensemble.pdb"
        ensemble.write_text(''.join(['MODEL        1\n'] + chain_a + ['ENDMDL\n', 'MODEL        2\n'] + chain_a + ['ENDMDL\n']))
        s = read_structure(str(ensemble))
        assert s.sequences == [kdpg.sequences[0]]
        # one chain per model with blank chain ids (as written by Modeller/Chimera): all chains are kept, even
        # identical copies of a homo-oligomer, because they do not overlap in space
        per_chain = ['MODEL        1\n'] + [l[:21] + ' ' + l[22:] for l in chain_a] + ['ENDMDL\n', 'MODEL        2\n'] + \
                    [l[:21] + ' ' + l[22:] for l in atoms if l[21] == 'B'] + ['ENDMDL\n']
        split = tmp_path / "split.pdb"
        split.write_text(''.join(per_chain))
        s = read_structure(str(split))
        assert s.sequences == kdpg.sequences[:2]
        assert s.chain_ids == [' ', '2']
        assert np.allclose(s.coords, kdpg.coords[:404])

    def test_skips_residues_without_ca_and_non_amino_acids(self, kdpg, tmp_path):
        atoms = [l for l in open(KDPG_PDB) if l.startswith('ATOM')]
        removed = kdpg.residue_ids.index('5')  # numbering starts at 2
        no_ca = [l for l in atoms if not (l[21] == 'A' and l[22:26].strip() == '5' and l[12:16] == ' CA ')]
        water = 'ATOM   9999  O   HOH A 999      0.000   0.000   0.000  1.00 10.00           O\n'
        f = tmp_path / "gap.pdb"
        f.write_text(''.join(no_ca) + water)
        s = read_structure(str(f))
        assert s.n_residues == N_KDPG - 1
        assert s.sequences[0] == kdpg.sequences[0][:removed] + kdpg.sequences[0][removed + 1:]
        assert '5' not in s.residue_ids[:201]

    def test_bfactor_file(self, tmp_path):
        f = tmp_path / "b.txt"
        f.write_text("# B-factors\n1.5 2.5\n3.5 # comment\n")
        assert np.array_equal(read_bfactor_file(str(f)), [1.5, 2.5, 3.5])


# =============================================================================
# Network and normal modes
# =============================================================================

class TestNormalModes:

    def test_hessian_and_fluctuations(self, kdpg):
        coords = kdpg.coords[:150]
        i, j = contact_pairs(coords, 13.)
        d = np.linalg.norm(coords[:, None] - coords[None], axis=-1)
        ref_i, ref_j = np.where(np.triu(d <= 13., 1))
        assert np.array_equal(i, ref_i) and np.array_equal(j, ref_j)

        k = np.linspace(0.5, 2., len(i))
        hess = build_hessian(coords, i, j, k)
        # Direct sum of the spring contributions
        ref = np.zeros_like(hess)
        for p in range(len(i)):
            u = (coords[j[p]] - coords[i[p]]) / d[i[p], j[p]]
            b = k[p] * np.outer(u, u)
            for x, y, s in ((i[p], i[p], 1), (j[p], j[p], 1), (i[p], j[p], -1), (j[p], i[p], -1)):
                ref[3*x:3*x+3, 3*y:3*y+3] += s * b
        assert np.allclose(hess, ref)
        # the six rigid-body motions cost no energy
        w = np.linalg.eigvalsh(hess)
        assert np.all(np.abs(w[:6]) < 1e-10) and w[6] > 1e-4

        pinv = np.linalg.pinv(hess, hermitian=True, rcond=1e-10)
        expected = BOLTZMANN * 300. * np.einsum('iaia->i', pinv.reshape(150, 3, 150, 3))
        for method in ('auto', 'cholesky', 'eigh'):
            assert np.allclose(mean_square_fluctuations(hess, coords, 300., method), expected, rtol=1e-8)

    def test_floppy_network(self, kdpg):
        # a residue with only two springs has a zero-energy mode beyond the six rigid-body ones
        coords = np.vstack([kdpg.coords[:100], kdpg.coords[99] + [0., 0., 15.]])
        i, j = contact_pairs(coords, 13.)
        far = j == 100
        i = np.concatenate([i[~far], [98, 99]])
        j = np.concatenate([j[~far], [100, 100]])
        hess = build_hessian(coords, i, j)
        with pytest.raises(np.linalg.LinAlgError):
            mean_square_fluctuations(hess, coords, method='cholesky')
        auto = mean_square_fluctuations(hess, coords, method='auto')
        assert np.allclose(auto, mean_square_fluctuations(hess, coords, method='eigh'))
        # whether round-off leaves the zero mode slightly positive depends on the BLAS; a tiny shift makes the
        # factorisation succeed everywhere, and the singular matrix must still be detected
        shifted = hess + 1e-13 * np.trace(hess) / len(hess) * np.eye(len(hess))
        with pytest.raises(np.linalg.LinAlgError):
            mean_square_fluctuations(shifted, coords, method='cholesky')
        assert np.allclose(mean_square_fluctuations(shifted, coords, method='auto'), auto)

    def test_linear_fit(self):
        rng = np.random.default_rng(1)
        msf = rng.uniform(0.5, 2., 200)
        exp_b = BFACTOR_CONVERSION * msf / 0.4
        exp_b[[3, 50]] *= 5  # noise on the experiment does not move the fit much
        msf[7] = 40.           # a floppy residue is excluded
        gamma, used = fit_spring_constant(msf, exp_b)
        assert not used[7] and used.sum() == 199
        assert gamma == pytest.approx(1 / (np.dot(msf[used], exp_b[used]) / np.dot(msf[used], msf[used]) / BFACTOR_CONVERSION))
        gamma, used = fit_spring_constant(msf, exp_b, method='grid')
        assert used.all() and 0.001 <= gamma < 5


# =============================================================================
# ANMT: bending and torsion in the normal mode analysis
# =============================================================================

def _rotate(theta, v):
    angle = np.linalg.norm(theta)
    if angle == 0:
        return v
    k = theta / angle
    return v * np.cos(angle) + np.cross(k, v) * np.sin(angle) + k * np.dot(k, v) * (1 - np.cos(angle))


class TestANMT:

    def test_angular_hessian_matches_finite_differences(self, kdpg):
        """The analytic Hessian of the engine's _protein_ang_pot, with frame rotations theta (da = theta x a)."""
        n = 6
        x = kdpg.coords[:n]
        from oxDNA_analysis_tools.anm_fluctuation_solver import orientation_vectors
        a1, a3 = (v[:n] for v in orientation_vectors(kdpg))
        kb, kt = 2., 3.
        hxx, hxt, htt = angular_hessian(x, np.arange(n - 1), a1, a3, kb, kt)
        full = np.block([[hxx, hxt], [hxt.T, htt]])

        def energy(dof):
            pos, th = x + dof[:3 * n].reshape(n, 3), dof[3 * n:].reshape(n, 3)
            e = 0.
            for p in range(n - 1):
                q = p + 1
                r0 = (x[q] - x[p]) / np.linalg.norm(x[q] - x[p])
                r = (pos[q] - pos[p]) / np.linalg.norm(pos[q] - pos[p])
                ap, bq = _rotate(th[p], a1[p]), _rotate(th[q], a1[q])
                Ap, Bq = _rotate(th[p], a3[p]), _rotate(th[q], a3[q])
                e += kb / 2 * ((r @ ap - r0 @ a1[p])**2 + (-r @ bq + r0 @ a1[q])**2)
                e += kt / 2 * ((ap @ bq - a1[p] @ a1[q])**2 + (Ap @ Bq - a3[p] @ a3[q])**2)
            return e

        h = 1e-4
        m = 6 * n
        fd = np.zeros((m, m))
        for a in range(m):
            for b in range(a, m):
                ea, eb = np.zeros(m), np.zeros(m)
                ea[a], eb[b] = h, h
                fd[a, b] = fd[b, a] = (energy(ea + eb) - energy(ea - eb) - energy(eb - ea) + energy(-ea - eb)) / (4 * h * h)
        assert np.allclose(full, fd, atol=1e-5 * np.abs(full).max())

        # Schur complement against a dense pseudo-inverse
        hs = build_hessian(x, *contact_pairs(x, 13.))
        heff = eliminate_rotations(hs + hxx, hxt, htt)
        ref = hs + hxx - hxt @ np.linalg.pinv(htt, rcond=1e-12, hermitian=True) @ hxt.T
        assert np.allclose(heff, ref, atol=1e-8 * np.abs(ref).max())

    def test_oxview_reference(self, kdpg):
        """
        Fluctuations computed independently by the oxView solver's reference script for chain A of KDPG (202 beads, frames
        from orientation_vectors): H = sum k g g^T with central-difference gradients, rotations eliminated with a pseudo-inverse.
        """
        from oxDNA_analysis_tools.anm_fluctuation_solver import orientation_vectors
        with open(RESOURCES / "oxview_kdpgA_ref.json") as f:
            ref = json.load(f)
        x = np.array(ref['coords']).reshape(-1, 3) * 8.518
        a1 = np.array(ref['a1']).reshape(-1, 3)
        a3 = np.array(ref['a3']).reshape(-1, 3)
        n = len(x)
        # the reference uses the frames of this module
        oat_a1, oat_a3 = orientation_vectors(kdpg)
        assert np.allclose(x, kdpg.coords[:n], atol=1e-6)
        assert np.allclose(a1, oat_a1[:n], atol=1e-12) and np.allclose(a3, oat_a3[:n], atol=1e-12)
        i, j = contact_pairs(x, 13.)
        assert np.array_equal(i, ref['p1']) and np.array_equal(j, ref['p2'])
        backbone = np.array([p for p, _ in ref['pairs']])
        assert np.array_equal(backbone, i[j == i + 1])

        springs = build_hessian(x, i, j, SIM_FORCE_CONSTANT)  # k = 1 oxDNA unit
        assert np.allclose(mean_square_fluctuations(springs, x, 300.), ref['msf_anm_1'], rtol=1e-8)
        hxx, hxt, htt = angular_hessian(x, backbone, a1, a3, 50 * SIM_ENERGY, 50 * SIM_ENERGY)
        angular = eliminate_rotations(hxx, hxt, htt)
        for scale, key in ((1., 'msf_anmt_1'), (3.7, 'msf_anmt_37')):
            msf = mean_square_fluctuations(springs * scale + angular, x, 300.)
            # the frames that follow the backbone give soft rotational modes (~1e-10 of the largest), which a large
            # regularisation of H_tt would stiffen: with the smallest shift the two solvers agree to ~1e-5
            assert np.allclose(msf, ref[key], rtol=2e-5)
            assert np.median(np.abs(msf / np.array(ref[key]) - 1)) < 1e-7

    def test_orientation_vectors(self, kdpg):
        from oxDNA_analysis_tools.anm_fluctuation_solver import orientation_vectors
        a1, a3 = orientation_vectors(kdpg)
        assert np.allclose(np.linalg.norm(a1, axis=1), 1) and np.allclose(np.linalg.norm(a3, axis=1), 1)
        assert np.allclose(np.einsum('ij,ij->i', a1, a3), 0)
        # a1 points to the center of mass of the side chain (CA included, backbone N, C, O excluded)
        res = [l for l in open(KDPG_PDB) if l.startswith('ATOM') and l[21] == 'A' and l[22:26].strip() == kdpg.residue_ids[0]]
        xyz = {l[12:16].strip(): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in res}
        com = np.mean([v for k, v in xyz.items() if k not in ('N', 'C', 'O')], axis=0)
        assert np.allclose(a1[0], (com - xyz['CA']) / np.linalg.norm(com - xyz['CA']))
        # a3 is the N -> C direction of the residue, made perpendicular to a1
        nc = xyz['C'] - xyz['N']
        nc_perp = nc - np.dot(nc, a1[0]) * a1[0]
        assert np.allclose(a3[0], nc_perp / np.linalg.norm(nc_perp))
        # and therefore it follows the chain: it points towards the next C-alpha
        x = kdpg.coords
        chain = np.repeat(np.arange(len(kdpg.sequences)), [len(q) for q in kdpg.sequences])
        same = chain[1:] == chain[:-1]
        along = np.einsum('ij,ij->i', a3[:-1][same], x[1:][same] - x[:-1][same])
        assert np.all(along > 0)
        # glycines have no side chain: a1 points to the ideal C-beta, which is close to where the real C-beta of the
        # other residues is
        seq = ''.join(kdpg.sequences)
        gly = [r for r in range(kdpg.n_residues) if seq[r] == 'G']
        assert len(gly) > 0
        for r in gly:
            nv, cv = kdpg.n_vectors[r], kdpg.c_vectors[r]
            assert np.dot(a1[r], np.cross(-nv, cv)) < 0  # the C-beta is on the opposite side of the N-CA-C normal
        ideal = []
        for r in range(kdpg.n_residues):
            if seq[r] != 'G':
                b, c = -kdpg.n_vectors[r], kdpg.c_vectors[r]
                cb = -0.58273431 * np.cross(b, c) + 0.56802827 * b - 0.54067466 * c
                ideal.append(np.dot(cb / np.linalg.norm(cb), a1[r]))
        assert np.median(ideal) > 0.8  # the side chain points roughly the way of the C-beta

    def test_frames_do_not_depend_on_the_input_orientation(self, kdpg):
        # rotating the structure must rotate the frames with it and leave every parameter unchanged
        from dataclasses import replace
        rng = np.random.default_rng(7)
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(q) < 0:
            q[:, 0] *= -1
        rotated = replace(kdpg, coords=kdpg.coords @ q.T + np.array([3., -7., 11.]), com_vectors=kdpg.com_vectors @ q.T,
                          n_vectors=kdpg.n_vectors @ q.T, c_vectors=kdpg.c_vectors @ q.T)
        a = anm_fluctuation_solver(kdpg, model='anmt', spring_constant=5.)
        b = anm_fluctuation_solver(rotated, model='anmt', spring_constant=5.)
        assert np.allclose(b.a1s, a.a1s @ q.T, atol=1e-10) and np.allclose(b.a3s, a.a3s @ q.T, atol=1e-10)
        assert np.allclose(b.angles, a.angles, atol=1e-10)
        # the B-factors only differ by round-off, amplified (up to ~3e-5) by the nearly singular rotation block
        assert np.allclose(b.r0, a.r0) and np.allclose(b.pred_bfactors, a.pred_bfactors, rtol=1e-4)
        assert np.median(np.abs(b.pred_bfactors / a.pred_bfactors - 1)) < 1e-7

    def test_anmt_fit(self, kdpg):
        anm = anm_fluctuation_solver(kdpg, model='anm')
        anmt = anm_fluctuation_solver(kdpg, model='anmt', kb=50., kt=50.)
        # the angular terms stiffen the network, so the fitted springs are softer
        assert anmt.gamma < anm.gamma
        assert anmt.kb == 50. and anmt.kt == 50. and not anmt.write_kbkt and anm.kb is None
        assert np.array_equal(anmt.fitted, anm.fitted)
        for net in (anm, anmt):
            m, e = net.pred_bfactors[net.fitted], net.exp_bfactors[net.fitted]
            assert np.dot(e, m) / np.dot(m, m) == pytest.approx(1., abs=1e-6)
        # the angular contribution to the predicted B-factors is not a uniform rescaling
        ratio = anmt.pred_bfactors / anm.pred_bfactors
        assert np.ptp(ratio) > 0.05
        # without bending and torsion ANMT is the ANM
        flat = anm_fluctuation_solver(kdpg, model='anmt', kb=0., kt=0.)
        assert flat.gamma == pytest.approx(anm.gamma) and np.allclose(flat.pred_bfactors, anm.pred_bfactors)
        # a fixed spring constant gives the same B-factors as the fit at that value
        fixed = anm_fluctuation_solver(kdpg, model='anmt', spring_constant=anmt.k[0])
        assert np.allclose(fixed.pred_bfactors, anmt.pred_bfactors, rtol=1e-8)


# =============================================================================
# Reproduction of the ANM-oxDNA reference files
# =============================================================================

class TestReference:

    @pytest.mark.parametrize("model,ref", [("anm", "kdpg_ref.par.gz"), ("anmt", "kdpgt_ref.par.gz")])
    def test_reference_par(self, kdpg, tmp_path, model, ref):
        net = anm_fluctuation_solver(kdpg, model=model, cutoff=13., temperature=300., fit='grid')
        assert net.gamma == pytest.approx(0.424)
        out = tmp_path / "out.par"
        write_par(str(out), net)
        n_ref, ref_rows = read_par(RESOURCES / ref)
        n_out, out_rows = read_par(out)
        assert n_out == n_ref == N_KDPG
        assert len(out_rows) == len(ref_rows) == 12716
        assert [r[:2] for r in out_rows] == [r[:2] for r in ref_rows]
        assert [len(r) for r in out_rows] == [len(r) for r in ref_rows]
        assert all(r[3] == 's' for r in out_rows)
        ref_num = np.array([[float(r[2]), float(r[4])] for r in ref_rows])
        out_num = np.array([[float(r[2]), float(r[4])] for r in out_rows])
        # the reference used float32 coordinates (Biopython)
        assert np.allclose(out_num[:, 0], ref_num[:, 0], rtol=1e-5)
        assert np.allclose(out_num[:, 1], ref_num[:, 1], rtol=1e-12)
        if model == 'anmt':
            # the cosines are not compared: the bead frames are not those of models.py (see orientation_vectors)
            out_ang = np.array([[float(x) for x in r[5:]] for r in out_rows if len(r) > 5])
            assert len(out_ang) == 602  # every j == i + 1 spring (the chain ends are too far apart for a spring)
            assert np.all(np.abs(out_ang) <= 1)

    def test_linear_fit_matches_models_py(self, kdpg):
        # values from models.py ANM.calc_ANM_unitary() (linear fit) on the same structure.  models.py used
        # kB = 0.00138064852 (N/m) A^2 / K, here kT = T / 3000 oxDNA energy units, which is 1e-4 larger
        net = anm_fluctuation_solver(kdpg, model='anm', cutoff=13.)
        assert net.gamma == pytest.approx(0.389629155, rel=2e-4)
        assert np.where(~net.fitted)[0].tolist() == [0, 55, 58, 59, 100, 185, 200, 201, 202, 257, 260, 261, 302,
                                                     385, 387, 402, 458, 461, 462, 503, 588, 603, 604]
        assert np.allclose(net.k, net.gamma / SIM_FORCE_CONSTANT)
        assert 0.6 < net.correlation < 0.7
        assert np.allclose(net.rmsf, np.sqrt(net.pred_bfactors / BFACTOR_CONVERSION))


# =============================================================================
# Consistency of the outputs with the oxDNA engine
# =============================================================================

def _engine_cosines(conf, i, j):
    """Angular deviations as computed by DNANMInteraction::_protein_ang_pot after the engine reads a configuration."""
    a1 = conf.a1s / np.linalg.norm(conf.a1s, axis=1)[:, None]
    a3 = conf.a3s / np.linalg.norm(conf.a3s, axis=1)[:, None]
    a1 = a1 - a3 * np.einsum('ij,ij->i', a1, a3)[:, None]
    a1 /= np.linalg.norm(a1, axis=1)[:, None]
    r = conf.positions[j] - conf.positions[i]
    r /= np.linalg.norm(r, axis=1)[:, None]
    return np.stack([np.einsum('ij,ij->i', r, a1[i]), np.einsum('ij,ij->i', -r, a1[j]),
                     np.einsum('ij,ij->i', a1[i], a1[j]), np.einsum('ij,ij->i', a3[i], a3[j])], axis=1)


class TestCLI:

    def test_cli_outputs(self, kdpg, tmp_path):
        out = tmp_path / "kdpg"
        plot = tmp_path / "bfactors.png"
        argv = ['anm_fluctuation_solver.py', str(KDPG_PDB), '-o', str(out), '-m', 'anmt', '--kb', '1.5', '--kt', '2.5', '--write-kbkt',
                '-i', '100', '--old-top', '--plot', str(plot), '-q']
        with patch.object(sys, 'argv', argv):
            main()

        # parameter file
        n, rows = read_par(str(out) + '.par')
        assert n == N_KDPG
        idx = np.array([[int(r[0]), int(r[1])] for r in rows])
        assert idx.min() == 100 and idx.max() == 100 + N_KDPG - 1
        assert np.all(idx[:, 0] < idx[:, 1])
        backbone = idx[:, 1] - idx[:, 0] == 1
        assert all(len(r) == 11 and r[9:] == ['1.5', '2.5'] for r, b in zip(rows, backbone) if b)
        assert all(len(r) == 5 for r, b in zip(rows, backbone) if not b)
        net = anm_fluctuation_solver(kdpg, model='anmt', kb=1.5, kt=2.5)
        assert float(rows[0][4]) == pytest.approx(net.gamma / SIM_FORCE_CONSTANT)

        # new-style topology
        top = (Path(str(out) + '.top')).read_text().splitlines()
        assert top[0] == f'{N_KDPG} 3 0 {N_KDPG} 0 5->3'
        assert [l.split()[0] for l in top[1:]] == kdpg.sequences
        assert all('type=peptide' in l for l in top[1:])
        system, _ = strand_describe(str(out) + '.top')
        assert [s.get_sequence() for s in system.strands] == kdpg.sequences
        assert all(s.is_peptide() for s in system.strands)

        # old-style topology: <strand> <aa> <nside> <cside> only
        old = (Path(str(out) + '_old.top')).read_text().splitlines()
        assert old[0] == f'{N_KDPG} 3 0 {N_KDPG} 0'
        cols = [l.split() for l in old[1:]]
        assert all(len(c) == 4 for c in cols)
        assert [c[0] for c in cols] == ['-1'] * 202 + ['-2'] * 202 + ['-3'] * 201
        assert ''.join(c[1] for c in cols) == ''.join(kdpg.sequences)
        assert cols[0][2:] == ['-1', '101'] and cols[1][2:] == ['100', '102']
        assert cols[201][2:] == ['300', '-1'] and cols[202][2:] == ['-1', '303']
        assert old == get_old_top_string(kdpg, 100).splitlines()

        # configuration: same geometry as the structure, and the angular potential is at its minimum
        top_info, traj_info = describe(str(out) + '.top', str(out) + '.dat')
        conf = get_confs(top_info, traj_info, 0, 1)[0]
        assert conf.positions.shape == (N_KDPG, 3)
        assert np.allclose(np.linalg.norm(conf.positions[1] - conf.positions[0]) * 8.518,
                           np.linalg.norm(kdpg.coords[1] - kdpg.coords[0]))
        assert np.all(conf.box > np.ptp(conf.positions, axis=0))
        assert np.allclose(np.linalg.norm(conf.a1s, axis=1), 1) and np.allclose(np.linalg.norm(conf.a3s, axis=1), 1)
        assert np.allclose(np.einsum('ij,ij->i', conf.a1s, conf.a3s), 0, atol=1e-12)
        bb = idx[backbone] - 100
        written = np.array([[float(x) for x in r[5:9]] for r, b in zip(rows, backbone) if b])
        assert np.allclose(_engine_cosines(conf, bb[:, 0], bb[:, 1]), written, atol=1e-8)

        assert plot.exists() and plot.stat().st_size > 10000

    def test_fixed_spring_constant_and_bfactor_file(self, kdpg, tmp_path):
        bfile = tmp_path / "b.txt"
        bfile.write_text('\n'.join(str(b) for b in kdpg.bfactors * 2))
        out = tmp_path / "fixed"
        with patch.object(sys, 'argv', ['anm_fluctuation_solver.py', str(KDPG_PDB), '-o', str(out), '-k', '5', '-c', '10', '-q']):
            main()
        n, rows = read_par(str(out) + '.par')
        assert all(float(r[4]) == pytest.approx(5.) for r in rows) and all(len(r) == 5 for r in rows)
        assert all(float(r[2]) * 8.518 <= 10. for r in rows)

        # doubling the experimental B-factors halves the fitted spring constant
        ref = anm_fluctuation_solver(kdpg)
        doubled = anm_fluctuation_solver(kdpg, bfactors=read_bfactor_file(str(bfile)))
        assert doubled.gamma == pytest.approx(ref.gamma / 2)

    def test_errors(self, kdpg):
        with pytest.raises(RuntimeError):
            anm_fluctuation_solver(kdpg, model='anmt', kb=-1.)
        with pytest.raises(RuntimeError):
            anm_fluctuation_solver(kdpg, model='anm', write_kbkt=True)
        with pytest.raises(RuntimeError):
            anm_fluctuation_solver(kdpg, bfactors=np.ones(N_KDPG))
        with pytest.raises(RuntimeError):
            anm_fluctuation_solver(kdpg, bfactors=np.ones(10))
        with pytest.raises(RuntimeError):
            anm_fluctuation_solver(kdpg, model='hanm')

    def test_parser(self):
        args = cli_parser().parse_args(['x.pdb'])
        assert args.model == 'anm' and args.cutoff == 13. and args.temperature == 300. and args.fit == 'linear'
        assert args.index_offset == 0 and not args.old_top and args.plot is None
        assert args.kb == 50. and args.kt == 50. and not args.write_kbkt
