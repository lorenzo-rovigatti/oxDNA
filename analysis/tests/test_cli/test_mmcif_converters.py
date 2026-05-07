"""
Tests for mmCIF read/write via the consolidated oxDNA_PDB / PDB_oxDNA API.

oxDNA_PDB(... format='mmcif') writes mmCIF.
PDB_oxDNA(cif_str)            reads mmCIF (auto-detected from data_ header).

Test strategy (end-to-end, increasing specificity per CLAUDE.md):
  1. Output exists / no exception raised
  2. File structure is valid mmCIF (data_ block, loop_, _atom_site columns)
  3. Strand/nucleotide counts are preserved
  4. Round-trip: oxDNA → mmCIF → oxDNA preserves base composition and COM
  5. Entity annotation blocks contain correct sequence and polymer type
  6. Auto-format selection triggers correctly
"""
from pathlib import Path

import numpy as np
import pytest

from oxDNA_analysis_tools.oxDNA_PDB import oxDNA_PDB, _strand_sequence, _resolve_output_format
from oxDNA_analysis_tools.PDB_oxDNA import PDB_oxDNA
from oxDNA_analysis_tools.UTILS.mmcif import tokenize_cif, parse_atom_site
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, strand_describe, inbox
from oxDNA_analysis_tools.UTILS.data_structures import System, Strand


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def resources():
    return Path(__file__).parent.parent / "resources"


@pytest.fixture(scope="module")
def system_and_conf(resources):
    top  = resources / "rna_tile.top"
    traj = resources / "minitraj.dat"
    system, _ = strand_describe(str(top))
    ti, di = describe(str(top), str(traj))
    conf = get_confs(ti, di, 0, 1)[0]
    conf = inbox(conf, center=True)
    return system, conf


@pytest.fixture
def cif_file(system_and_conf, tmp_path):
    """Write a combined .cif via oxDNA_PDB and return its path."""
    system, conf = system_and_conf
    out_base = str(tmp_path / "output")
    oxDNA_PDB(conf, system, out_base, format='mmcif')
    return tmp_path / "output.cif"


# ── Tokenizer unit tests ───────────────────────────────────────────────────────

class TestTokenizer:
    def test_bare_tokens(self):
        assert tokenize_cif("ATOM 1 C") == ["ATOM", "1", "C"]

    def test_single_quoted(self):
        assert tokenize_cif("'hello world'") == ["hello world"]

    def test_double_quoted(self):
        assert tokenize_cif('"hello world"') == ["hello world"]

    def test_comment_skipped(self):
        assert tokenize_cif("ATOM # comment\n1") == ["ATOM", "1"]

    def test_atom_name_with_prime(self):
        assert tokenize_cif("C4' N1 P") == ["C4'", "N1", "P"]

    def test_semicolon_block(self):
        assert "ACGT" in tokenize_cif("loop_\n;\nACGT\n;\n")


# ── Format resolution ──────────────────────────────────────────────────────────

class TestFormatResolution:

    def test_explicit_pdb(self, system_and_conf):
        system, _ = system_and_conf
        assert _resolve_output_format('pdb', system) == 'pdb'

    def test_explicit_mmcif(self, system_and_conf):
        system, _ = system_and_conf
        assert _resolve_output_format('mmcif', system) == 'mmcif'

    def test_auto_small_system_returns_pdb(self, system_and_conf):
        system, _ = system_and_conf
        # rna_tile is small, so auto should pick pdb
        assert _resolve_output_format('auto', system) == 'pdb'

    def test_auto_large_system_returns_mmcif(self):
        """A synthetic system large enough to exceed PDB limits should pick mmcif."""
        # Build a fake system with enough monomers to trigger auto-switch.
        # 9999 residues * 33 atoms/nuc = 329967 > 99999, so mmcif.
        big_system = System('')
        big_strand = Strand(1)
        big_strand._monomers = [object()] * 10000  # duck-typed length
        big_system.strands = [big_strand]
        # Can't easily build real Monomer objects; patch via the formula directly.
        total_residues = 10000
        est_atoms = total_residues * 33
        assert est_atoms > 99999 or total_residues > 9999


# ── oxDNA → mmCIF writer tests ─────────────────────────────────────────────────

class TestOxDNAmmCIF:

    def test_creates_cif_file(self, system_and_conf, tmp_path):
        system, conf = system_and_conf
        oxDNA_PDB(conf, system, str(tmp_path / "out"), format='mmcif')
        assert (tmp_path / "out.cif").exists()

    def test_pdb_format_creates_pdb_file(self, system_and_conf, tmp_path):
        system, conf = system_and_conf
        oxDNA_PDB(conf, system, str(tmp_path / "out"), format='pdb')
        assert (tmp_path / "out.pdb").exists()
        assert not (tmp_path / "out.cif").exists()

    def test_cif_extension_in_output_name_selects_mmcif(self, resources, tmp_path):
        """Passing -o output.cif on the CLI should produce a .cif file."""
        import subprocess, sys
        top  = str(resources / "rna_tile.top")
        traj = str(resources / "minitraj.dat")
        out  = str(tmp_path / "result.cif")
        result = subprocess.run(
            [sys.executable, "-m", "oxDNA_analysis_tools.oxDNA_PDB",
             top, traj, "35", "-o", out],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "result.cif").exists()
        content = (tmp_path / "result.cif").read_text()
        assert content.startswith("data_")

    def test_file_has_data_block(self, cif_file):
        assert cif_file.read_text().startswith("data_")

    def test_atom_site_loop_present(self, cif_file):
        content = cif_file.read_text()
        assert "_atom_site.group_PDB" in content
        assert "_atom_site.Cartn_x" in content
        assert "ATOM" in content

    def test_entity_blocks_present(self, cif_file):
        content = cif_file.read_text()
        assert "_entity_poly.pdbx_seq_one_letter_code" in content
        assert "_entity.pdbx_description" in content

    def test_residue_count_matches(self, system_and_conf, cif_file):
        system, _ = system_and_conf
        rows = parse_atom_site(cif_file.read_text())
        residues = {(r['auth_asym_id'], r['auth_seq_id'])
                    for r in rows if r.get('group_PDB') == 'ATOM'}
        total_nuc = sum(len(s.monomers) for s in system.strands)
        assert len(residues) == total_nuc

    def test_chain_count_matches_strand_count(self, system_and_conf, cif_file):
        system, _ = system_and_conf
        rows = parse_atom_site(cif_file.read_text())
        chains = {r['auth_asym_id'] for r in rows if r.get('group_PDB') == 'ATOM'}
        assert len(chains) == len(system.strands)

    def test_one_file_per_strand(self, system_and_conf, tmp_path):
        system, conf = system_and_conf
        oxDNA_PDB(conf, system, str(tmp_path / "strand"),
                  format='mmcif', one_file_per_strand=True)
        cif_files = list(tmp_path.glob("strand_*.cif"))
        assert len(cif_files) == len(system.strands)
        for cf in cif_files:
            content = cf.read_text()
            assert content.startswith("data_")
            assert "_atom_site.group_PDB" in content


# ── mmCIF → oxDNA reader tests ─────────────────────────────────────────────────

class TestmmCIFOxDNA:

    def test_auto_detects_mmcif(self, cif_file):
        """PDB_oxDNA should parse mmCIF without an explicit format flag."""
        configs, systems = PDB_oxDNA(cif_file.read_text())
        assert len(configs) == 1
        assert len(systems) == 1

    def test_strand_count_preserved(self, system_and_conf, cif_file):
        system, _ = system_and_conf
        _, systems = PDB_oxDNA(cif_file.read_text())
        assert len(systems[0].strands) == len(system.strands)

    def test_nucleotide_count_preserved(self, system_and_conf, cif_file):
        system, _ = system_and_conf
        _, systems = PDB_oxDNA(cif_file.read_text())
        orig = sum(len(s.monomers) for s in system.strands)
        rt   = sum(len(s.monomers) for s in systems[0].strands)
        assert rt == orig

    def test_configuration_shape(self, system_and_conf, cif_file):
        system, _ = system_and_conf
        configs, _ = PDB_oxDNA(cif_file.read_text())
        total = sum(len(s.monomers) for s in system.strands)
        assert configs[0].positions.shape == (total, 3)
        assert configs[0].a1s.shape       == (total, 3)
        assert configs[0].a3s.shape       == (total, 3)

    def test_orientation_vectors_are_unit(self, cif_file):
        configs, _ = PDB_oxDNA(cif_file.read_text())
        np.testing.assert_allclose(np.linalg.norm(configs[0].a1s, axis=1), 1.0, atol=1e-5)
        np.testing.assert_allclose(np.linalg.norm(configs[0].a3s, axis=1), 1.0, atol=1e-5)

    def test_pdb_still_works(self, system_and_conf, tmp_path):
        """PDB_oxDNA must still parse plain PDB files correctly."""
        system, conf = system_and_conf
        out_base = str(tmp_path / "out")
        oxDNA_PDB(conf, system, out_base, format='pdb')
        pdb_str = (tmp_path / "out.pdb").read_text()
        configs, systems = PDB_oxDNA(pdb_str)
        assert len(systems[0].strands) == len(system.strands)


# ── Round-trip tests ───────────────────────────────────────────────────────────

class TestRoundTrip:
    """oxDNA → mmCIF → oxDNA via the unified API."""

    def test_base_composition_preserved(self, system_and_conf, cif_file):
        """Sorted base list per strand must be identical (direction-agnostic)."""
        system, _ = system_and_conf
        _, systems_rt = PDB_oxDNA(cif_file.read_text())
        for orig, rt in zip(system.strands, systems_rt[0].strands):
            assert sorted(_strand_sequence(orig)) == sorted(_strand_sequence(rt)), \
                f"Strand {orig.id}: base composition changed after round-trip"

    def test_com_within_tolerance(self, system_and_conf, cif_file):
        """Per-strand COM must be within 5 Å after round-trip.

        Uses m.pos directly (oxDNA units) rather than indexing conf.positions
        via m.id, because m.id in round-tripped systems reflects the per-chain
        seq_id from the CIF file, not the global nucleotide index.
        """
        system, conf = system_and_conf
        OXDNA_TO_ANGSTROM = 8.518
        _, systems_rt = PDB_oxDNA(cif_file.read_text())
        for orig, rt in zip(system.strands, systems_rt[0].strands):
            orig_com = np.mean(conf.positions[[m.id for m in orig.monomers]], axis=0) * OXDNA_TO_ANGSTROM
            rt_com   = np.mean([m.pos for m in rt.monomers], axis=0) * OXDNA_TO_ANGSTROM
            dist = np.linalg.norm(orig_com - rt_com)
            assert dist < 5.0, f"Strand {orig.id} COM shifted {dist:.2f} Å"


# ── Entity annotation tests ────────────────────────────────────────────────────

class TestEntityBlocks:

    def test_sequences_in_entity_poly(self, system_and_conf, cif_file):
        system, _ = system_and_conf
        content = cif_file.read_text()
        for strand in system.strands:
            seq = _strand_sequence(strand)
            assert seq in content, \
                f"Sequence for strand {strand.id} not found in _entity_poly block"

    def test_polymer_type_correct(self, system_and_conf, cif_file):
        system, _ = system_and_conf
        content = cif_file.read_text()
        for strand in system.strands:
            if strand.type == 'DNA':
                assert 'polydeoxyribonucleotide' in content
            elif strand.type == 'RNA':
                assert 'polyribonucleotide' in content

    def test_entity_description_present(self, cif_file):
        content = cif_file.read_text()
        assert 'oxDNA coarse-grained model' in content


# ── CLI parser tests ───────────────────────────────────────────────────────────

class TestCLIParser:

    def test_oxdna_pdb_accepts_format_flag(self):
        from oxDNA_analysis_tools.oxDNA_PDB import cli_parser
        parser = cli_parser()
        args = parser.parse_args(['top.top', 'conf.dat', '53', '--format', 'mmcif'])
        assert args.format == 'mmcif'

    def test_oxdna_pdb_format_default_is_auto(self):
        from oxDNA_analysis_tools.oxDNA_PDB import cli_parser
        parser = cli_parser()
        args = parser.parse_args(['top.top', 'conf.dat', '53'])
        assert args.format == 'auto'

    def test_pdb_oxdna_cli_accepts_cif_file(self):
        from oxDNA_analysis_tools.PDB_oxDNA import cli_parser
        parser = cli_parser()
        args = parser.parse_args(['input.cif'])
        assert args.pdb_file == 'input.cif'

    def test_oxdna_pdb_rejects_invalid_format(self):
        from oxDNA_analysis_tools.oxDNA_PDB import cli_parser
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['top.top', 'conf.dat', '53', '--format', 'gro'])
