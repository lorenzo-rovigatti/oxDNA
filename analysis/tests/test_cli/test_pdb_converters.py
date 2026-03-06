"""
Tests for oxDNA_PDB and PDB_oxDNA modules.

Tests cover:
- oxDNA_PDB() API function
- PDB_oxDNA() API function
- Round-trip conversion (oxDNA -> PDB -> oxDNA)
- CLI argument parsing for both modules
"""
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.oxDNA_PDB import oxDNA_PDB, _format_atom_serial, cli_parser as oxDNA_PDB_cli_parser
from oxDNA_analysis_tools.PDB_oxDNA import PDB_oxDNA, _parse_atom_serial, cli_parser as PDB_oxDNA_cli_parser
from oxDNA_analysis_tools.UTILS.RyeReader import (
    describe,
    get_confs,
    strand_describe,
    inbox,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def test_resources():
    """Get the path to test resources directory."""
    return Path(__file__).parent.parent / "resources"


@pytest.fixture(scope="module")
def topology_path(test_resources):
    """Path to the RNA tile topology file."""
    return test_resources / "rna_tile.top"


@pytest.fixture(scope="module")
def trajectory_path(test_resources):
    """Path to the mini trajectory file."""
    return test_resources / "minitraj.dat"


@pytest.fixture(scope="module")
def system_and_conf(topology_path, trajectory_path):
    """Load system and first configuration from test files."""
    system, _ = strand_describe(str(topology_path))
    top_info, traj_info = describe(str(topology_path), str(trajectory_path))
    conf = get_confs(top_info, traj_info, 0, 1)[0]
    conf = inbox(conf, center=True)
    return system, conf


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# oxDNA_PDB Tests
# =============================================================================

class TestOxDNAPDB:
    """Tests for the oxDNA_PDB() function."""

    def test_oxdna_pdb_creates_file(self, system_and_conf, temp_output_dir):
        """Test that oxDNA_PDB creates a PDB file."""
        system, conf = system_and_conf
        out_basename = str(temp_output_dir / "output")

        oxDNA_PDB(conf, system, out_basename)

        pdb_file = temp_output_dir / "output.pdb"
        assert pdb_file.exists(), "PDB file should be created"

        # Check file has content
        content = pdb_file.read_text()
        assert len(content) > 0, "PDB file should not be empty"
        assert "ATOM" in content, "PDB file should contain ATOM records"
        assert "END" in content, "PDB file should contain END record"

    def test_oxdna_pdb_contains_all_residues(self, system_and_conf, temp_output_dir):
        """Test that the PDB file contains all expected residues."""
        system, conf = system_and_conf
        out_basename = str(temp_output_dir / "output")

        oxDNA_PDB(conf, system, out_basename)

        pdb_file = temp_output_dir / "output.pdb"
        content = pdb_file.read_text()

        # Count TER records (one per strand)
        ter_count = content.count("TER")
        assert ter_count == len(system.strands), \
            f"Should have {len(system.strands)} TER records, got {ter_count}"

    def test_oxdna_pdb_one_file_per_strand(self, system_and_conf, temp_output_dir):
        """Test one_file_per_strand option creates separate files."""
        system, conf = system_and_conf
        out_basename = str(temp_output_dir / "output")

        oxDNA_PDB(conf, system, out_basename, one_file_per_strand=True)

        # Check that files were created for each strand
        for strand in system.strands:
            strand_file = temp_output_dir / f"output_{strand.id}.pdb"
            assert strand_file.exists(), f"File for strand {strand.id} should exist"


# =============================================================================
# PDB_oxDNA Tests
# =============================================================================

class TestPDBOxDNA:
    """Tests for the PDB_oxDNA() function."""

    def test_pdb_oxdna_returns_correct_types(self, system_and_conf, temp_output_dir):
        """Test that PDB_oxDNA returns correct types."""
        system, conf = system_and_conf
        out_basename = str(temp_output_dir / "output")

        # First create a PDB file
        oxDNA_PDB(conf, system, out_basename)
        pdb_file = temp_output_dir / "output.pdb"

        # Convert back
        with open(pdb_file) as f:
            pdb_str = f.read()

        configs, systems = PDB_oxDNA(pdb_str)

        assert isinstance(configs, list), "Should return list of configs"
        assert isinstance(systems, list), "Should return list of systems"
        assert len(configs) == 1, "Should have one config"
        assert len(systems) == 1, "Should have one system"

    def test_pdb_oxdna_preserves_strand_count(self, system_and_conf, temp_output_dir):
        """Test that conversion preserves strand count."""
        system, conf = system_and_conf
        out_basename = str(temp_output_dir / "output")

        oxDNA_PDB(conf, system, out_basename)
        pdb_file = temp_output_dir / "output.pdb"

        with open(pdb_file) as f:
            pdb_str = f.read()

        configs, systems = PDB_oxDNA(pdb_str)

        assert len(systems[0].strands) == len(system.strands), \
            "Should preserve strand count"


# =============================================================================
# Round-Trip Conversion Tests
# =============================================================================

class TestRoundTripConversion:
    """Tests for round-trip conversion (oxDNA -> PDB -> oxDNA)."""

    def test_roundtrip_correctness(self, system_and_conf, temp_output_dir):
        """
        End-to-end test: convert oxDNA to PDB and back, verifying correctness.

        Since the input topology is old-style (3'-5'), we use:
        - reverse=True in oxDNA_PDB to write PDB in standard 5'-3' order
        - old_top=True in PDB_oxDNA to convert back to 3'-5' order

        This test verifies:
        1. Sequences are preserved exactly
        2. Nucleotide count is preserved
        3. Positions are approximately preserved
        """
        system, conf = system_and_conf
        out_basename = str(temp_output_dir / "output")

        # Get original data
        original_sequences = [strand.get_sequence() for strand in system.strands]
        original_nucleotide_count = sum(len(strand) for strand in system.strands)
        original_com = np.mean(conf.positions, axis=0)

        # Convert to PDB with reverse=True since input is old-style (3'-5')
        oxDNA_PDB(conf, system, out_basename, reverse=True)
        pdb_file = temp_output_dir / "output.pdb"

        # Convert back to oxDNA with old_top=True to match original format
        with open(pdb_file) as f:
            pdb_str = f.read()

        configs, systems = PDB_oxDNA(pdb_str, old_top=True)
        converted_system = systems[0]
        converted_conf = configs[0]

        # 1. Verify sequences match exactly
        converted_sequences = [strand.get_sequence() for strand in converted_system.strands]

        assert len(original_sequences) == len(converted_sequences), \
            "Should have same number of strands"

        assert all([original_sequences[i] == converted_sequences[i] for i in range(len(original_sequences))]), \
            f"Sequences should match exactly.\nOriginal: {original_sequences}\nConverted: {converted_sequences}"

        # 2. Verify nucleotide count
        converted_nucleotide_count = sum(len(strand) for strand in converted_system.strands)
        assert original_nucleotide_count == converted_nucleotide_count, \
            f"Nucleotide count mismatch. Original: {original_nucleotide_count}, Converted: {converted_nucleotide_count}"

        # 3. Verify positions are approximately preserved
        # (won't be exact due to backmapping approximations)
        converted_com = np.mean(converted_conf.positions, axis=0)
        distance = np.linalg.norm(original_com - converted_com)
        assert distance < 1.0, \
            f"Center of mass should be preserved within 1 oxDNA unit, got {distance}"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestOxDNAPDBCLIParser:
    """Tests for the oxDNA_PDB CLI argument parser."""

    def test_parser_requires_positional_args(self):
        """Test that parser requires topology, configuration, and direction."""
        parser = oxDNA_PDB_cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_accepts_all_options(self):
        """Test parser accepts all options."""
        parser = oxDNA_PDB_cli_parser()
        args = parser.parse_args([
            "topology.top",
            "config.dat",
            "35",
            "-o", "output.pdb",
            "-d", "53",
            "-H",
            "-u",
            "-1",
            "-q"
        ])

        assert args.topology == "topology.top"
        assert args.configuration == "config.dat"
        assert args.direction == "35"
        assert args.output == "output.pdb"
        assert args.output_direction == "53"
        assert args.hydrogen is True
        assert args.uniform_residue_names is True
        assert args.one_file_per_strand is True
        assert args.quiet is True


class TestPDBOxDNACLIParser:
    """Tests for the PDB_oxDNA CLI argument parser."""

    def test_parser_requires_pdb_file(self):
        """Test that parser requires PDB file argument."""
        parser = PDB_oxDNA_cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_accepts_all_options(self):
        """Test parser accepts all options."""
        parser = PDB_oxDNA_cli_parser()
        args = parser.parse_args([
            "input.pdb",
            "-o", "output",
            "-b",
            "-q"
        ])

        assert args.pdb_file == "input.pdb"
        assert args.output == "output"
        assert args.backward is True
        assert args.quiet is True

    def test_parser_defaults(self):
        """Test parser default values."""
        parser = PDB_oxDNA_cli_parser()
        args = parser.parse_args(["input.pdb"])

        assert args.pdb_file == "input.pdb"
        assert args.output is None
        assert args.backward is False
        assert args.quiet is False


# =============================================================================
# Hybrid36 Encoding Tests
# =============================================================================

# Single-nucleotide PDB with hybrid36 atom serials (serials 164067-164100).
# Used to verify that PDB_oxDNA correctly parses atom serials beyond 99999.
HYBRID36_PDB = """\
ATOM  A1DFN  N9  DG3 b  17     579.947 852.831 955.384  1.00  1.00
ATOM  A1DFO  C8  DG3 b  17     578.625 852.499 955.308  1.00  1.00
ATOM  A1DFP  N7  DG3 b  17     578.362 851.268 954.982  1.00  1.00
ATOM  A1DFQ  C5  DG3 b  17     579.637 850.726 954.823  1.00  1.00
ATOM  A1DFR  C6  DG3 b  17     580.023 849.411 954.468  1.00  1.00
ATOM  A1DFS  O6  DG3 b  17     579.316 848.438 954.225  1.00  1.00
ATOM  A1DFT  N1  DG3 b  17     581.398 849.278 954.418  1.00  1.00
ATOM  A1DFU  C2  DG3 b  17     582.306 850.279 954.677  1.00  1.00
ATOM  A1DFV  N2  DG3 b  17     583.687 849.936 954.574  1.00  1.00
ATOM  A1DFW  N3  DG3 b  17     581.951 851.521 955.012  1.00  1.00
ATOM  A1DFX  C4  DG3 b  17     580.606 851.671 955.067  1.00  1.00
ATOM  A1DFY  H8  DG3 b  17     577.836 853.221 955.510  1.00  1.00
ATOM  A1DFZ  H1  DG3 b  17     581.718 848.353 954.167  1.00  1.00
ATOM  A1DG0 H21  DG3 b  17     583.942 849.000 954.294  1.00  1.00
ATOM  A1DG1 H22  DG3 b  17     584.392 850.629 954.782  1.00  1.00
ATOM  A1DG2  P   DG3 b  17     576.205 856.057 957.236  1.00  1.00
ATOM  A1DG3 OP1  DG3 b  17     576.173 857.443 957.758  1.00  1.00
ATOM  A1DG4 OP2  DG3 b  17     575.448 855.724 956.007  1.00  1.00
ATOM  A1DG5 O5'  DG3 b  17     577.745 855.651 956.980  1.00  1.00
ATOM  A1DG6 C5'  DG3 b  17     578.792 856.097 957.857  1.00  1.00
ATOM  A1DG7 C4'  DG3 b  17     580.190 855.788 957.310  1.00  1.00
ATOM  A1DG8 O4'  DG3 b  17     580.396 854.374 957.137  1.00  1.00
ATOM  A1DG9 C3'  DG3 b  17     580.471 856.419 955.940  1.00  1.00
ATOM  A1DGA O3'  DG3 b  17     581.871 856.664 955.759  1.00  1.00
ATOM  A1DGB C2'  DG3 b  17     579.956 855.319 955.032  1.00  1.00
ATOM  A1DGC C1'  DG3 b  17     580.582 854.116 955.736  1.00  1.00
ATOM  A1DGD H5'  DG3 b  17     578.682 855.603 958.823  1.00  1.00
ATOM  A1DGE H5'' DG3 b  17     578.703 857.173 958.001  1.00  1.00
ATOM  A1DGF H4'  DG3 b  17     580.926 856.158 958.024  1.00  1.00
ATOM  A1DGG H3'  DG3 b  17     579.938 857.363 955.823  1.00  1.00
ATOM  A1DGH HO3' DG3 b  17     582.162 857.194 956.504  1.00  1.00
ATOM  A1DGI H2'  DG3 b  17     578.873 855.260 955.136  1.00  1.00
ATOM  A1DGJ H2'' DG3 b  17     580.255 855.403 953.987  1.00  1.00
ATOM  A1DGK H1'  DG3 b  17     581.647 854.086 955.508  1.00  1.00
TER
END
"""


class TestHybrid36:
    """Tests for hybrid36 atom serial encoding/decoding beyond 99999."""

    def test_format_and_parse_roundtrip(self):
        """
        End-to-end encoding round-trip: verify _format_atom_serial and
        _parse_atom_serial are exact inverses across the decimal/hybrid36
        boundary and within the hybrid36 range.

        Key values exercised:
          99999  -> "99999"  (last pure-decimal serial)
          100000 -> "A0000"  (first hybrid36 serial)
          100009 -> "A0009"  (last digit still numeric)
          100010 -> "A000A"  (first carry into base-36 territory)
          100035 -> "A000Z"  (Z = 35, last value in units digit)
          100036 -> "A0010"  (carry into tens digit)
          164067 -> "A1DFN"  (first atom serial in HYBRID36_PDB fixture)
          164100 -> "A1DGK"  (last atom serial in HYBRID36_PDB fixture)
        """
        cases = [
            (99999,  "99999"),
            (100000, "A0000"),
            (100009, "A0009"),
            (100010, "A000A"),
            (100035, "A000Z"),
            (100036, "A0010"),
            (164067, "A1DFN"),
            (164100, "A1DGK"),
        ]
        for value, encoded in cases:
            assert _format_atom_serial(value).strip() == encoded, \
                f"_format_atom_serial({value}) should be '{encoded}'"
            assert _parse_atom_serial(encoded) == value, \
                f"_parse_atom_serial('{encoded}') should be {value}"

    def test_pdb_oxdna_hybrid36_single_nucleotide(self):
        """
        End-to-end conversion of a PDB with hybrid36 atom serials (>100000).

        The fixture is a single DG3 (DNA guanine, 3' terminal) residue whose
        atom serials start at A1DFN (164067) and span into base-36 territory
        (e.g. A1DGA, A1DGB ...).  This exercises both the letter-prefix and the
        all-base-36-digits-past-9 cases (A000A style) within one parse call.

        Checks (in increasing specificity):
          1. No exception is raised (hybrid36 serials are parsed without ValueError)
          2. Returns exactly one configuration and one system
          3. The system contains exactly one strand with exactly one nucleotide
          4. The nucleotide base type is 'G'
          5. The nucleotide position is within 5 oxDNA units of the expected
             centre-of-mass computed from the all-atom coordinates
          6. The orientation vectors a1 and a3 are unit vectors
        """
        configs, systems = PDB_oxDNA(HYBRID36_PDB)

        # 1+2: structure
        assert len(configs) == 1 and len(systems) == 1

        system = systems[0]
        conf   = configs[0]

        # 3: one strand, one nucleotide
        assert len(system.strands) == 1
        strand = system.strands[0]
        assert len(strand.monomers) == 1

        monomer = strand.monomers[0]

        # 4: correct base identity
        assert monomer.btype == 'G', \
            f"Expected base 'G', got '{monomer.btype}'"

        # 5: position close to the all-atom COM (angstrom -> oxDNA units)
        ANGSTROM_TO_OXDNA = 1.0 / 8.518
        atom_coords = np.array([
            [579.947, 852.831, 955.384], [578.625, 852.499, 955.308],
            [578.362, 851.268, 954.982], [579.637, 850.726, 954.823],
            [580.023, 849.411, 954.468], [579.316, 848.438, 954.225],
            [581.398, 849.278, 954.418], [582.306, 850.279, 954.677],
            [583.687, 849.936, 954.574], [581.951, 851.521, 955.012],
            [580.606, 851.671, 955.067], [577.836, 853.221, 955.510],
            [581.718, 848.353, 954.167], [583.942, 849.000, 954.294],
            [584.392, 850.629, 954.782], [576.205, 856.057, 957.236],
            [576.173, 857.443, 957.758], [575.448, 855.724, 956.007],
            [577.745, 855.651, 956.980], [578.792, 856.097, 957.857],
            [580.190, 855.788, 957.310], [580.396, 854.374, 957.137],
            [580.471, 856.419, 955.940], [581.871, 856.664, 955.759],
            [579.956, 855.319, 955.032], [580.582, 854.116, 955.736],
            [578.682, 855.603, 958.823], [578.703, 857.173, 958.001],
            [580.926, 856.158, 958.024], [579.938, 857.363, 955.823],
            [582.162, 857.194, 956.504], [578.873, 855.260, 955.136],
            [580.255, 855.403, 953.987], [581.647, 854.086, 955.508],
        ])
        expected_pos = np.mean(atom_coords, axis=0) * ANGSTROM_TO_OXDNA
        # conf.positions is a flat array ordered by system position (not by monomer.id)
        actual_pos = conf.positions[0]
        assert np.linalg.norm(actual_pos - expected_pos) < 5.0, \
            f"Position {actual_pos} too far from expected COM {expected_pos}"

        # 6: orientation vectors are unit vectors
        a1 = conf.a1s[0]
        a3 = conf.a3s[0]
        assert abs(np.linalg.norm(a1) - 1.0) < 1e-6, "a1 should be a unit vector"
        assert abs(np.linalg.norm(a3) - 1.0) < 1e-6, "a3 should be a unit vector"
