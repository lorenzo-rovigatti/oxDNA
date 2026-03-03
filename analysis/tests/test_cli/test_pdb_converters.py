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

from oxDNA_analysis_tools.oxDNA_PDB import oxDNA_PDB, cli_parser as oxDNA_PDB_cli_parser
from oxDNA_analysis_tools.PDB_oxDNA import PDB_oxDNA, cli_parser as PDB_oxDNA_cli_parser
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
