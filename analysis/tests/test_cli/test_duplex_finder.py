"""
Tests for oxDNA_analysis_tools.duplex_finder module.

Tests cover:
- Duplex dataclass
- find_duplex() function
- duplex_finder() API function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.duplex_finder import (
    Duplex,
    find_duplex,
    duplex_finder,
    cli_parser,
    main
)
from oxDNA_analysis_tools.UTILS.RyeReader import describe, strand_describe
from oxDNA_analysis_tools.UTILS.data_structures import Monomer


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def test_resources():
    """Get the path to test resources directory."""
    return Path(__file__).parent.parent / "resources"


@pytest.fixture(scope="module")
def mini_traj_path(test_resources):
    """Path to the mini trajectory file."""
    return test_resources / "minitraj.dat"


@pytest.fixture(scope="module")
def topology_path(test_resources):
    """Path to the topology file."""
    return test_resources / "rna_tile.top"


@pytest.fixture(scope="module")
def input_file_path(test_resources):
    """Path to the input file."""
    return test_resources / "input_rna"


@pytest.fixture(scope="module")
def trajectory_info(topology_path, mini_traj_path):
    """Get topology and trajectory info."""
    top_info, traj_info = describe(str(topology_path), str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture(scope="module")
def system_and_monomers(topology_path):
    """Get system and monomers from topology."""
    system, monomers = strand_describe(str(topology_path))
    return system, monomers


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


@pytest.fixture
def simple_monomers():
    """Create a simple list of monomers for testing find_duplex."""
    # Create 8 monomers: 4 in one strand (0-3), 4 in another (4-7)
    monomers = []
    for i in range(8):
        m = Monomer(id=i, btype='A')
        m.strand = 0 if i < 4 else 1
        m.n3 = i - 1 if (i > 0 and i != 4) else None
        m.n5 = i + 1 if (i < 3 or (i >= 4 and i < 7)) else None
        m.pair = None
        monomers.append(m)
    return monomers


# =============================================================================
# Unit Tests - Duplex dataclass
# =============================================================================

class TestDuplexDataclass:
    """Tests for the Duplex dataclass."""

    def test_duplex_behavior(self):
        """Test Duplex creation, field access, and mutability."""
        d = Duplex(
            time=0,
            index=1,
            start1=0,
            end1=10,
            start2=20,
            end2=30,
            axis=np.array([1, 0, 0]),
            pos=np.array([0, 0, 0])
        )

        # Check fields
        assert d.time == 0
        assert d.index == 1
        assert d.start1 == 0
        assert d.end1 == 10
        assert d.start2 == 20
        assert d.end2 == 30
        np.testing.assert_array_equal(d.axis, [1, 0, 0])
        np.testing.assert_array_equal(d.pos, [0, 0, 0])

        # Check mutability
        d.end1 = 15
        d.axis = np.array([0, 1, 0])
        assert d.end1 == 15
        np.testing.assert_array_equal(d.axis, [0, 1, 0])


# =============================================================================
# Unit Tests - find_duplex()
# =============================================================================

class TestFindDuplex:
    """Tests for the find_duplex() function."""

    def test_find_duplex_behavior(self, simple_monomers):
        """Test find_duplex returns list, handles no pairs, and edge cases."""
        # Test returns list with unpaired monomers
        result = find_duplex(simple_monomers)
        assert isinstance(result, list), "Should return list"
        assert len(result) == 0, "No duplexes without pairs"

        # Test with some pairs (3 nucleotides - below minimum)
        simple_monomers[0].pair = 7
        simple_monomers[7].pair = 0
        simple_monomers[1].pair = 6
        simple_monomers[6].pair = 1
        simple_monomers[2].pair = 5
        simple_monomers[5].pair = 2

        result_with_pairs = find_duplex(simple_monomers)
        assert isinstance(result_with_pairs, list), "Should return list"

        # Test edge cases
        assert find_duplex([]) == [], "Empty input should return empty list"

        single_m = Monomer(id=0, btype='A')
        single_m.strand = 0
        single_m.n3 = None
        single_m.n5 = None
        single_m.pair = None
        assert find_duplex([single_m]) == [], "Single monomer should return empty list"


# =============================================================================
# API Tests - duplex_finder() function
# =============================================================================

class TestDuplexFinderFunction:
    """Tests for the duplex_finder() API function."""

    def test_duplex_finder_behavior(self, trajectory_info, system_and_monomers, input_file_path, test_resources):
        """Test duplex_finder() returns correct structure with one list per configuration."""
        top_info, traj_info = trajectory_info
        system, monomers = system_and_monomers

        import os
        original_dir = os.getcwd()
        try:
            os.chdir(test_resources)
            result = duplex_finder(traj_info, top_info, str(input_file_path), monomers, ncpus=1)
        finally:
            os.chdir(original_dir)

        # Check return type
        assert isinstance(result, list), "Should return list"
        assert len(result) == traj_info.nconfs, "Should have one entry per configuration"

        # Each element should be a list of duplexes for that timestep
        for step_duplexes in result:
            assert isinstance(step_duplexes, list), "Each step should have list of duplexes"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_behavior(self):
        """Test parser requires arguments, accepts all options, and has correct defaults."""
        parser = cli_parser()

        # Requires input and trajectory
        with pytest.raises(SystemExit):
            parser.parse_args([])
        with pytest.raises(SystemExit):
            parser.parse_args(["input.txt"])

        # All options
        args = parser.parse_args([
            "-p", "4",
            "-o", "angles.txt",
            "-q",
            "input.txt",
            "trajectory.dat"
        ])

        assert args.input == ["input.txt"], "Input not parsed"
        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.output == ["angles.txt"], "Output option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Defaults
        args_defaults = parser.parse_args(["input.txt", "traj.dat"])
        assert args_defaults.parallel is None
        assert args_defaults.output is None
        assert args_defaults.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output(self, test_resources, temp_output_dir, monkeypatch):
        """Test main() creates output file with correct format."""
        monkeypatch.chdir(test_resources)

        output_file = temp_output_dir / "angles.txt"

        test_args = [
            "duplex_finder.py",
            "-o", str(output_file),
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output file should be created"

        # Check header format
        with open(output_file) as f:
            header = f.readline()

        expected_columns = ["time", "duplex", "start1", "end1", "start2", "end2",
                          "axisX", "axisY", "axisZ", "hel_pos"]
        for col in expected_columns:
            assert col in header, f"Header should contain '{col}'"

    def test_main_default_output(self, test_resources, monkeypatch):
        """Test main() uses default output filename."""
        monkeypatch.chdir(test_resources)

        test_args = ["duplex_finder.py", "input_rna", "minitraj.dat"]

        with patch.object(sys, 'argv', test_args):
            main()

        default_output = test_resources / "angles.txt"
        assert default_output.exists(), "Default output 'angles.txt' should be created"
        # Clean up
        default_output.unlink()
