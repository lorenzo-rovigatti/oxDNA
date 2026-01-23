"""
Tests for oxDNA_analysis_tools.output_bonds module.

Tests cover:
- output_bonds() API function
- parse_header() helper function
- get_potentials() helper function
- CLI argument parsing
- main() CLI entry point

Note: These tests require oxpy and a valid input file.
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.UTILS.RyeReader import describe

# Skip all tests if oxpy is not available
pytest.importorskip("oxpy")

from oxDNA_analysis_tools.output_bonds import (
    output_bonds,
    parse_header,
    cli_parser,
    main
)


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
def input_file_path(test_resources):
    """Path to the input file."""
    return test_resources / "input_rna"


@pytest.fixture(scope="module")
def trajectory_info(mini_traj_path, test_resources):
    """Get topology and trajectory info for the mini trajectory."""
    top_file = test_resources / "rna_tile.top"
    top_info, traj_info = describe(str(top_file), str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# Helper Function Tests
# =============================================================================

class TestParseHeader:
    """Tests for the parse_header() helper function."""

    def test_parse_header_extracts_potentials(self):
        """Test parse_header extracts potential names from header string."""
        header = "# id1 id2 fene bexc stack nexc hb cr_stack cx_stack total, t = 1000000"
        result = parse_header(header)

        assert isinstance(result, list), "Should return list"
        assert "fene" in result, "Should contain 'fene'"
        assert "total" in result, "Should contain 'total'"
        # id1 and id2 should be excluded
        assert "id1" not in result, "Should not contain 'id1'"
        assert "id2" not in result, "Should not contain 'id2'"


# =============================================================================
# API Tests - output_bonds() function
# =============================================================================

class TestOutputBondsFunction:
    """Tests for the output_bonds() API function."""

    def test_output_bonds_visualize_mode(
        self, trajectory_info, input_file_path, test_resources, monkeypatch
    ):
        """Test output_bonds() with visualize=True returns energies and potential names."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info

        energies, pot_names = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=1, ncpus=1
        )

        # Check return types
        assert isinstance(energies, np.ndarray), "Energies should be numpy array"
        assert isinstance(pot_names, list), "Potential names should be list"

        # Check shapes
        assert energies.shape[0] == top_info.nbases, "Energies should have nbases rows"
        assert energies.shape[1] == len(pot_names), "Energies columns should match potential count"

        # Check potential names are strings
        assert all(isinstance(p, str) for p in pot_names), "Potential names should be strings"

    def test_output_bonds_conversion_factor(
        self, trajectory_info, input_file_path, test_resources, monkeypatch
    ):
        """Test output_bonds() respects conversion factor."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info

        # Get energies with conversion factor 1
        energies_1, _ = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=1, ncpus=1
        )

        # Get energies with conversion factor 2
        energies_2, _ = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=2, ncpus=1
        )

        # Energies with factor 2 should be ~2x those with factor 1
        # (may not be exact due to accumulation)
        ratio = np.mean(np.abs(energies_2)) / np.mean(np.abs(energies_1)) if np.mean(np.abs(energies_1)) > 0 else 1
        assert 1.5 < ratio < 2.5, f"Conversion factor should scale energies, got ratio {ratio}"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_arguments(self):
        """Test that parser requires inputfile and trajectory arguments."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-v", "output.json",
            "-p", "4",
            "-u", "pNnm",
            "-q",
            "input.txt",
            "trajectory.dat"
        ])

        assert args.inputfile == ["input.txt"], "Inputfile not parsed"
        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.outfile == ["output.json"], "Output option not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.units == ["pNnm"], "Units option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["input.txt", "traj.dat"])
        assert args_defaults.quiet is False, "Quiet should default to False"
        assert args_defaults.outfile is None, "Outfile should default to None"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_visualize_creates_files(
        self, mini_traj_path, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """Test main() with -v creates oxView overlay JSON files."""
        monkeypatch.chdir(test_resources)

        output_base = temp_output_dir / "energies.json"

        test_args = [
            "output_bonds.py",
            "-v", str(output_base),
            str(input_file_path),
            str(mini_traj_path)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Should create one file per potential type
        json_files = list(temp_output_dir.glob("energies_*.json"))
        assert len(json_files) > 0, "Should create at least one JSON file per potential"

    def test_main_with_units_option(
        self, mini_traj_path, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """Test main() with different unit options."""
        monkeypatch.chdir(test_resources)

        # Test with pNnm units
        output_pn = temp_output_dir / "energies_pn.json"
        test_args_pn = [
            "output_bonds.py",
            "-v", str(output_pn),
            "-u", "pNnm",
            str(input_file_path),
            str(mini_traj_path)
        ]

        with patch.object(sys, 'argv', test_args_pn):
            main()

        # Test with oxDNA units
        output_ox = temp_output_dir / "energies_ox.json"
        test_args_ox = [
            "output_bonds.py",
            "-v", str(output_ox),
            "-u", "oxDNA",
            str(input_file_path),
            str(mini_traj_path)
        ]

        with patch.object(sys, 'argv', test_args_ox):
            main()

        # Both should create files
        pn_files = list(temp_output_dir.glob("energies_pn_*.json"))
        ox_files = list(temp_output_dir.glob("energies_ox_*.json"))
        assert len(pn_files) > 0, "pNnm units should create files"
        assert len(ox_files) > 0, "oxDNA units should create files"
