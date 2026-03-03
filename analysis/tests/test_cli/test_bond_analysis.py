"""
Tests for oxDNA_analysis_tools.bond_analysis module.

Tests cover:
- bond_analysis() API function
- oxView_overlay() output function
- plot_trajectories() output function
- CLI argument parsing
- main() CLI entry point

Note: These tests require oxpy and a valid input file.
"""
import json
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.UTILS.RyeReader import describe

# Skip all tests if oxpy is not available
pytest.importorskip("oxpy")

from oxDNA_analysis_tools.bond_analysis import (
    bond_analysis,
    compute,
    ComputeContext,
    oxView_overlay,
    plot_trajectories,
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
def designed_pairs_path(test_resources):
    """Path to the designed pairs file."""
    return test_resources / "designed_pairs.txt"


@pytest.fixture(scope="module")
def trajectory_info(mini_traj_path, test_resources):
    """Get topology and trajectory info for the mini trajectory."""
    top_file = test_resources / "rna_tile.top"
    top_info, traj_info = describe(str(top_file), str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture(scope="module")
def designed_pairs(designed_pairs_path):
    """Load designed pairs as dict."""
    with open(designed_pairs_path, 'r') as f:
        pairs_txt = f.readlines()
    return {int(p[0]): int(p[1]) for p in [p.split() for p in pairs_txt]}


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# API Tests - bond_analysis() function
# =============================================================================

class TestBondAnalysisFunction:
    """Tests for the bond_analysis() API function."""

    def test_bond_analysis_returns_correct_shapes(
        self, trajectory_info, designed_pairs, input_file_path, test_resources, monkeypatch
    ):
        """Test bond_analysis() returns arrays with correct shapes and non-negative values."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info

        total_bonds, correct_bonds, incorrect_bonds, nt_array = bond_analysis(
            traj_info, top_info, designed_pairs, str(input_file_path), ncpus=1
        )

        # Check shapes
        assert total_bonds.shape == (traj_info.nconfs,), "total_bonds should have nconfs elements"
        assert correct_bonds.shape == (traj_info.nconfs,), "correct_bonds should have nconfs elements"
        assert incorrect_bonds.shape == (traj_info.nconfs,), "incorrect_bonds should have nconfs elements"
        assert nt_array.shape == (top_info.nbases,), "nt_array should have nbases elements"

        # Check non-negative values
        assert np.all(total_bonds >= 0), "total_bonds should be non-negative"
        assert np.all(correct_bonds >= 0), "correct_bonds should be non-negative"
        assert np.all(incorrect_bonds >= 0), "incorrect_bonds should be non-negative"
        assert np.all(nt_array >= 0), "nt_array should be non-negative"

        # Total should be sum of correct and incorrect
        np.testing.assert_array_equal(
            total_bonds, correct_bonds + incorrect_bonds,
            err_msg="Total bonds should equal correct + incorrect"
        )

        # Directly test the compute function
        ctx = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            designed_pairs=designed_pairs,
            input_file=str(input_file_path)
        )

        tot, corr, incorr, nt_arr = compute(ctx, chunk_size=1, chunk_id=0)

        # Check return types
        assert isinstance(tot, np.ndarray), "compute() should return numpy arrays"
        assert isinstance(corr, np.ndarray), "correct_bonds should be numpy array"
        assert isinstance(incorr, np.ndarray), "incorrect_bonds should be numpy array"
        assert isinstance(nt_arr, np.ndarray), "nt_array should be numpy array"

        # Check shapes (single chunk)
        assert len(tot) == 1, "Single chunk should return 1 value"
        assert len(corr) == 1, "Single chunk should return 1 value"
        assert len(incorr) == 1, "Single chunk should return 1 value"
        assert len(nt_arr) == top_info.nbases, "nt_array should have nbases elements"

        # Check values are non-negative integers
        assert tot[0] >= 0, "Total bonds should be non-negative"
        assert corr[0] >= 0, "Correct bonds should be non-negative"
        assert incorr[0] >= 0, "Incorrect bonds should be non-negative"
        assert np.all(nt_arr >= 0), "nt_array values should be non-negative"

        # Check that total = correct + incorrect (core logic validation)
        assert tot[0] == corr[0] + incorr[0], "Total bonds should equal correct + incorrect"

        # Check nt_arr only has values where designed pairs exist
        # For each designed pair (a, b), if there's a correct bond, both a and b should be incremented
        assert nt_arr.dtype == np.int64 or nt_arr.dtype == np.int32, "nt_array should be integer type"

        # Test a second chunk to ensure consistency
        if traj_info.nconfs > 1:
            tot2, corr2, incorr2, nt_arr2 = compute(ctx, chunk_size=1, chunk_id=1)
            assert len(tot2) == 1, "Second chunk should also return 1 value"
            assert tot2[0] == corr2[0] + incorr2[0], "Logic should be consistent across chunks"


# =============================================================================
# Output Function Tests
# =============================================================================

class TestOutputFunctions:
    """Tests for the output helper functions."""

    def test_oxView_overlay_creates_file(self, temp_output_dir, trajectory_info):
        """Test oxView_overlay creates valid JSON file."""
        top_info, _ = trajectory_info
        nt_array = np.random.rand(top_info.nbases)

        outfile = str(temp_output_dir / "bonds.json")
        oxView_overlay(nt_array, outfile)

        assert Path(outfile).exists(), "Output file should be created"

        # Verify it's valid JSON
        with open(outfile) as f:
            data = json.load(f)
        assert "occupancy" in data, "JSON should have 'occupancy' key"
        assert len(data["occupancy"]) == top_info.nbases, "Should have nbases values"

    def test_plot_trajectories_creates_file(self, temp_output_dir):
        """Test plot_trajectories creates PNG file."""
        correct = np.array([5, 6, 5, 7, 6])
        incorrect = np.array([1, 0, 1, 0, 1])
        designed = 10
        plotname = str(temp_output_dir / "bonds.png")

        plot_trajectories(correct, incorrect, designed, plotname)

        assert Path(plotname).exists(), "Plot file should be created"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_arguments(self):
        """Test that parser requires inputfile, trajectory, and designed_pairs arguments."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-o", "output.json",
            "-t", "plot.png",
            "-d", "data.json",
            "-p", "4",
            "-q",
            "input.txt",
            "trajectory.dat",
            "pairs.txt"
        ])

        assert args.inputfile == "input.txt", "Inputfile not parsed"
        assert args.trajectory == "trajectory.dat", "Trajectory not parsed"
        assert args.designed_pairs == "pairs.txt", "Designed pairs not parsed"
        assert args.outfile == "output.json", "Output option not parsed"
        assert args.traj_plot == "plot.png", "Plot option not parsed"
        assert args.data_file == "data.json", "Data option not parsed"
        assert args.parallel == 4, "Parallel option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["input.txt", "traj.dat", "pairs.txt"])
        assert args_defaults.quiet is False, "Quiet should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output_files(
        self, mini_traj_path, input_file_path, designed_pairs_path, test_resources, temp_output_dir, monkeypatch
    ):
        """Test main() creates all output files."""
        monkeypatch.chdir(test_resources)

        output_file = temp_output_dir / "test_bonds.json"
        plot_file = temp_output_dir / "test_bonds.png"

        test_args = [
            "bond_analysis.py",
            "-o", str(output_file),
            "-t", str(plot_file),
            str(input_file_path),
            str(mini_traj_path),
            str(designed_pairs_path)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output JSON should be created"
        assert plot_file.exists(), "Plot file should be created"
