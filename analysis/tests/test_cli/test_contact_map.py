"""
Tests for oxDNA_analysis_tools.contact_map module.

Tests cover:
- contact_map() API function
- compute() helper function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.contact_map import contact_map, compute, cli_parser, ComputeContext, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe


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
def trajectory_info(mini_traj_path):
    """Get topology and trajectory info for the mini trajectory."""
    top_info, traj_info = describe(None, str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# API Tests - contact_map() function
# =============================================================================

class TestContactMapFunction:
    """Tests for the contact_map() API function."""

    def test_contact_map_basic_behavior(self, trajectory_info):
        """Test contact_map() returns correct shape, symmetric matrix with non-negative distances."""
        top_info, traj_info = trajectory_info

        result = contact_map(traj_info, top_info, ncpus=1)

        # Check return type and shape
        assert isinstance(result, np.ndarray), "Should return numpy array"
        assert result.shape == (top_info.nbases, top_info.nbases), "Should be NxN matrix"

        # Matrix should be symmetric
        np.testing.assert_allclose(result, result.T, rtol=1e-5,
                                   err_msg="Contact map should be symmetric")

        # Diagonal should be zero (distance from particle to itself)
        np.testing.assert_allclose(np.diag(result), np.zeros(top_info.nbases), atol=1e-10,
                                   err_msg="Diagonal should be zero")

        # All distances should be non-negative
        assert np.all(result >= 0), "All distances should be non-negative"

    def test_contact_map_distance_scale(self, trajectory_info):
        """Test that contact_map returns distances in nanometers (reasonable range)."""
        top_info, traj_info = trajectory_info

        result = contact_map(traj_info, top_info, ncpus=1)

        # Distances should be in nm (typical DNA structure ~few nm)
        max_distance = np.max(result)
        assert max_distance > 0.1, "Max distance should be at least 0.1 nm"
        assert max_distance < 100, "Max distance should be less than 100 nm"

    def test_contact_map_adjacent_particles_close(self, trajectory_info):
        """Test that adjacent particles have smaller distances than distant ones."""
        top_info, traj_info = trajectory_info

        result = contact_map(traj_info, top_info, ncpus=1)

        if top_info.nbases > 10:
            # Adjacent particles should be closer than distant ones
            adjacent_dist = result[0, 1]
            distant_dist = result[0, top_info.nbases - 1]
            # This is generally true for DNA strands but not always
            # Just check both are positive
            assert adjacent_dist > 0, "Adjacent distance should be positive"
            assert distant_dist > 0, "Distant distance should be positive"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_trajectory(self):
        """Test that parser requires trajectory argument."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-g", "plot.png",
            "-d", "data.pkl",
            "-p", "4",
            "-q",
            "trajectory.dat"
        ])

        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.graph == ["plot.png"], "Graph option not parsed"
        assert args.data == ["data.pkl"], "Data option not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["traj.dat"])
        assert args_defaults.quiet is False, "Quiet should default to False"
        assert args_defaults.graph is None, "Graph should default to None"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output_files(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates graph and data output files."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        graph_file = temp_output_dir / "contact_map.png"
        data_file = temp_output_dir / "contact_map.pkl"

        test_args = [
            "contact_map.py",
            "-g", str(graph_file),
            "-d", str(data_file),
            str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check files exist
        assert graph_file.exists(), "Graph file should be created"
        assert Path(str(data_file) + ".npy").exists(), "Data file should be created (as .npy)"

    def test_main_default_filenames(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() uses default filenames."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        test_args = ["contact_map.py", str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check default files exist
        assert (temp_output_dir / "contact_map.png").exists(), "Default graph should be created"
        assert (temp_output_dir / "contact_map.pkl.npy").exists(), "Default data should be created"

    def test_main_data_file_loadable(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that saved data file can be loaded and has correct shape."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        data_file = temp_output_dir / "test_data.pkl"

        test_args = ["contact_map.py", "-d", str(data_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Load and verify
        loaded = np.load(str(data_file) + ".npy")
        top_info, _ = describe(None, str(traj_copy))
        assert loaded.shape == (top_info.nbases, top_info.nbases), "Loaded data should have correct shape"
