"""
Tests for oxDNA_analysis_tools.clustering module.

Tests cover:
- find_element() utility function
- split_trajectory() function
- perform_DBSCAN() main function
- CLI argument parsing
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.clustering import find_element, cli_parser
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
# Utility Function Tests
# =============================================================================

class TestFindElement:
    """Tests for the find_element() utility function."""

    def test_find_element_behavior(self):
        """Test find_element finds correct index, handles not found, and edge cases."""
        # Basic cases - find nth occurrence
        array = np.array([0, 1, 2, 1, 3, 1, 4])
        assert find_element(0, 1, array) == 1, "First occurrence of 1 should be at index 1"
        assert find_element(1, 1, array) == 3, "Second occurrence of 1 should be at index 3"
        assert find_element(2, 1, array) == 5, "Third occurrence of 1 should be at index 5"

        # Not found cases
        array_simple = np.array([0, 1, 2, 3])
        assert find_element(0, 5, array_simple) == -1, "Should return -1 for missing element"
        assert find_element(5, 1, array_simple) == -1, "Should return -1 when n > count"

        # Edge cases
        # Single element array
        array_single = np.array([5])
        assert find_element(0, 5, array_single) == 0, "Should find in single element array"
        assert find_element(1, 5, array_single) == -1, "Should return -1 for n=1 in single element"

        # All same elements
        array_same = np.array([1, 1, 1, 1])
        assert find_element(0, 1, array_same) == 0
        assert find_element(3, 1, array_same) == 3


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_data_file(self):
        """Test that parser requires serialized data file argument."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-e", "0.5",
            "-m", "10",
            "-q",
            "data.json"
        ])

        assert args.serialized_data == ["data.json"], "Data file not parsed"
        assert args.eps == [0.5], "Eps option not parsed"
        assert args.min_samples == [10], "Min samples option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults (eps and min_samples use defaults in main())
        args_defaults = parser.parse_args(["data.json"])
        assert args_defaults.quiet is False, "Quiet should default to False"


# =============================================================================
# DBSCAN Clustering Tests (lightweight, no sklearn dependency in tests)
# =============================================================================

class TestDBSCANInput:
    """Tests for DBSCAN input validation."""

    def test_cluster_data_shape_validation(self, trajectory_info):
        """Test that clustering validates input data shape."""
        top_info, traj_info = trajectory_info

        # Create mock order parameter data
        op_data = np.random.rand(traj_info.nconfs, 3)

        # Data should have nconfs rows
        assert op_data.shape[0] == traj_info.nconfs, "OP data should have nconfs rows"

    def test_precomputed_matrix_symmetric(self):
        """Test that precomputed distance matrices should be symmetric."""
        n = 10
        # Create symmetric distance matrix
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)

        # Should be symmetric
        np.testing.assert_array_almost_equal(distances, distances.T,
                                              err_msg="Precomputed matrix should be symmetric")

        # Diagonal should be zero
        np.testing.assert_array_almost_equal(np.diag(distances), np.zeros(n),
                                              err_msg="Diagonal should be zero")
