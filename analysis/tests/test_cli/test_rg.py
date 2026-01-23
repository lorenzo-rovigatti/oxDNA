"""
Tests for oxDNA_analysis_tools.rg module (radius of gyration).

Tests cover:
- rg() API function
- compute() helper function
- CLI argument parsing
- main() CLI entry point
"""
import json
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.rg import rg, compute, cli_parser, ComputeContext, main
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
# Unit Tests - compute()
# =============================================================================

class TestCompute:
    """Tests for the compute() helper function."""

    def test_compute_behavior(self, trajectory_info):
        """Test compute() returns correct array with positive Rg values."""
        top_info, traj_info = trajectory_info

        ctx = ComputeContext(
            top_info=top_info,
            traj_info=traj_info
        )

        # Process first chunk
        result = compute(ctx, chunk_size=1, chunk_id=0)

        # Check return type and shape
        assert isinstance(result, np.ndarray), "compute() should return numpy array"
        assert len(result) == 1, "Single chunk should return one Rg value"

        # Rg should be positive
        assert result[0] > 0, "Radius of gyration should be positive"

        # Process second chunk
        result2 = compute(ctx, chunk_size=1, chunk_id=1)
        assert isinstance(result2, np.ndarray), "Second chunk should return array"
        assert result2[0] > 0, "Second Rg should be positive"


# =============================================================================
# API Tests - rg() function
# =============================================================================

class TestRgFunction:
    """Tests for the rg() API function."""

    def test_rg_behavior(self, trajectory_info):
        """Test rg() returns correct shape, positive values, reasonable range, and stability."""
        top_info, traj_info = trajectory_info

        result = rg(top_info, traj_info, ncpus=1)

        # Check return type and shape
        assert isinstance(result, np.ndarray), "Should return numpy array"
        assert result.shape == (traj_info.nconfs,), "Should have one Rg per configuration"

        # All Rg values should be positive
        assert np.all(result > 0), "Radius of gyration should be positive"

        # Rg should be in reasonable range for DNA (typically ~1-50 nm in oxDNA units)
        mean_rg = np.mean(result)
        assert mean_rg > 0.1, "Mean Rg should be at least 0.1"
        assert mean_rg < 1000, "Mean Rg should be less than 1000"

        # Standard deviation should be much smaller than mean (stable structure)
        std_rg = np.std(result)
        assert std_rg < mean_rg, "Rg variation should be less than mean (stable structure)"


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
            "-p", "4",
            "-o", "output.json",
            "-q",
            "trajectory.dat"
        ])

        assert args.trajectory == "trajectory.dat", "Trajectory not parsed"
        assert args.parallel == 4, "Parallel option not parsed"
        assert args.output == "output.json", "Output option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["traj.dat"])
        assert args_defaults.quiet is False, "Quiet should default to False"
        assert args_defaults.parallel is None, "Parallel should default to None"
        assert args_defaults.output is None, "Output should default to None"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_behavior(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates output file with correct structure, default naming, and positive values."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Test with explicit output file
        output_file = temp_output_dir / "test_rg.json"
        test_args = ["rg.py", "-o", str(output_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output file should be created"

        # Verify JSON structure
        with open(output_file) as f:
            data = json.load(f)

        assert "rg" in data, "JSON should have 'rg' key"

        top_info, traj_info = describe(None, str(traj_copy))
        assert len(data["rg"]) == traj_info.nconfs, "Should have one Rg per configuration"

        # All values should be positive
        assert all(v > 0 for v in data["rg"]), "All Rg values should be positive"

        # Test default output filename
        test_args_default = ["rg.py", str(traj_copy)]
        with patch.object(sys, 'argv', test_args_default):
            main()

        default_output = temp_output_dir / "rg.json"
        assert default_output.exists(), "Default output 'rg.json' should be created"
