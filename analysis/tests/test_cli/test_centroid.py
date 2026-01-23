"""
Tests for oxDNA_analysis_tools.centroid module.

Tests cover:
- centroid() API function
- compute_centroid() helper function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.centroid import centroid, compute_centroid, cli_parser, ComputeContext, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, inbox
from oxDNA_analysis_tools.UTILS.data_structures import Configuration


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
def mean_conf_path(test_resources):
    """Path to the mean configuration file."""
    return test_resources / "mean.dat"


@pytest.fixture(scope="module")
def trajectory_info(mini_traj_path):
    """Get topology and trajectory info for the mini trajectory."""
    top_info, traj_info = describe(None, str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture(scope="module")
def reference_conf(trajectory_info, mean_conf_path):
    """Get reference configuration for centroid calculation."""
    top_info, _ = trajectory_info
    _, mean_info = describe(None, str(mean_conf_path))
    return get_confs(top_info, mean_info, 0, 1)[0]


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# API Tests - centroid() function
# =============================================================================

class TestCentroidFunction:
    """Tests for the centroid() API function."""

    def test_centroid_basic_behavior(self, trajectory_info, reference_conf):
        """Test centroid() returns Configuration with correct shape and positive RMSD."""
        top_info, traj_info = trajectory_info

        result_conf, rmsd = centroid(traj_info, top_info, reference_conf, ncpus=1)

        # Check return types
        assert isinstance(result_conf, Configuration), "Should return Configuration"
        assert isinstance(rmsd, float), "RMSD should be float"

        # Check shapes
        assert result_conf.positions.shape == (top_info.nbases, 3), "Wrong positions shape"
        assert result_conf.a1s.shape == (top_info.nbases, 3), "Wrong a1s shape"
        assert result_conf.a3s.shape == (top_info.nbases, 3), "Wrong a3s shape"

        # RMSD should be non-negative and in reasonable range (nm)
        assert rmsd >= 0, "RMSD should be non-negative"
        assert rmsd < 100, "RMSD should be in reasonable range"

    def test_centroid_with_custom_indexes(self, trajectory_info, reference_conf):
        """Test centroid() with subset of particle indexes for alignment."""
        top_info, traj_info = trajectory_info

        # Use only first half of particles for alignment
        indexes = list(range(top_info.nbases // 2))

        result_conf, rmsd = centroid(traj_info, top_info, reference_conf, indexes=indexes, ncpus=1)

        # Should still return full configuration
        assert result_conf.positions.shape == (top_info.nbases, 3), "Should return full system"
        assert rmsd >= 0, "RMSD should be non-negative"

    def test_centroid_has_valid_time(self, trajectory_info, reference_conf):
        """Test that centroid returns configuration with valid time."""
        top_info, traj_info = trajectory_info

        result_conf, _ = centroid(traj_info, top_info, reference_conf, ncpus=1)

        # Time should be set (from one of the trajectory configurations)
        assert result_conf.time is not None, "Time should be set"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_arguments(self):
        """Test that parser requires reference and trajectory arguments."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-p", "4",
            "-o", "output.dat",
            "-i", "index.txt",
            "-q",
            "reference.dat",
            "trajectory.dat"
        ])

        assert args.reference_structure == ["reference.dat"], "Reference not parsed"
        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.output == ["output.dat"], "Output option not parsed"
        assert args.index_file == ["index.txt"], "Index file option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["ref.dat", "traj.dat"])
        assert args_defaults.quiet is False, "Quiet should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output_file(self, mini_traj_path, mean_conf_path, temp_output_dir, monkeypatch):
        """Test main() creates output file with correct properties."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        ref_copy = temp_output_dir / "ref.dat"
        copy(mini_traj_path, traj_copy)
        copy(mean_conf_path, ref_copy)

        output_file = temp_output_dir / "test_centroid.dat"

        test_args = ["centroid.py", "-o", str(output_file), str(ref_copy), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check file exists and has content
        assert output_file.exists(), "Output file should be created"
        assert output_file.stat().st_size > 0, "Output file should have content"

        # Verify it's a valid configuration
        top_info, traj_info = describe(None, str(output_file))
        assert traj_info.nconfs == 1, "Output should have exactly one configuration"

    def test_main_default_output_filename(self, mini_traj_path, mean_conf_path, temp_output_dir, monkeypatch):
        """Test that main() uses default output filename 'centroid.dat'."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        ref_copy = temp_output_dir / "ref.dat"
        copy(mini_traj_path, traj_copy)
        copy(mean_conf_path, ref_copy)

        test_args = ["centroid.py", str(ref_copy), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        default_output = temp_output_dir / "centroid.dat"
        assert default_output.exists(), "Default output 'centroid.dat' should be created"

    def test_main_with_index_file(self, mini_traj_path, mean_conf_path, temp_output_dir, monkeypatch):
        """Test main() with index file for subset alignment."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        ref_copy = temp_output_dir / "ref.dat"
        copy(mini_traj_path, traj_copy)
        copy(mean_conf_path, ref_copy)

        # Create index file with subset of particles
        index_file = temp_output_dir / "index.txt"
        index_file.write_text("0 1 2 3 4 5 6 7 8 9")

        output_file = temp_output_dir / "indexed_centroid.dat"

        test_args = ["centroid.py", "-i", str(index_file), "-o", str(output_file), str(ref_copy), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Index file option should work"
