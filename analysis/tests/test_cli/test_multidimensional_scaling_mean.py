"""
Tests for oxDNA_analysis_tools.multidimensional_scaling_mean module.

Tests cover:
- make_heatmap() visualization function
- multidimensional_scaling_mean() API function
- distance_deviations() API function
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

from oxDNA_analysis_tools.multidimensional_scaling_mean import (
    multidimensional_scaling_mean,
    distance_deviations,
    cli_parser,
    main
)
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs
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
def trajectory_info(mini_traj_path):
    """Get topology and trajectory info for the mini trajectory."""
    top_info, traj_info = describe(None, str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# API Tests - multidimensional_scaling_mean() function
# =============================================================================

class TestMdsFunction:
    """Tests for the multidimensional_scaling_mean() API function."""

    def test_mds_returns_configuration_and_mask(self, trajectory_info):
        """Test that MDS returns a Configuration and masked mean."""
        top_info, traj_info = trajectory_info

        mean_conf, masked_mean = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)

        assert isinstance(mean_conf, Configuration), "Should return Configuration"
        assert isinstance(masked_mean, np.ma.MaskedArray), "Should return masked array"

    def test_mds_configuration_shape(self, trajectory_info):
        """Test that MDS configuration has correct shape."""
        top_info, traj_info = trajectory_info

        mean_conf, _ = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)

        assert mean_conf.positions.shape == (top_info.nbases, 3), "Positions shape mismatch"
        assert mean_conf.a1s.shape == (top_info.nbases, 3), "a1s shape mismatch"
        assert mean_conf.a3s.shape == (top_info.nbases, 3), "a3s shape mismatch"

    def test_mds_masked_mean_shape(self, trajectory_info):
        """Test that masked mean has correct shape."""
        top_info, traj_info = trajectory_info

        _, masked_mean = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)

        assert masked_mean.shape == (top_info.nbases, top_info.nbases), "Masked mean shape mismatch"

    def test_mds_positions_finite(self, trajectory_info):
        """Test that MDS positions are finite values."""
        top_info, traj_info = trajectory_info

        mean_conf, _ = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)

        assert np.all(np.isfinite(mean_conf.positions)), "All positions should be finite"


# =============================================================================
# API Tests - distance_deviations() function
# =============================================================================

class TestDistanceDeviationsFunction:
    """Tests for the distance_deviations() API function."""

    def test_deviations_returns_array(self, trajectory_info):
        """Test that distance_deviations returns an array."""
        top_info, traj_info = trajectory_info

        _, masked_mean = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)
        devs = distance_deviations(traj_info, top_info, masked_mean, ncpus=1)

        assert isinstance(devs, np.ndarray), "Should return numpy array"

    def test_deviations_shape(self, trajectory_info):
        """Test that deviations have correct shape."""
        top_info, traj_info = trajectory_info

        _, masked_mean = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)
        devs = distance_deviations(traj_info, top_info, masked_mean, ncpus=1)

        assert len(devs) == top_info.nbases, "Should have one deviation per particle"

    def test_deviations_non_negative(self, trajectory_info):
        """Test that deviations are non-negative (RMSD values)."""
        top_info, traj_info = trajectory_info

        _, masked_mean = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)
        devs = distance_deviations(traj_info, top_info, masked_mean, ncpus=1)

        # Deviations should be non-negative (they're squared then sqrt'd)
        assert np.all(devs >= 0), "Deviations should be non-negative"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_trajectory(self):
        """Test that parser requires trajectory."""
        parser = cli_parser()

        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options."""
        parser = cli_parser()

        args = parser.parse_args([
            "-o", "mean.dat",
            "-d", "devs.json",
            "-p", "4",
            "-q",
            "trajectory.dat"
        ])

        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.output == ["mean.dat"], "Output option not parsed"
        assert args.dev_file == ["devs.json"], "Dev file option not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

    def test_parser_defaults(self):
        """Test parser default values."""
        parser = cli_parser()

        args = parser.parse_args(["traj.dat"])

        assert args.output is None
        assert args.dev_file is None
        assert args.parallel is None
        assert args.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_mean_file(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates mean configuration file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "mean.dat"

        test_args = [
            "multidimensional_scaling_mean.py",
            "-o", str(output_file),
            str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Mean file should be created"

    def test_main_default_outputs(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() uses default output filenames."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        test_args = [
            "multidimensional_scaling_mean.py",
            str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        default_mean = temp_output_dir / "mean_mds.dat"
        default_devs = temp_output_dir / "devs_mds.json"

        assert default_mean.exists(), "Default mean 'mean_mds.dat' should be created"
        assert default_devs.exists(), "Default devs 'devs_mds.json' should be created"

    def test_main_creates_devs_file(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates deviations JSON file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        devs_file = temp_output_dir / "devs.json"

        test_args = [
            "multidimensional_scaling_mean.py",
            "-d", str(devs_file),
            str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert devs_file.exists(), "Devs file should be created"

        # Verify JSON structure
        with open(devs_file) as f:
            data = json.load(f)

        assert "contact deviation" in data, "JSON should have 'contact deviation' key"
        assert isinstance(data["contact deviation"], list), "Deviations should be list"

    def test_main_mean_readable(self, mini_traj_path, temp_output_dir, monkeypatch, trajectory_info):
        """Test that mean output file is readable."""
        monkeypatch.chdir(temp_output_dir)
        top_info, _ = trajectory_info

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "mean.dat"

        test_args = [
            "multidimensional_scaling_mean.py",
            "-o", str(output_file),
            str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Verify the output is readable
        mean_top, mean_traj = describe(None, str(output_file))
        assert mean_traj.nconfs == 1, "Mean should have 1 configuration"

        mean_conf = get_confs(mean_top, mean_traj, 0, 1)[0]
        assert mean_conf.positions.shape == (top_info.nbases, 3), "Mean should have correct particle count"

    def test_main_quiet_mode(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() with quiet mode."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "quiet_mean.dat"

        test_args = [
            "multidimensional_scaling_mean.py",
            "-q",
            "-o", str(output_file),
            str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Quiet mode should still create output"


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_mds_preserves_box(self, trajectory_info):
        """Test that MDS preserves box from example configuration."""
        top_info, traj_info = trajectory_info

        mean_conf, _ = multidimensional_scaling_mean(traj_info, top_info, ncpus=1)

        # Box should be preserved from example conf
        assert mean_conf.box is not None, "Box should be set"
        assert len(mean_conf.box) == 3, "Box should have 3 dimensions"
