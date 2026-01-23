"""
Tests for oxDNA_analysis_tools.deviations module.

Tests cover:
- deviations() API function
- compute() helper function
- output() function
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

from oxDNA_analysis_tools.deviations import deviations, compute, output, cli_parser, ComputeContext, main
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
def mean_conf(trajectory_info, mean_conf_path):
    """Get mean configuration for deviation calculation."""
    top_info, _ = trajectory_info
    _, mean_info = describe(None, str(mean_conf_path))
    return get_confs(top_info, mean_info, 0, 1)[0]


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# API Tests - deviations() function
# =============================================================================

class TestDeviationsFunction:
    """Tests for the deviations() API function."""

    def test_deviations_behavior(self, trajectory_info, mean_conf):
        """Test deviations() returns correct shapes, non-negative values, reasonable scale, and handles indexes."""
        top_info, traj_info = trajectory_info

        # Basic call
        RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=1)

        # Check return types
        assert isinstance(RMSDs, np.ndarray), "RMSDs should be numpy array"
        assert isinstance(RMSFs, np.ndarray), "RMSFs should be numpy array"

        # Check shapes
        assert RMSDs.shape == (traj_info.nconfs,), "RMSDs should have one value per conf"
        assert RMSFs.shape == (top_info.nbases,), "RMSFs should have one value per particle"

        # All values should be non-negative
        assert np.all(RMSDs >= 0), "RMSDs should be non-negative"
        assert np.all(RMSFs >= 0), "RMSFs should be non-negative"

        # Values should be in nm scale (typical DNA fluctuations)
        assert np.mean(RMSDs) < 50, "Mean RMSD should be less than 50 nm"
        assert np.mean(RMSFs) < 50, "Mean RMSF should be less than 50 nm"

        # Test with custom indexes
        indexes = list(range(top_info.nbases // 2))
        RMSDs_idx, RMSFs_idx = deviations(traj_info, top_info, mean_conf, indexes=indexes, ncpus=1)

        # Should still return arrays for all configs and particles
        assert RMSDs_idx.shape == (traj_info.nconfs,), "RMSDs shape unchanged with indexes"
        assert RMSFs_idx.shape == (top_info.nbases,), "RMSFs shape unchanged with indexes"


# =============================================================================
# Output Function Tests
# =============================================================================

class TestOutputFunction:
    """Tests for the output() function."""

    def test_output_behavior(self, trajectory_info, mean_conf, temp_output_dir):
        """Test output() creates files with valid JSON structure and expected keys."""
        import matplotlib
        matplotlib.use('Agg')

        top_info, traj_info = trajectory_info
        RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=1)

        outfile = str(temp_output_dir / "devs.json")
        plot_name = str(temp_output_dir / "rmsd.png")
        data_file = str(temp_output_dir / "rmsd_op.json")

        output(RMSDs, RMSFs, outfile, plot_name, data_file)

        # Check files exist
        assert Path(outfile).exists(), "RMSF JSON file should be created"
        assert Path(plot_name).exists(), "RMSD plot should be created"
        assert Path(data_file).exists(), "RMSD data file should be created"

        # Check RMSF file structure
        with open(outfile) as f:
            rmsf_data = json.load(f)
        assert "RMSF (nm)" in rmsf_data, "RMSF file should have 'RMSF (nm)' key"
        assert len(rmsf_data["RMSF (nm)"]) == top_info.nbases, "RMSF should have nbases values"

        # Check RMSD file structure
        with open(data_file) as f:
            rmsd_data = json.load(f)
        assert "RMSD (nm)" in rmsd_data, "RMSD file should have 'RMSD (nm)' key"
        assert len(rmsd_data["RMSD (nm)"]) == traj_info.nconfs, "RMSD should have nconfs values"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_arguments(self):
        """Test that parser requires mean structure and trajectory arguments."""
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
            "-i", "index.txt",
            "-r", "rmsd.png",
            "-d", "data.json",
            "-q",
            "mean.dat",
            "trajectory.dat"
        ])

        assert args.mean_structure == ["mean.dat"], "Mean structure not parsed"
        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.output == ["output.json"], "Output option not parsed"
        assert args.index_file == ["index.txt"], "Index file option not parsed"
        assert args.rmsd_plot == ["rmsd.png"], "RMSD plot option not parsed"
        assert args.rmsd_data == ["data.json"], "RMSD data option not parsed"
        assert args.quiet is True, "Quiet option not parsed"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output_files(self, mini_traj_path, mean_conf_path, temp_output_dir, monkeypatch):
        """Test main() creates all output files."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        mean_copy = temp_output_dir / "mean.dat"
        copy(mini_traj_path, traj_copy)
        copy(mean_conf_path, mean_copy)

        output_file = temp_output_dir / "test_devs.json"
        plot_file = temp_output_dir / "test_rmsd.png"
        data_file = temp_output_dir / "test_data.json"

        test_args = [
            "deviations.py",
            "-o", str(output_file),
            "-r", str(plot_file),
            "-d", str(data_file),
            str(mean_copy),
            str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "RMSF output file should be created"
        assert plot_file.exists(), "RMSD plot should be created"
        assert data_file.exists(), "RMSD data file should be created"

    def test_main_default_filenames(self, mini_traj_path, mean_conf_path, temp_output_dir, monkeypatch):
        """Test main() uses default filenames."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        mean_copy = temp_output_dir / "mean.dat"
        copy(mini_traj_path, traj_copy)
        copy(mean_conf_path, mean_copy)

        test_args = ["deviations.py", str(mean_copy), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert (temp_output_dir / "devs.json").exists(), "Default RMSF file should be created"
        assert (temp_output_dir / "rmsd.png").exists(), "Default plot should be created"
        assert (temp_output_dir / "rmsd_op.json").exists(), "Default data file should be created"
