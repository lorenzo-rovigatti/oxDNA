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
from copy import deepcopy
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


# =============================================================================
# Known-Value Validation Tests
# =============================================================================

class TestKnownValueValidation:
    """Tests validating actual RMSD/RMSF values against expected results."""

    def test_rmsd_calculation_correctness(self, trajectory_info, mean_conf):
        """
        Test RMSD calculation produces correct values.

        RMSD should be the root mean square deviation of each configuration from the mean.
        Validates calculation by computing RMSD manually and comparing.
        """
        from oxDNA_analysis_tools.align import svd_align

        top_info, traj_info = trajectory_info

        # Get deviations output
        RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=1)

        # Manually compute RMSD for verification
        from oxDNA_analysis_tools.UTILS.RyeReader import inbox, get_confs

        # Prepare mean configuration the same way deviations() does
        mean_centered = inbox(mean_conf, center=True)
        indexes = np.arange(top_info.nbases)
        ref_cms = np.mean(mean_centered.positions[indexes], axis=0)
        mean_centered.positions -= ref_cms

        # Load and align configurations manually
        confs = get_confs(top_info, traj_info, 0, traj_info.nconfs)
        manual_rmsds = []

        for conf in confs:
            # Inbox and center
            conf_centered = inbox(conf, center=True)

            # Align to mean
            aligned_coords = svd_align(
                mean_centered.positions[indexes],
                np.array([conf_centered.positions, conf_centered.a1s, conf_centered.a3s]),
                indexes,
                ref_center=np.zeros(3)
            )[0]

            # Compute squared deviations
            sq_devs = np.power(np.linalg.norm(aligned_coords - mean_centered.positions, axis=1), 2)

            # RMSD: sqrt(mean(squared deviations)) * 0.8518 nm
            rmsd = np.sqrt(np.mean(sq_devs)) * 0.8518
            manual_rmsds.append(rmsd)

        manual_rmsds = np.array(manual_rmsds)

        # Compare with deviations() output
        np.testing.assert_allclose(
            RMSDs, manual_rmsds, rtol=1e-5,
            err_msg="RMSD calculation does not match manual computation"
        )

    def test_rmsf_calculation_correctness(self, trajectory_info, mean_conf):
        """
        Test RMSF calculation produces correct values.

        RMSF should be the root mean square fluctuation per nucleotide.
        Validates calculation by computing RMSF manually and comparing.
        """
        from oxDNA_analysis_tools.align import svd_align

        top_info, traj_info = trajectory_info

        # Get deviations output
        RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=1)
        # Manually compute RMSF for verification
        from oxDNA_analysis_tools.UTILS.RyeReader import inbox, get_confs

        # Prepare mean configuration the same way deviations() does
        mean_centered = inbox(mean_conf, center=True)
        indexes = np.arange(top_info.nbases)
        ref_cms = np.mean(mean_centered.positions[indexes], axis=0)
        mean_centered.positions -= ref_cms

        # Load and align configurations manually
        confs = get_confs(top_info, traj_info, 0, traj_info.nconfs)
        all_sq_devs = []

        for conf in confs:
            # Inbox and center
            conf_centered = inbox(conf, center=True)

            # Align to mean
            aligned_coords = svd_align(
                mean_centered.positions[indexes],
                np.array([conf_centered.positions, conf_centered.a1s, conf_centered.a3s]),
                indexes,
                ref_center=np.zeros(3)
            )[0]

            # Compute squared deviations per particle
            sq_devs = np.power(np.linalg.norm(aligned_coords - mean_centered.positions, axis=1), 2)
            all_sq_devs.append(sq_devs)

        all_sq_devs = np.array(all_sq_devs)

        # RMSF: sqrt(mean across configs of squared deviations) * 0.8518 nm
        manual_rmsfs = np.sqrt(np.mean(all_sq_devs, axis=0)) * 0.8518

        # Compare with deviations() output
        np.testing.assert_allclose(
            RMSFs, manual_rmsfs, rtol=1e-5,
            err_msg="RMSF calculation does not match manual computation"
        )

# =============================================================================
# Output Content Validation Tests
# =============================================================================

class TestOutputContentValidation:
    """Tests validating the structure and content of JSON output files."""

    def test_rmsf_json_content_structure(self, trajectory_info, mean_conf, temp_output_dir):
        """Test RMSF JSON file contains correct structure and valid values."""
        top_info, traj_info = trajectory_info

        RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=1)

        outfile = str(temp_output_dir / "rmsf_content.json")
        output(RMSDs, RMSFs, outfile,
               str(temp_output_dir / "plot.png"),
               str(temp_output_dir / "data.json"))

        # Load and validate JSON structure
        try:
            with open(outfile) as f:
                rmsf_data = json.load(f)
        except json.JSONDecodeError as e:
            pytest.fail(f"RMSF file is not valid JSON: {e}")

        # Check keys
        assert "RMSF (nm)" in rmsf_data, "RMSF JSON should have 'RMSF (nm)' key"

        # Check data is a list
        assert isinstance(rmsf_data["RMSF (nm)"], list), "RMSF data should be a list"

        # Check length matches number of bases
        assert len(rmsf_data["RMSF (nm)"]) == top_info.nbases, \
            f"RMSF should have {top_info.nbases} values, got {len(rmsf_data['RMSF (nm)'])}"

        # Check all values are numeric and non-negative
        rmsf_values = np.array(rmsf_data["RMSF (nm)"])
        assert np.all(np.isfinite(rmsf_values)), "All RMSF values should be finite"
        assert np.all(rmsf_values >= 0), "All RMSF values should be non-negative"

        # Check values match computed RMSFs
        np.testing.assert_allclose(
            rmsf_values, RMSFs, rtol=1e-10,
            err_msg="RMSF values in JSON don't match computed values"
        )

    def test_rmsd_json_content_structure(self, trajectory_info, mean_conf, temp_output_dir):
        """Test RMSD JSON file contains correct structure and valid values."""
        top_info, traj_info = trajectory_info

        # Slice trajectory because we can
        traj_info_copy = deepcopy(traj_info)
        traj_info_copy.idxs = traj_info.idxs[0:4]
        traj_info_copy.nconfs = 4

        RMSDs, RMSFs = deviations(traj_info_copy, top_info, mean_conf, ncpus=1)

        data_file = str(temp_output_dir / "rmsd_content.json")
        output(RMSDs, RMSFs,
               str(temp_output_dir / "rmsf.json"),
               str(temp_output_dir / "plot.png"),
               data_file)

        # Load and validate JSON structure
        try:
            with open(data_file) as f:
                rmsd_data = json.load(f)
        except json.JSONDecodeError as e:
            pytest.fail(f"RMSD file is not valid JSON: {e}")

        # Check keys
        assert "RMSD (nm)" in rmsd_data, "RMSD JSON should have 'RMSD (nm)' key"

        # Check data is a list
        assert isinstance(rmsd_data["RMSD (nm)"], list), "RMSD data should be a list"

        # Check length matches number of configurations
        assert len(rmsd_data["RMSD (nm)"]) == traj_info_copy.nconfs, \
            f"RMSD should have {traj_info_copy.nconfs} values, got {len(rmsd_data['RMSD (nm)'])}"

        # Check all values are numeric and non-negative
        rmsd_values = np.array(rmsd_data["RMSD (nm)"])
        assert np.all(np.isfinite(rmsd_values)), "All RMSD values should be finite"
        assert np.all(rmsd_values >= 0), "All RMSD values should be non-negative"

        # Check values match computed RMSDs
        np.testing.assert_allclose(
            rmsd_values, RMSDs, rtol=1e-10,
            err_msg="RMSD values in JSON don't match computed values"
        )

# =============================================================================
# Parallel Processing Tests
# =============================================================================

class TestParallelProcessing:
    """Tests for parallel processing functionality."""

    def test_parallel_execution_completes(self, trajectory_info, mean_conf):
        """Test that parallel processing completes without errors."""
        top_info, traj_info = trajectory_info

        # Run with ncpus=2 (should not crash)
        RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=2)

        # Check that we get results with correct shapes
        assert isinstance(RMSDs, np.ndarray), "Parallel execution should return RMSD array"
        assert isinstance(RMSFs, np.ndarray), "Parallel execution should return RMSF array"
        assert RMSDs.shape == (traj_info.nconfs,), "Parallel RMSD should have correct shape"
        assert RMSFs.shape == (top_info.nbases,), "Parallel RMSF should have correct shape"

    def test_serial_execution_stability(self, trajectory_info, mean_conf):
        """Test that serial execution produces consistent results when run twice."""
        top_info, traj_info = trajectory_info

        # Run twice with ncpus=1
        RMSDs_1, RMSFs_1 = deviations(traj_info, top_info, mean_conf, ncpus=1)
        RMSDs_2, RMSFs_2 = deviations(traj_info, top_info, mean_conf, ncpus=1)

        # Serial execution should be completely deterministic
        np.testing.assert_array_equal(
            RMSDs_1, RMSDs_2,
            err_msg="Serial execution should be deterministic for RMSDs"
        )

        np.testing.assert_array_equal(
            RMSFs_1, RMSFs_2,
            err_msg="Serial execution should be deterministic for RMSFs"
        )

    def test_parallel_with_different_cpu_counts(self, trajectory_info, mean_conf):
        """Test that parallel processing works with different CPU counts."""
        top_info, traj_info = trajectory_info

        # Get baseline with ncpus=1
        RMSDs_1, RMSFs_1 = deviations(traj_info, top_info, mean_conf, ncpus=1)

        # Test with ncpus=2, and 4
        for ncpus in [2, 4]:
            RMSDs, RMSFs = deviations(traj_info, top_info, mean_conf, ncpus=ncpus)
            # Check shapes are consistent regardless of CPU count
            assert RMSDs.shape == (traj_info.nconfs,), \
                f"RMSD shape should be consistent with ncpus={ncpus}"
            assert RMSFs.shape == (top_info.nbases,), \
                f"RMSF shape should be consistent with ncpus={ncpus}"
            
            np.testing.assert_array_equal(
                RMSDs_1, RMSDs,
            err_msg=f"Parallel execution with ncpus={ncpus} should be deterministic for RMSDs"
            )

            np.testing.assert_array_equal(
                RMSFs_1, RMSFs,
            err_msg=f"Parallel execution with ncpus={ncpus} should be deterministic for RMSFs"
            )
