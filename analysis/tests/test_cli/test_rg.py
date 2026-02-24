"""
Tests for oxDNA_analysis_tools.rg module (radius of gyration).

Tests cover:
- rg() API function
- compute() helper function
- CLI argument parsing
- main() CLI entry point
- Mathematical correctness of Rg calculation
- Scaling factor validation
- Index subsetting functionality
- Parallel vs serial consistency
"""
import json
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.rg import rg, compute, cli_parser, ComputeContext, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, inbox
from oxDNA_analysis_tools.UTILS.data_structures import TrajInfo


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

    def test_compute_basic_behavior(self, trajectory_info):
        """Test compute() returns correct array with positive Rg values."""
        top_info, traj_info = trajectory_info

        # Create context with empty indexes (use all particles)
        ctx = ComputeContext(
            top_info=top_info,
            traj_info=traj_info,
            indexes=np.arange(top_info.nbases)
        )

        # Process first chunk
        result = compute(ctx, chunk_size=1, chunk_id=0)

        # Progressive validation
        # 1. Correct type
        assert isinstance(result, np.ndarray), "compute() should return numpy array"

        # 2. Correct shape
        assert len(result) == 1, "Single chunk should return one Rg value"

        # 3. Reasonable values - Rg should be positive
        assert result[0] > 0, "Radius of gyration should be positive"

        # Additional sanity check - for DNA structures, Rg should be in reasonable range
        # Typical DNA structures are a few nanometers
        assert result[0] < 100, "Rg should be less than 100 nm for typical structures"

        # Process second chunk to verify consistency
        result2 = compute(ctx, chunk_size=1, chunk_id=1)
        assert isinstance(result2, np.ndarray), "Second chunk should return array"
        assert result2[0] > 0, "Second Rg should be positive"

    def test_compute_with_indexes_subset(self, trajectory_info):
        """Test compute() with subset of particles via indexes parameter."""
        top_info, traj_info = trajectory_info

        # Create context with subset of particles (first 10)
        subset_indexes = np.arange(10)
        ctx_subset = ComputeContext(
            top_info=top_info,
            traj_info=traj_info,
            indexes=subset_indexes
        )

        result_subset = compute(ctx_subset, chunk_size=1, chunk_id=0)

        # Should still return valid Rg value
        assert isinstance(result_subset, np.ndarray), "Should return numpy array"
        assert len(result_subset) == 1, "Should return one Rg value"
        assert result_subset[0] > 0, "Rg should be positive"

        # Create context with all particles
        ctx_all = ComputeContext(
            top_info=top_info,
            traj_info=traj_info,
            indexes=np.arange(top_info.nbases)
        )

        result_all = compute(ctx_all, chunk_size=1, chunk_id=0)

        # Subset Rg should be different from full system Rg
        assert not np.isclose(result_subset[0], result_all[0]), \
            "Subset Rg should differ from full system Rg"


# =============================================================================
# Mathematical Correctness Tests
# =============================================================================

class TestMathematicalCorrectness:
    """Tests that verify the mathematical correctness of Rg calculation."""

    def test_rg_formula_correctness(self, trajectory_info):
        """
        Test that Rg is correctly computed as sqrt(mean(dist^2)) from center of mass.

        Manually compute Rg for a single configuration and compare to function output.
        Formula: Rg = sqrt(mean((r_i - r_com)^2)) * scaling_factor
        where r_com is the center of mass.
        """
        top_info, traj_info = trajectory_info

        # Use a single configuration for deterministic test
        traj_info_single = TrajInfo(
            traj_info.path,
            nconfs=1,
            idxs=traj_info.idxs[0:1],
            incl_v=traj_info.incl_v
        )

        # Get the configuration
        conf = get_confs(top_info, traj_info_single, 0, 1)[0]
        conf = inbox(conf)

        # Manually compute Rg
        positions = conf.positions
        com = np.mean(positions, axis=0)
        distances = np.linalg.norm(positions - com, axis=1)
        sq_distances = np.power(distances, 2)
        mean_sq_distance = np.sum(sq_distances) / top_info.nbases
        manual_rg = np.sqrt(mean_sq_distance) * 0.8518  # Apply scaling factor

        # Compute using function
        result = rg(top_info, traj_info_single, indexes=np.array([]), ncpus=1)

        # 4. Correct content - compare with manually computed value
        np.testing.assert_allclose(
            result[0],
            manual_rg,
            rtol=1e-10,
            err_msg="Function Rg doesn't match manually computed Rg"
        )

    def test_rg_scaling_factor(self, trajectory_info):
        """
        Test that the 0.8518 scaling factor is correctly applied.

        Verify that output is in nanometers by comparing scaled vs unscaled values.
        """
        top_info, traj_info = trajectory_info

        # Use single configuration
        traj_info_single = TrajInfo(
            traj_info.path,
            nconfs=1,
            idxs=traj_info.idxs[0:1],
            incl_v=traj_info.incl_v
        )

        # Get configuration and compute Rg in simulation units (without scaling)
        conf = get_confs(top_info, traj_info_single, 0, 1)[0]
        conf = inbox(conf)
        positions = conf.positions
        com = np.mean(positions, axis=0)
        distances = np.linalg.norm(positions - com, axis=1)
        sq_distances = np.power(distances, 2)
        mean_sq_distance = np.sum(sq_distances) / top_info.nbases
        rg_sim_units = np.sqrt(mean_sq_distance)

        # Compute using function (should be in nm)
        result = rg(top_info, traj_info_single, indexes=np.array([]), ncpus=1)

        # Verify scaling: result should equal sim_units * 0.8518
        expected_nm = rg_sim_units * 0.8518
        np.testing.assert_allclose(
            result[0],
            expected_nm,
            rtol=1e-10,
            err_msg="Scaling factor 0.8518 not correctly applied"
        )

        # Also verify the ratio
        ratio = result[0] / rg_sim_units
        np.testing.assert_allclose(
            ratio,
            0.8518,
            rtol=1e-10,
            err_msg="Ratio of nm to sim units should be 0.8518"
        )

    def test_rg_with_indexes_mathematical_correctness(self, trajectory_info):
        """
        Test that Rg computed with subset of particles is mathematically correct.

        Manually compute Rg for a subset and compare to function output.
        """
        top_info, traj_info = trajectory_info

        # Use single configuration
        traj_info_single = TrajInfo(
            traj_info.path,
            nconfs=1,
            idxs=traj_info.idxs[0:1],
            incl_v=traj_info.incl_v
        )

        # Define subset (first 20 particles)
        subset_indexes = np.arange(20)

        # Get configuration
        conf = get_confs(top_info, traj_info_single, 0, 1)[0]
        conf = inbox(conf)

        # Manually compute Rg for subset
        subset_positions = conf.positions[subset_indexes]
        com_subset = np.mean(subset_positions, axis=0)
        distances_subset = np.linalg.norm(subset_positions - com_subset, axis=1)
        sq_distances_subset = np.power(distances_subset, 2)
        mean_sq_distance_subset = np.sum(sq_distances_subset) / len(subset_indexes)
        manual_rg_subset = np.sqrt(mean_sq_distance_subset) * 0.8518

        # Compute using function
        result = rg(top_info, traj_info_single, indexes=subset_indexes, ncpus=1)

        # 1. Correct type
        assert isinstance(result, np.ndarray), "Should return numpy array"

        # 2. Correct shape
        assert result.shape == (traj_info_single.nconfs,), \
            f"Should have one Rg per configuration, got shape {result.shape}"

        # 3. Reasonable values
        assert np.all(result > 0), "All Rg values should be positive"
        assert np.all(result < 100), "Rg values should be less than 100 nm"
        assert not np.any(np.isnan(result)), "Should not contain NaN values"
        assert not np.any(np.isinf(result)), "Should not contain infinite values"

        # 4. Exact values
        np.testing.assert_allclose(
            result[0],
            manual_rg_subset,
            rtol=1e-10,
            err_msg="Subset Rg doesn't match manually computed value"
        )


# =============================================================================
# API Tests - rg() function
# =============================================================================

class TestRgFunction:
    """Tests for the rg() API function."""

    def test_rg_with_indexes_subset(self, trajectory_info):
        """Test rg() with subset of particles specified via indexes parameter."""
        top_info, traj_info = trajectory_info

        # Use first quarter of particles
        subset_indexes = np.array(list(range(top_info.nbases // 4)))

        # Execute
        result_subset = rg(top_info, traj_info, indexes=subset_indexes, ncpus=1)

        # Validate
        assert isinstance(result_subset, np.ndarray), "Should return numpy array"
        assert result_subset.shape == (traj_info.nconfs,), "Should have one Rg per conf"
        assert np.all(result_subset > 0), "All Rg values should be positive"

        # Compare with full system
        result_full = rg(top_info, traj_info, indexes=np.array([]), ncpus=1)

        # Subset should generally be smaller (but not guaranteed for all structures)
        # Just verify they're different
        assert not np.allclose(result_subset, result_full), \
            "Subset Rg should differ from full system Rg"

    def test_rg_empty_indexes_uses_all_particles(self, trajectory_info):
        """Test that empty indexes array uses all particles."""
        top_info, traj_info = trajectory_info

        # Use deterministic subset
        traj_info_subset = TrajInfo(
            traj_info.path,
            nconfs=3,
            idxs=traj_info.idxs[0:3],
            incl_v=traj_info.incl_v
        )

        result_empty = rg(top_info, traj_info_subset, indexes=np.array([]), ncpus=1)
        result_all = rg(top_info, traj_info_subset,
                       indexes=np.array(list(range(top_info.nbases))), ncpus=1)

        # Both should produce identical results
        np.testing.assert_allclose(
            result_empty,
            result_all,
            rtol=1e-10,
            err_msg="Empty indexes should behave same as all indexes"
        )


# =============================================================================
# Parallel Processing Consistency Tests
# =============================================================================

class TestParallelConsistency:
    """Tests that verify parallel processing produces same results as serial."""

    def test_parallel_vs_serial_consistency(self, trajectory_info):
        """
        Test that ncpus=1 and ncpus=2 produce identical results.

        This ensures parallelization doesn't introduce numerical differences.
        """
        top_info, traj_info = trajectory_info

        # Use deterministic subset for faster testing
        traj_info_subset = TrajInfo(
            traj_info.path,
            nconfs=6,  # Divisible by 2 for clean chunking
            idxs=traj_info.idxs[0:6],
            incl_v=traj_info.incl_v
        )

        # Execute both serial and parallel
        result_serial = rg(top_info, traj_info_subset, indexes=np.array([]), ncpus=1)
        result_parallel = rg(top_info, traj_info_subset, indexes=np.array([]), ncpus=2)

        # Results should be identical within floating point precision
        np.testing.assert_allclose(
            result_serial,
            result_parallel,
            rtol=1e-12,
            atol=1e-14,
            err_msg="Parallel and serial processing should produce identical results"
        )

    def test_parallel_consistency_with_indexes(self, trajectory_info):
        """Test parallel consistency when using subset of indexes."""
        top_info, traj_info = trajectory_info

        # Use subset of configurations
        traj_info_subset = TrajInfo(
            traj_info.path,
            nconfs=4,
            idxs=traj_info.idxs[0:4],
            incl_v=traj_info.incl_v
        )

        # Use subset of particles
        subset_indexes = np.array(list(range(top_info.nbases // 2)))

        result_serial = rg(top_info, traj_info_subset, indexes=subset_indexes, ncpus=1)
        result_parallel = rg(top_info, traj_info_subset, indexes=subset_indexes, ncpus=2)

        np.testing.assert_allclose(
            result_serial,
            result_parallel,
            rtol=1e-12,
            atol=1e-14,
            err_msg="Parallel consistency should hold with subset indexes"
        )


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
            "-o", "output.png",
            "-d", "data.json",
            "-i", "index.txt",
            "-q",
            "trajectory.dat"
        ])

        assert args.trajectory == "trajectory.dat", "Trajectory not parsed"
        assert args.parallel == 4, "Parallel option not parsed"
        assert args.output == "output.png", "Output option not parsed"
        assert args.data == "data.json", "Data option not parsed"
        assert args.index_file == ["index.txt"], "Index file option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["traj.dat"])
        assert args_defaults.quiet is False, "Quiet should default to False"
        assert args_defaults.parallel is None, "Parallel should default to None"
        assert args_defaults.output is None, "Output should default to None"
        assert args_defaults.data is None, "Data should default to None"
        assert args_defaults.index_file is None, "Index file should default to None"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_basic_output(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates output files with correct structure and values."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Test with explicit output files
        data_file = temp_output_dir / "test_rg.json"
        plot_file = temp_output_dir / "test_rg.png"
        test_args = ["rg.py", "-d", str(data_file), "-o", str(plot_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Progressive validation of data file
        # 1. File exists
        assert data_file.exists(), "Data file should be created"

        # 2. Valid JSON structure
        with open(data_file) as f:
            data = json.load(f)

        assert "rg" in data, "JSON should have 'rg' key"

        # 3. Correct shape
        top_info, traj_info = describe(None, str(traj_copy))
        assert len(data["rg"]) == traj_info.nconfs, \
            "Should have one Rg per configuration"

        # 4. Correct content - all values should be positive and reasonable
        rg_values = np.array(data["rg"])
        assert np.all(rg_values > 0), "All Rg values should be positive"
        assert np.all(rg_values < 100), "All Rg values should be reasonable"

        # Verify plot file was created
        assert plot_file.exists(), "Plot file should be created"

    def test_main_default_output_filenames(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that main() uses default output filenames."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Run without specifying data file (but must specify plot file due to bug)
        # TODO: Fix bug in rg.py where plotfile is not defined when -o is omitted
        plot_file = temp_output_dir / "plot.png"
        test_args = ["rg.py", "-o", str(plot_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check default data filename
        default_data = temp_output_dir / "rg.json"
        assert default_data.exists(), "Default data file 'rg.json' should be created"

        # Verify content
        with open(default_data) as f:
            data = json.load(f)
        assert "rg" in data, "Default output should have correct structure"

    def test_main_with_index_file(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() with index file for subset of particles."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Create index file with subset of particles
        index_file = temp_output_dir / "index.txt"
        index_file.write_text("0 1 2 3 4 5 6 7 8 9")

        data_file = temp_output_dir / "indexed_rg.json"
        plot_file = temp_output_dir / "indexed_rg.png"
        test_args = ["rg.py", "-i", str(index_file), "-d", str(data_file), "-o", str(plot_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert data_file.exists(), "Index file option should work"

        # Verify output
        with open(data_file) as f:
            data = json.load(f)
        assert "rg" in data, "Output should have correct structure"
        assert all(v > 0 for v in data["rg"]), "All Rg values should be positive"

    def test_main_quiet_mode(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that main() respects quiet mode."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        data_file = temp_output_dir / "quiet_rg.json"
        plot_file = temp_output_dir / "quiet_rg.png"
        test_args = ["rg.py", "-q", "-d", str(data_file), "-o", str(plot_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert data_file.exists(), "Quiet mode should still create output"

    def test_main_invalid_index_file_raises(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that main() raises error for invalid index file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Create invalid index file (non-integer content)
        index_file = temp_output_dir / "bad_index.txt"
        index_file.write_text("not an integer")

        test_args = ["rg.py", "-i", str(index_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            with pytest.raises(RuntimeError, match="index file must be a space-seperated list"):
                main()


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_rg_single_configuration(self, trajectory_info):
        """Test rg() works with single configuration."""
        top_info, traj_info = trajectory_info

        # Create single-conf trajectory info
        single_conf_traj = TrajInfo(
            traj_info.path,
            nconfs=1,
            idxs=traj_info.idxs[0:1],
            incl_v=traj_info.incl_v
        )

        result = rg(top_info, single_conf_traj, indexes=np.array([]), ncpus=1)

        assert isinstance(result, np.ndarray), "Should return numpy array"
        assert result.shape == (1,), "Should have one Rg value"
        assert result[0] > 0, "Rg should be positive"

    def test_rg_high_ncpus(self, trajectory_info):
        """Test that rg() handles ncpus > nconfs gracefully."""
        top_info, traj_info = trajectory_info

        # Use more CPUs than configurations
        result = rg(top_info, traj_info, indexes=np.array([]), ncpus=100)

        assert isinstance(result, np.ndarray), "Should handle high ncpus gracefully"
        assert result.shape == (traj_info.nconfs,), "Should return correct number of values"
        assert np.all(result > 0), "All values should be positive"

    def test_rg_small_subset_indexes(self, trajectory_info):
        """Test rg() with very small subset of particles."""
        top_info, traj_info = trajectory_info

        # Use just 3 particles
        small_subset = np.array([0, 1, 2])

        result = rg(top_info, traj_info, indexes=small_subset, ncpus=1)

        assert isinstance(result, np.ndarray), "Should return numpy array"
        assert result.shape == (traj_info.nconfs,), "Should have one Rg per conf"
        assert np.all(result > 0), "All Rg values should be positive"
