"""
Tests for oxDNA_analysis_tools.mean module.

Tests cover:
- mean() API function
- CLI argument parsing
- main() CLI entry point
- compute() helper function
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.mean import mean, compute, cli_parser, ComputeContext, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, inbox, write_conf
from oxDNA_analysis_tools.UTILS.data_structures import Configuration, TrajInfo
from oxDNA_analysis_tools.align import svd_align


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
def expected_mean_path(test_resources):
    """Path to the expected mean configuration."""
    return test_resources / "mean.dat"


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
        """Test compute() returns correct array with proper shape for full and subset indexes."""
        top_info, traj_info = trajectory_info

        # Get a reference configuration for alignment
        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]
        ref_conf = inbox(ref_conf, center=True)

        # Test with all indexes
        indexes = list(range(top_info.nbases))
        reference_coords = ref_conf.positions[indexes]
        ref_cms = np.mean(reference_coords, axis=0)
        centered_ref_coords = reference_coords - ref_cms

        ctx = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            centered_ref_coords=centered_ref_coords,
            indexes=indexes
        )

        result = compute(ctx, chunk_size=1, chunk_id=0)

        # Should return array of shape [3, nbases, 3] (positions, a1s, a3s)
        assert isinstance(result, np.ndarray), "compute() should return numpy array"
        assert result.shape == (3, top_info.nbases, 3), "Shape should be (3, nbases, 3)"

        # Test with subset indexes - shape should still be full system
        subset_indexes = list(range(min(10, top_info.nbases)))
        subset_reference_coords = ref_conf.positions[subset_indexes]
        subset_ref_cms = np.mean(subset_reference_coords, axis=0)
        subset_centered_ref_coords = subset_reference_coords - subset_ref_cms

        ctx_subset = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            centered_ref_coords=subset_centered_ref_coords,
            indexes=subset_indexes
        )

        result_subset = compute(ctx_subset, chunk_size=1, chunk_id=0)
        assert result_subset.shape == (3, top_info.nbases, 3), "Subset alignment should still return full system"


# =============================================================================
# API Tests - mean() function
# =============================================================================

class TestMeanFunction:
    """Tests for the mean() API function."""

    def test_mean_basic_behavior(self, trajectory_info):
        """Test mean() returns correct Configuration with proper shapes and normalized orientations."""
        top_info, traj_info = trajectory_info

        result = mean(traj_info, top_info, ncpus=1)

        # Check return type
        assert isinstance(result, Configuration), "mean() should return Configuration"

        # Check shapes
        assert result.positions.shape == (top_info.nbases, 3), "Positions shape mismatch"
        assert result.a1s.shape == (top_info.nbases, 3), "a1s shape mismatch"
        assert result.a3s.shape == (top_info.nbases, 3), "a3s shape mismatch"

        # Check that a1 and a3 vectors are unit vectors (normalized)
        a1_norms = np.linalg.norm(result.a1s, axis=1)
        a3_norms = np.linalg.norm(result.a3s, axis=1)
        np.testing.assert_allclose(a1_norms, np.ones(top_info.nbases), rtol=1e-5,
                                   err_msg="a1 vectors should be normalized")
        np.testing.assert_allclose(a3_norms, np.ones(top_info.nbases), rtol=1e-5,
                                   err_msg="a3 vectors should be normalized")

    def test_mean_matches_expected_output(self, trajectory_info, expected_mean_path):
        """Test that mean() produces expected output (regression test)."""
        top_info, traj_info = trajectory_info

        # Compute mean with fixed reference (conf 0) for reproducibility
        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]
        ref_conf = inbox(ref_conf)

        result = mean(traj_info, top_info, ref_conf=ref_conf, ncpus=1)

        # Load expected result
        expected_ti, expected_di = describe(None, str(expected_mean_path))
        expected = get_confs(expected_ti, expected_di, 0, 1)[0]

        # Compare positions
        np.testing.assert_allclose(result.positions, expected.positions, rtol=1e-5,
                                   err_msg="Mean positions don't match expected")

    def test_mean_with_custom_ref_conf_and_box(self, trajectory_info):
        """Test mean() with custom reference configuration preserves box."""
        top_info, traj_info = trajectory_info

        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]
        ref_conf = inbox(ref_conf)

        result = mean(traj_info, top_info, ref_conf=ref_conf, ncpus=1)

        assert isinstance(result, Configuration), "Should return Configuration"
        assert result.positions.shape == (top_info.nbases, 3), "Wrong shape"
        np.testing.assert_array_equal(result.box, ref_conf.box,
                                      err_msg="Box should be preserved from reference")

    def test_mean_with_custom_indexes(self, trajectory_info):
        """Test mean() with subset of particle indexes for alignment."""
        top_info, traj_info = trajectory_info

        # Use only first half of particles for alignment
        indexes = list(range(top_info.nbases // 2))

        result = mean(traj_info, top_info, indexes=indexes, ncpus=1)

        # Should still return full configuration
        assert result.positions.shape == (top_info.nbases, 3), "Should return full system"

    def test_mean_empty_indexes_uses_all(self, trajectory_info):
        """Test that empty indexes list uses all particles."""
        top_info, traj_info = trajectory_info

        result_empty = mean(traj_info, top_info, indexes=[], ncpus=1)
        result_none = mean(traj_info, top_info, indexes=None, ncpus=1)

        # Both should produce same shape (using all particles)
        assert result_empty.positions.shape == result_none.positions.shape, \
            "Empty and None indexes should behave the same"


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
            "-o", "output.dat",
            "-d", "devs.json",
            "-i", "index.txt",
            "-a", "0",
            "-q",
            "trajectory.dat"
        ])

        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.output == ["output.dat"], "Output option not parsed"
        assert args.deviations == ["devs.json"], "Deviations option not parsed"
        assert args.index_file == ["index.txt"], "Index file option not parsed"
        assert args.align == ["0"], "Align option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["test.dat"])
        assert args_defaults.trajectory == ["test.dat"]
        assert args_defaults.quiet is False, "Quiet should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_basic_output(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates output file with correct properties."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "test_mean.dat"

        test_args = ["mean.py", "-o", str(output_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check file exists and has content
        assert output_file.exists(), "Output file should be created"
        assert output_file.stat().st_size > 0, "Output file should have content"

        # Verify it's a valid configuration with single conf
        top_info, traj_info = describe(None, str(output_file))
        assert traj_info.nconfs == 1, "Output should have exactly one configuration"

        conf = get_confs(top_info, traj_info, 0, 1)[0]
        assert conf is not None, "Should be able to read configuration"
        assert conf.positions is not None, "Configuration should have positions"

    def test_main_default_output_filename(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that main() uses default output filename 'mean.dat'."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        test_args = ["mean.py", str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        default_output = temp_output_dir / "mean.dat"
        assert default_output.exists(), "Default output 'mean.dat' should be created"

    def test_main_with_parallel(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() with parallel processing option."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "parallel_mean.dat"

        test_args = ["mean.py", "-p", "2", "-o", str(output_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Parallel processing should create output"

    def test_main_with_index_and_align(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() with index file and alignment options."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Create index file with subset of particles
        index_file = temp_output_dir / "index.txt"
        index_file.write_text("0 1 2 3 4 5 6 7 8 9")

        output_indexed = temp_output_dir / "indexed_mean.dat"
        output_aligned = temp_output_dir / "aligned_mean.dat"

        # Test with index file
        test_args_index = ["mean.py", "-i", str(index_file), "-o", str(output_indexed), str(traj_copy)]
        with patch.object(sys, 'argv', test_args_index):
            main()
        assert output_indexed.exists(), "Index file option should work"

        # Test with align option
        test_args_align = ["mean.py", "-a", "0", "-o", str(output_aligned), str(traj_copy)]
        with patch.object(sys, 'argv', test_args_align):
            main()
        assert output_aligned.exists(), "Align option should work"

    def test_main_quiet_mode(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that main() respects quiet mode."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "quiet_mean.dat"

        test_args = ["mean.py", "-q", "-o", str(output_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Quiet mode should still create output"

    def test_main_invalid_index_file_raises(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that main() raises error for invalid index file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Create invalid index file (non-integer content)
        index_file = temp_output_dir / "bad_index.txt"
        index_file.write_text("not an integer")

        test_args = ["mean.py", "-i", str(index_file), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            with pytest.raises(RuntimeError, match="index file must be a space-seperated list"):
                main()


# =============================================================================
# Mean Calculation Validation Tests
# =============================================================================

class TestMeanCalculationValidation:
    """Tests that verify the correctness of the mean calculation."""

    def test_mean_calculation_accuracy(self, trajectory_info):
        """
        Test that mean calculation produces correct arithmetic mean of aligned positions.

        This test manually computes the expected mean for a small number of configurations
        and compares it to the mean() function output.
        """
        top_info, traj_info = trajectory_info

        # Use a small, deterministic subset of configurations
        traj_info_subset = TrajInfo(
            traj_info.path,
            nconfs=3,
            idxs=traj_info.idxs[0:3],  # First 3 configurations
            incl_v=traj_info.incl_v
        )

        # Use first configuration as reference for alignment
        ref_conf = get_confs(top_info, traj_info_subset, 0, 1)[0]
        ref_conf = inbox(ref_conf, center=True)

        # Compute mean using the function
        result = mean(traj_info_subset, top_info, ref_conf=ref_conf, ncpus=1)

        # Manually compute expected mean
        # Get all configurations
        confs = get_confs(top_info, traj_info_subset, 0, traj_info_subset.nconfs)

        # Apply same processing as mean() function
        confs = [inbox(c, center=True) for c in confs]

        # Align to reference
        indexes = np.arange(top_info.nbases)
        reference_coords = ref_conf.positions
        ref_cms = np.mean(reference_coords, axis=0)
        centered_ref_coords = reference_coords - ref_cms

        # Manually align each configuration and accumulate
        manual_sum = np.zeros([3, top_info.nbases, 3])
        for conf in confs:
            conf_array = np.array([conf.positions, conf.a1s, conf.a3s])
            aligned = svd_align(centered_ref_coords, conf_array, indexes, ref_center=np.zeros(3))
            manual_sum += aligned

        # Compute manual mean
        manual_mean = manual_sum / traj_info_subset.nconfs
        manual_positions = manual_mean[0]

        # Compare positions (main validation)
        np.testing.assert_allclose(
            result.positions,
            manual_positions,
            rtol=1e-6,
            err_msg="Mean positions don't match manually computed mean"
        )

        # Verify the mean is actually different from any single configuration
        # (sanity check that we're not just returning the reference)
        for i, conf in enumerate(confs):
            position_diff = np.linalg.norm(result.positions - conf.positions)
            assert position_diff > 0.01, \
                f"Mean should differ from individual configuration {i}"

    def test_mean_accumulation_across_chunks(self, trajectory_info):
        """
        Test that mean accumulation works correctly when data is split across chunks.

        This verifies that the callback function correctly accumulates results
        from multiple compute() calls.
        """
        top_info, traj_info = trajectory_info

        # Use a small deterministic subset
        traj_info_subset = TrajInfo(
            traj_info.path,
            nconfs=4,
            idxs=traj_info.idxs[0:4],
            incl_v=traj_info.incl_v
        )

        # Fixed reference for reproducibility
        ref_conf = get_confs(top_info, traj_info_subset, 0, 1)[0]
        ref_conf = inbox(ref_conf, center=True)

        # Prepare context for manual chunk processing
        indexes = list(range(top_info.nbases))
        reference_coords = ref_conf.positions[indexes]
        ref_cms = np.mean(reference_coords, axis=0)
        centered_ref_coords = reference_coords - ref_cms

        ctx = ComputeContext(
            traj_info=traj_info_subset,
            top_info=top_info,
            centered_ref_coords=centered_ref_coords,
            indexes=indexes
        )

        # Manually process two chunks (2 confs each)
        chunk_size = 2
        chunk_0_result = compute(ctx, chunk_size=chunk_size, chunk_id=0)
        chunk_1_result = compute(ctx, chunk_size=chunk_size, chunk_id=1)

        # Manually accumulate
        manual_accumulated = chunk_0_result + chunk_1_result
        manual_mean = manual_accumulated / traj_info_subset.nconfs

        # Compare with mean() function result
        result = mean(traj_info_subset, top_info, ref_conf=ref_conf, ncpus=1)

        np.testing.assert_allclose(
            result.positions,
            manual_mean[0],  # positions are first element
            rtol=1e-6,
            err_msg="Accumulated mean across chunks doesn't match mean() output"
        )


# =============================================================================
# Parallel Processing Consistency Tests
# =============================================================================

class TestParallelConsistency:
    """Tests that verify parallel processing produces same results as serial."""

    def test_parallel_vs_serial_consistency(self, trajectory_info):
        """Test that ncpus=1 and ncpus=2 produce identical results.

        This ensures the parallelization doesn't introduce any numerical
        differences or race conditions.
        """
        top_info, traj_info = trajectory_info

        # Use fixed reference for reproducibility
        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]
        ref_conf = inbox(ref_conf, center=True)

        # Compute with serial processing
        result_serial = mean(traj_info, top_info, ref_conf=ref_conf, ncpus=1)

        # Compute with parallel processing
        result_parallel = mean(traj_info, top_info, ref_conf=ref_conf, ncpus=2)

        # Results should be identical (within floating point tolerance)
        np.testing.assert_allclose(
            result_serial.positions,
            result_parallel.positions,
            rtol=1e-10,
            atol=1e-12,
            err_msg="Parallel and serial processing should produce identical positions"
        )

        np.testing.assert_allclose(
            result_serial.a1s,
            result_parallel.a1s,
            rtol=1e-10,
            atol=1e-12,
            err_msg="Parallel and serial processing should produce identical a1 vectors"
        )

        np.testing.assert_allclose(
            result_serial.a3s,
            result_parallel.a3s,
            rtol=1e-10,
            atol=1e-12,
            err_msg="Parallel and serial processing should produce identical a3 vectors"
        )

    def test_parallel_consistency_with_subset_indexes(self, trajectory_info):
        """Test parallel consistency when using subset of indexes for alignment."""
        top_info, traj_info = trajectory_info

        # Use first half of particles for alignment
        indexes = list(range(top_info.nbases // 2))

        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]
        ref_conf = inbox(ref_conf, center=True)

        result_serial = mean(traj_info, top_info, ref_conf=ref_conf, indexes=indexes, ncpus=1)
        result_parallel = mean(traj_info, top_info, ref_conf=ref_conf, indexes=indexes, ncpus=2)

        np.testing.assert_allclose(
            result_serial.positions,
            result_parallel.positions,
            rtol=1e-10,
            atol=1e-12,
            err_msg="Parallel consistency should hold with subset indexes"
        )


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_mean_single_conf_and_high_ncpus(self, temp_output_dir, trajectory_info):
        """Test mean() works with single configuration and handles high ncpus gracefully."""
        top_info, traj_info = trajectory_info

        # Create a single-conf trajectory
        single_conf = get_confs(top_info, traj_info, 0, 1)[0]
        single_traj = temp_output_dir / "single.dat"
        write_conf(str(single_traj), single_conf)

        single_top, single_traj_info = describe(None, str(single_traj))

        # Test single conf
        result = mean(single_traj_info, single_top, ncpus=1)
        assert result.positions.shape == single_conf.positions.shape, \
            "Mean of single conf should match shape"

        # Test with more CPUs than configurations
        result_high_cpu = mean(traj_info, top_info, ncpus=100)
        assert isinstance(result_high_cpu, Configuration), \
            "Should handle ncpus > nconfs gracefully"
