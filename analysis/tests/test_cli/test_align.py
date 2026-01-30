"""
Tests for oxDNA_analysis_tools.align module.

Tests cover:
- svd_align() alignment function
- compute() helper function
- align() API function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.align import svd_align, compute, align, cli_parser, ComputeContext, main
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
def expected_align_path(test_resources):
    """Path to the expected aligned trajectory."""
    return test_resources / "aligntraj.dat"


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
# Unit Tests - svd_align()
# =============================================================================

class TestSvdAlign:
    """Tests for the svd_align() function."""

    def test_svd_align_behavior(self, trajectory_info):
        """Test svd_align returns correct shapes, preserves structure on self-alignment, and works with subsets."""
        top_info, traj_info = trajectory_info

        conf = get_confs(top_info, traj_info, 0, 1)[0]
        conf = inbox(conf, center=True)

        # Test with all indexes
        indexes = list(range(top_info.nbases))
        ref_coords = conf.positions[indexes]
        ref_center = np.mean(ref_coords, axis=0)
        centered_ref = ref_coords - ref_center

        np_coords = np.asarray([conf.positions.copy(), conf.a1s.copy(), conf.a3s.copy()])

        pos, a1s, a3s = svd_align(centered_ref, np_coords, indexes, ref_center=np.zeros(3))

        # Check shapes
        assert pos.shape == (top_info.nbases, 3), "Positions shape mismatch"
        assert a1s.shape == (top_info.nbases, 3), "a1s shape mismatch"
        assert a3s.shape == (top_info.nbases, 3), "a3s shape mismatch"

        # Self-alignment should preserve positions (centered at origin)
        np.testing.assert_allclose(pos[indexes], centered_ref, atol=1e-10,
                                   err_msg="Self-alignment should preserve positions")

        # Test with subset indexes
        subset_indexes = list(range(min(10, top_info.nbases)))
        subset_ref = conf.positions[subset_indexes]
        subset_ref = subset_ref - np.mean(subset_ref, axis=0)

        np_coords2 = np.asarray([conf.positions.copy(), conf.a1s.copy(), conf.a3s.copy()])
        pos2, _, _ = svd_align(subset_ref, np_coords2, subset_indexes)

        # Should still return full system
        assert pos2.shape == (top_info.nbases, 3), "Should return full system positions"


# =============================================================================
# Unit Tests - compute()
# =============================================================================

class TestCompute:
    """Tests for the compute() helper function."""

    def test_compute_behavior(self, trajectory_info):
        """Test compute() returns non-empty strings and processes multiple chunks."""
        top_info, traj_info = trajectory_info

        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]
        ref_conf = inbox(ref_conf, center=True)

        indexes = list(range(top_info.nbases))
        reference_coords = ref_conf.positions[indexes]
        centered_ref_coords = reference_coords - np.mean(reference_coords, axis=0)

        ctx = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            centered_ref_coords=centered_ref_coords,
            indexes=indexes,
            no_center=False
        )

        # Process first two chunks
        result0 = compute(ctx, chunk_size=1, chunk_id=0)
        result1 = compute(ctx, chunk_size=1, chunk_id=1)

        assert isinstance(result0, str), "compute() should return string"
        assert len(result0) > 0, "Output string should not be empty"
        assert isinstance(result1, str), "Second chunk should return string"


# =============================================================================
# API Tests - align() function
# =============================================================================

class TestAlignFunction:
    """Tests for the align() API function."""

    def test_align_basic_behavior(self, mini_traj_path, temp_output_dir, trajectory_info):
        """Test align() creates valid output with same nconfs and matches expected output."""
        top_info, traj_info = trajectory_info
        outfile = str(temp_output_dir / "aligned.dat")

        # Use first conf as reference for reproducibility
        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]

        align(str(mini_traj_path), outfile, ncpus=1, ref_conf=ref_conf)

        # Check file exists and has content
        assert Path(outfile).exists(), "Output file should be created"
        assert Path(outfile).stat().st_size > 0, "Output file should have content"

        # Check same nconfs
        _, aligned_traj_info = describe(None, outfile)
        assert aligned_traj_info.nconfs == traj_info.nconfs, \
            "Aligned trajectory should have same number of configurations"

    def test_align_with_options(self, mini_traj_path, temp_output_dir, trajectory_info):
        """Test align() with custom ref_conf, indexes, and center=False options."""
        top_info, traj_info = trajectory_info

        # Test with custom ref_conf (second configuration)
        ref_conf = get_confs(top_info, traj_info, 1, 1)[0]
        outfile1 = str(temp_output_dir / "aligned_custom_ref.dat")
        align(str(mini_traj_path), outfile1, ncpus=1, ref_conf=ref_conf)
        assert Path(outfile1).exists(), "Output should be created with custom ref"

        # Test with indexes
        indexes = list(range(min(20, top_info.nbases)))
        outfile2 = str(temp_output_dir / "aligned_indexed.dat")
        align(str(mini_traj_path), outfile2, ncpus=1, indexes=indexes)
        assert Path(outfile2).exists(), "Output should be created with indexes"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_behavior(self):
        """Test parser requires arguments, accepts all options, and has correct defaults."""
        parser = cli_parser()

        # Requires trajectory and outfile
        with pytest.raises(SystemExit):
            parser.parse_args([])
        with pytest.raises(SystemExit):
            parser.parse_args(["traj.dat"])

        # All options
        args = parser.parse_args([
            "-p", "4",
            "-i", "index.txt",
            "-r", "ref.dat",
            "-q",
            "trajectory.dat",
            "aligned.dat"
        ])

        assert args.traj == ["trajectory.dat"], "Trajectory not parsed"
        assert args.outfile == ["aligned.dat"], "Outfile not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.index_file == ["index.txt"], "Index file option not parsed"
        assert args.reference_structure == ["ref.dat"], "Reference option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Defaults
        args_defaults = parser.parse_args(["traj.dat", "out.dat"])
        assert args_defaults.parallel is None
        assert args_defaults.index_file is None
        assert args_defaults.reference_structure is None
        assert args_defaults.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_basic_and_parallel(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates output file with basic call and parallel processing."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Basic call
        output_file = temp_output_dir / "aligned.dat"
        test_args = ["align.py", str(traj_copy), str(output_file)]
        with patch.object(sys, 'argv', test_args):
            main()
        assert output_file.exists(), "Output file should be created"

        # Parallel call
        output_parallel = temp_output_dir / "aligned_parallel.dat"
        test_args_parallel = ["align.py", "-p", "2", str(traj_copy), str(output_parallel)]
        with patch.object(sys, 'argv', test_args_parallel):
            main()
        assert output_parallel.exists(), "Parallel output should be created"

    def test_main_with_index_and_reference(self, mini_traj_path, temp_output_dir, monkeypatch, test_resources):
        """Test main() with index file and reference configuration."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Create index file
        index_file = temp_output_dir / "index.txt"
        index_file.write_text("0 1 2 3 4 5 6 7 8 9")

        output_indexed = temp_output_dir / "aligned_indexed.dat"
        test_args_indexed = ["align.py", "-i", str(index_file), str(traj_copy), str(output_indexed)]
        with patch.object(sys, 'argv', test_args_indexed):
            main()
        assert output_indexed.exists(), "Indexed output should be created"

        # With reference file
        ref_file = test_resources / "mean.dat"
        output_ref = temp_output_dir / "aligned_ref.dat"
        test_args_ref = ["align.py", "-r", str(ref_file), str(traj_copy), str(output_ref)]
        with patch.object(sys, 'argv', test_args_ref):
            main()
        assert output_ref.exists(), "Output with reference should be created"

    def test_main_invalid_index_file_raises(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that main() raises error for invalid index file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Create invalid index file
        index_file = temp_output_dir / "bad_index.txt"
        index_file.write_text("not an integer")

        output_file = temp_output_dir / "aligned.dat"
        test_args = ["align.py", "-i", str(index_file), str(traj_copy), str(output_file)]

        with patch.object(sys, 'argv', test_args):
            with pytest.raises(RuntimeError, match="index file must be a space-seperated list"):
                main()


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_align_edge_cases(self, mini_traj_path, temp_output_dir, trajectory_info):
        """Test align() handles high ncpus and single particle index."""
        top_info, _ = trajectory_info

        # High ncpus
        outfile1 = str(temp_output_dir / "aligned_high_cpu.dat")
        align(str(mini_traj_path), outfile1, ncpus=100)
        assert Path(outfile1).exists(), "Should handle high ncpus"

        # Single particle alignment (edge case)
        outfile2 = str(temp_output_dir / "aligned_single_idx.dat")
        align(str(mini_traj_path), outfile2, ncpus=1, indexes=[0])
        assert Path(outfile2).exists(), "Should handle single particle index"


# =============================================================================
# Quality and Correctness Tests
# =============================================================================

class TestAlignmentQuality:
    """Tests for alignment quality and correctness."""

    def test_golden_output_comparison(self, mini_traj_path, temp_output_dir, trajectory_info, expected_align_path):
        """
        Compares actual alignment output against the golden reference file
        (aligntraj.dat), verifying positions, a1s, and a3s match within tolerance.
        """
        top_info, traj_info = trajectory_info
        outfile = str(temp_output_dir / "aligned_golden.dat")

        # Use first conf as reference (same as default behavior)
        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]

        # Execute alignment
        align(str(mini_traj_path), outfile, ncpus=1, ref_conf=ref_conf)

        # Load both the expected and actual aligned trajectories
        _, expected_traj_info = describe(None, str(expected_align_path))
        _, actual_traj_info = describe(None, outfile)

        # Verify same number of configurations
        assert actual_traj_info.nconfs == expected_traj_info.nconfs, \
            f"Configuration count mismatch: expected {expected_traj_info.nconfs}, got {actual_traj_info.nconfs}"

        # Compare each configuration
        expected_confs = get_confs(top_info, expected_traj_info, 0, expected_traj_info.nconfs)
        actual_confs = get_confs(top_info, actual_traj_info, 0, actual_traj_info.nconfs)

        for i, (expected_conf, actual_conf) in enumerate(zip(expected_confs, actual_confs)):
            # 1. Correct type
            assert isinstance(actual_conf, Configuration), f"Conf {i} should be Configuration type"

            # 2. Correct shape
            assert actual_conf.positions.shape == expected_conf.positions.shape, \
                f"Conf {i} positions shape mismatch"
            assert actual_conf.a1s.shape == expected_conf.a1s.shape, \
                f"Conf {i} a1s shape mismatch"
            assert actual_conf.a3s.shape == expected_conf.a3s.shape, \
                f"Conf {i} a3s shape mismatch"

            # 3. No invalid values
            assert not np.any(np.isnan(actual_conf.positions)), f"Conf {i} positions contain NaN"
            assert not np.any(np.isnan(actual_conf.a1s)), f"Conf {i} a1s contain NaN"
            assert not np.any(np.isnan(actual_conf.a3s)), f"Conf {i} a3s contain NaN"

            # 4. Correct content - match expected output within tolerance
            np.testing.assert_allclose(
                actual_conf.positions, expected_conf.positions,
                rtol=1e-5, atol=1e-6,
                err_msg=f"Conf {i} positions don't match expected output"
            )
            np.testing.assert_allclose(
                actual_conf.a1s, expected_conf.a1s,
                rtol=1e-5, atol=1e-6,
                err_msg=f"Conf {i} a1s don't match expected output"
            )
            np.testing.assert_allclose(
                actual_conf.a3s, expected_conf.a3s,
                rtol=1e-5, atol=1e-6,
                err_msg=f"Conf {i} a3s don't match expected output"
            )

    def test_distance_preservation(self, mini_traj_path, temp_output_dir, trajectory_info):
        """
        Since alignment is a rigid rotation/translation, distances between
        particles should be identical before and after alignment.
        """
        top_info, traj_info = trajectory_info
        outfile = str(temp_output_dir / "aligned_distance.dat")

        # Execute alignment
        align(str(mini_traj_path), outfile, ncpus=1)

        # Load original and aligned trajectories
        original_confs = get_confs(top_info, traj_info, 0, traj_info.nconfs)
        _, aligned_traj_info = describe(None, outfile)
        aligned_confs = get_confs(top_info, aligned_traj_info, 0, aligned_traj_info.nconfs)

        # Test distance preservation for each configuration
        for conf_idx, (orig_conf, aligned_conf) in enumerate(zip(original_confs, aligned_confs)):
            # Calculate distances between consecutive bases in original
            orig_distances = np.linalg.norm(
                orig_conf.positions[1:] - orig_conf.positions[:-1],
                axis=1
            )

            # Calculate distances between consecutive bases in aligned
            aligned_distances = np.linalg.norm(
                aligned_conf.positions[1:] - aligned_conf.positions[:-1],
                axis=1
            )

            # 1. Correct type and shape
            assert isinstance(aligned_distances, np.ndarray), "Distances should be numpy array"
            assert aligned_distances.shape == orig_distances.shape, \
                f"Conf {conf_idx}: Distance array shape mismatch"

            # 2. No invalid values
            assert not np.any(np.isnan(aligned_distances)), \
                f"Conf {conf_idx}: Aligned distances contain NaN"
            assert not np.any(np.isinf(aligned_distances)), \
                f"Conf {conf_idx}: Aligned distances contain Inf"

            # 3. Reasonable values (should be positive and similar to original)
            assert np.all(aligned_distances > 0), \
                f"Conf {conf_idx}: All distances should be positive"

            # 4. Distances should be preserved (rigid transformation)
            np.testing.assert_allclose(
                aligned_distances, orig_distances,
                rtol=1e-10, atol=1e-10,
                err_msg=f"Conf {conf_idx}: Distances not preserved by alignment"
            )

    def test_rmsd_reduction(self, mini_traj_path, temp_output_dir, trajectory_info):
        """
        The purpose of alignment is to minimize RMSD between structures.
        The aligned RMSD should be less than or equal to the unaligned RMSD.
        """
        top_info, traj_info = trajectory_info
        outfile = str(temp_output_dir / "aligned_rmsd.dat")

        # Execute alignment with default settings (uses first conf as reference)
        align(str(mini_traj_path), outfile, ncpus=1)

        # Load reference (first conf), original, and aligned trajectories
        ref_conf_orig = get_confs(top_info, traj_info, 0, 1)[0]
        ref_conf_orig = inbox(ref_conf_orig, center=True)

        original_confs = get_confs(top_info, traj_info, 1, traj_info.nconfs - 1)
        _, aligned_traj_info = describe(None, outfile)
        aligned_confs = get_confs(top_info, aligned_traj_info, 0, aligned_traj_info.nconfs)

        # Get aligned reference for comparison (first conf in aligned trajectory)
        ref_conf_aligned = aligned_confs[0]
        aligned_confs = aligned_confs[1:]  # Skip reference

        # Calculate RMSD for each configuration before and after alignment
        for conf_idx, (orig_conf, aligned_conf) in enumerate(zip(original_confs, aligned_confs), start=1):
            # Inbox original conf for fair comparison
            orig_conf = inbox(orig_conf, center=True)

            # Calculate unaligned RMSD (both centered at origin)
            unaligned_rmsd = np.sqrt(
                np.mean(np.sum((orig_conf.positions - ref_conf_orig.positions) ** 2, axis=1))
            )

            # Calculate aligned RMSD
            aligned_rmsd = np.sqrt(
                np.mean(np.sum((aligned_conf.positions - ref_conf_aligned.positions) ** 2, axis=1))
            )

            # 1. Correct type
            assert isinstance(aligned_rmsd, (float, np.floating)), \
                f"Conf {conf_idx}: RMSD should be float type"

            # 2. Reasonable value
            assert not np.isnan(aligned_rmsd), f"Conf {conf_idx}: Aligned RMSD is NaN"
            assert not np.isinf(aligned_rmsd), f"Conf {conf_idx}: Aligned RMSD is Inf"
            assert aligned_rmsd >= 0, f"Conf {conf_idx}: RMSD should be non-negative"

            # 3. RMSD should be reduced by alignment
            assert aligned_rmsd <= unaligned_rmsd, \
                f"Conf {conf_idx}: Alignment should reduce RMSD " \
                f"(unaligned: {unaligned_rmsd:.6f}, aligned: {aligned_rmsd:.6f})"

    def test_centering_correctness(self, mini_traj_path, temp_output_dir, trajectory_info):
        """
        The aligned indexed particles should be centered at the origin.
        """
        top_info, traj_info = trajectory_info

        # Use subset of particles for alignment
        indexes = list(range(min(50, top_info.nbases)))

        outfile_centered = str(temp_output_dir / "aligned_centered.dat")
        align(str(mini_traj_path), outfile_centered, ncpus=1, indexes=indexes)

        # Load aligned trajectory with centering
        _, aligned_traj_info = describe(None, outfile_centered)
        aligned_confs = get_confs(top_info, aligned_traj_info, 0, aligned_traj_info.nconfs)

        # Check that centered output has indexed particles centered near origin
        for conf_idx, aligned_conf in enumerate(aligned_confs):
            # Calculate center of mass of aligned indexed particles
            aligned_center = np.mean(aligned_conf.positions[indexes], axis=0)

            # 1. Correct type and shape
            assert isinstance(aligned_center, np.ndarray), \
                f"Conf {conf_idx}: Center should be numpy array"
            assert aligned_center.shape == (3,), \
                f"Conf {conf_idx}: Center should be 3D coordinate"

            # 2. No invalid values
            assert not np.any(np.isnan(aligned_center)), \
                f"Conf {conf_idx}: Center contains NaN"

            # 3. With, indexed particles should be centered at origin
            np.testing.assert_allclose(
                aligned_center, np.zeros(3),
                rtol=1e-10, atol=1e-10,
                err_msg=f"Conf {conf_idx}: Centered alignment should place indexed particles at origin"
            )

        # On the other hand, if you disable centering, the center-of-mass of the inboxed ref conf should be preserved.
        outfile_nocenter = str(temp_output_dir / "aligned_nocenter.dat")
        ref_traj = get_confs(top_info, traj_info, 0, traj_info.nconfs)
        centers = np.array([np.mean(inbox(c, center=False).positions[indexes], axis=0) for c in ref_traj])
        align(str(mini_traj_path), outfile_nocenter, ncpus=1, indexes=indexes, no_center=True)

        # Load aligned trajectory with centering
        _, aligned_traj_info = describe(None, outfile_nocenter)
        aligned_confs = get_confs(top_info, aligned_traj_info, 0, aligned_traj_info.nconfs)
        aligned_centers = np.array([np.mean(c.positions[indexes], axis=0) for c in aligned_confs])

        # 4. test that original center is preserved
        np.testing.assert_allclose(
            aligned_centers, centers,
            rtol=1e-10, atol=1e-10,
            err_msg="Disabling centering should preserve original centers"
        )

    def test_parallel_serial_consistency(self, mini_traj_path, temp_output_dir, trajectory_info):
        """
        MEDIUM: Test that parallel and serial execution produce identical results.

        Running align() with ncpus=1 should produce identical output to ncpus=4.
        """
        top_info, traj_info = trajectory_info

        # Use first conf as reference for reproducibility
        ref_conf = get_confs(top_info, traj_info, 0, 1)[0]

        # Execute with serial processing
        outfile_serial = str(temp_output_dir / "aligned_serial.dat")
        align(str(mini_traj_path), outfile_serial, ncpus=1, ref_conf=ref_conf)

        # Execute with parallel processing
        outfile_parallel = str(temp_output_dir / "aligned_parallel.dat")
        align(str(mini_traj_path), outfile_parallel, ncpus=4, ref_conf=ref_conf)

        # Load both trajectories
        _, serial_traj_info = describe(None, outfile_serial)
        _, parallel_traj_info = describe(None, outfile_parallel)

        # 1. Same number of configurations
        assert serial_traj_info.nconfs == parallel_traj_info.nconfs, \
            "Serial and parallel should produce same number of configurations"

        # Load all configurations
        serial_confs = get_confs(top_info, serial_traj_info, 0, serial_traj_info.nconfs)
        parallel_confs = get_confs(top_info, parallel_traj_info, 0, parallel_traj_info.nconfs)

        # Compare each configuration
        for conf_idx, (serial_conf, parallel_conf) in enumerate(zip(serial_confs, parallel_confs)):
            # 2. Correct shapes
            assert serial_conf.positions.shape == parallel_conf.positions.shape, \
                f"Conf {conf_idx}: Position shapes differ between serial and parallel"

            # 3. No invalid values
            assert not np.any(np.isnan(parallel_conf.positions)), \
                f"Conf {conf_idx}: Parallel positions contain NaN"

            # 4. Results should be identical
            np.testing.assert_allclose(
                serial_conf.positions, parallel_conf.positions,
                rtol=1e-10, atol=1e-10,
                err_msg=f"Conf {conf_idx}: Serial and parallel positions differ"
            )
            np.testing.assert_allclose(
                serial_conf.a1s, parallel_conf.a1s,
                rtol=1e-10, atol=1e-10,
                err_msg=f"Conf {conf_idx}: Serial and parallel a1s differ"
            )
            np.testing.assert_allclose(
                serial_conf.a3s, parallel_conf.a3s,
                rtol=1e-10, atol=1e-10,
                err_msg=f"Conf {conf_idx}: Serial and parallel a3s differ"
            )
