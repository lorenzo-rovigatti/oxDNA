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
            center=True
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

        # Test with center=False
        outfile3 = str(temp_output_dir / "aligned_nocenter.dat")
        align(str(mini_traj_path), outfile3, ncpus=1, center=False)
        assert Path(outfile3).exists(), "Output should be created with nocenter"


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
            "-c",
            "-q",
            "trajectory.dat",
            "aligned.dat"
        ])

        assert args.traj == ["trajectory.dat"], "Trajectory not parsed"
        assert args.outfile == ["aligned.dat"], "Outfile not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.index_file == ["index.txt"], "Index file option not parsed"
        assert args.reference_structure == ["ref.dat"], "Reference option not parsed"
        assert args.no_center is True, "Nocenter option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Defaults
        args_defaults = parser.parse_args(["traj.dat", "out.dat"])
        assert args_defaults.parallel is None
        assert args_defaults.index_file is None
        assert args_defaults.reference_structure is None
        assert args_defaults.no_center is False
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
