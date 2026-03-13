"""
Tests for oxDNA_analysis_tools.minify module.

Tests cover:
- minify() API function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.minify import minify, cli_parser, main, compute, ComputeContext
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, _is_zstd_traj, inbox
from oxDNA_analysis_tools.mean import mean


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
# API Tests - minify() function
# =============================================================================

class TestMinifyFunction:
    """Tests for the minify() API function."""

    def test_minify_behavior(self, trajectory_info, mini_traj_path, temp_output_dir):
        """Test minify preserves structure, handles decimals, discards a-vectors, and preserves positions."""
        top_info, traj_info = trajectory_info

        # Basic minify - preserves structure
        output_basic = temp_output_dir / "minified.dat"
        minify(traj_info, top_info, str(output_basic), d=None, a=False, ncpus=1)

        assert output_basic.exists(), "Output file should be created"
        min_top, min_traj = describe(None, str(output_basic))
        assert min_traj.nconfs == traj_info.nconfs, "Should have same number of confs"
        assert min_top.nbases == top_info.nbases, "Should have same number of particles"

        # Compare positions are preserved
        orig_conf = get_confs(top_info, traj_info, 0, 1)[0]
        min_conf = get_confs(min_top, min_traj, 0, 1)[0]
        np.testing.assert_allclose(min_conf.positions, orig_conf.positions, rtol=1e-5,
                                   err_msg="Positions should be preserved")

        # With decimal rounding
        output_rounded = temp_output_dir / "minified_rounded.dat"
        minify(traj_info, top_info, str(output_rounded), d=2, a=False, ncpus=1)

        assert output_rounded.exists(), "Rounded output file should be created"
        round_top, round_traj = describe(None, str(output_rounded))
        assert round_traj.nconfs == traj_info.nconfs, "Should have same confs after rounding"

        # With a=True to discard orientation vectors
        output_no_a = temp_output_dir / "minified_no_a.dat"
        minify(traj_info, top_info, str(output_no_a), d=None, a=True, ncpus=1)

        assert output_no_a.exists(), "No-a output file should be created"
        no_a_top, no_a_traj = describe(None, str(output_no_a))
        no_a_conf = get_confs(no_a_top, no_a_traj, 0, 1)[0]

        # a1s and a3s should all be zero
        assert np.allclose(no_a_conf.a1s, 0), "a1 vectors should be zero"
        assert np.allclose(no_a_conf.a3s, 0), "a3 vectors should be zero"

        # Directly test the compute function (uncompressed path)
        ctx = ComputeContext(traj_info, top_info, d=2, a=True, compress=False)
        min_conf_str = compute(ctx, chunk_size=1, chunk_id=0)

        # Check that compute returns bytes
        assert isinstance(min_conf_str, bytes), "compute() should return bytes"

# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_arguments(self):
        """Test that parser requires trajectory and outfile arguments."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-p", "4",
            "-a",
            "-d", "3",
            "-q",
            "input.dat",
            "output.dat"
        ])

        assert args.trajectory == "input.dat", "Trajectory not parsed"
        assert args.outfile == "output.dat", "Outfile not parsed"
        assert args.parallel == 4, "Parallel option not parsed"
        assert args.no_a is True, "No-a option not parsed"
        assert args.decimals == 3, "Decimals option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["input.dat", "output.dat"])
        assert args_defaults.no_a is False, "No-a should default to False"
        assert args_defaults.decimals is None, "Decimals should default to None"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_behavior(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates output with basic call, decimals option, and -a flag."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Basic call
        output_basic = temp_output_dir / "minified_main.dat"
        test_args = ["minify.py", str(traj_copy), str(output_basic)]
        with patch.object(sys, 'argv', test_args):
            main()
        assert output_basic.exists(), "Output file should be created"

        # With decimals option
        output_dec = temp_output_dir / "minified_dec.dat"
        test_args_dec = ["minify.py", "-d", "2", str(traj_copy), str(output_dec)]
        with patch.object(sys, 'argv', test_args_dec):
            main()
        assert output_dec.exists(), "Output file should be created with decimals"

        # With -a flag
        output_no_a = temp_output_dir / "minified_no_a.dat"
        test_args_no_a = ["minify.py", "-a", str(traj_copy), str(output_no_a)]
        with patch.object(sys, 'argv', test_args_no_a):
            main()
        assert output_no_a.exists(), "Output file should be created with -a flag"

        # With -z flag
        output_z = temp_output_dir / "minified_z.dat"
        test_args_z = ["minify.py", "-z", str(traj_copy), str(output_z)]
        with patch.object(sys, 'argv', test_args_z):
            main()
        assert output_z.exists(), "Output file should be created with -z flag"
        assert _is_zstd_traj(str(output_z)), "Output should be zstd-compressed with -z flag"


# =============================================================================
# Compressed output / mean consistency Tests
# =============================================================================

class TestCompressedMeanConsistency:
    """Test that mean() produces equivalent results on a minify-compressed trajectory."""

    def test_mean_compressed_matches_original(self, trajectory_info, temp_output_dir):
        """
        Compress the test trajectory with minify(-z), run mean() on both the
        original and the compressed copy, and verify the results are identical.

        Validation is progressive:
          1. The compressed output is a valid zstd trajectory with the right shape.
          2. mean() on the compressed trajectory returns a Configuration with
             the correct shape and normalized orientation vectors.
          3. Positions, a1s, and a3s from both mean() calls are numerically
             identical (no rounding is applied, so values should be bit-exact
             up to floating-point round-trip through text).
        """
        top_info, traj_info = trajectory_info

        # --- Step 1: compress the trajectory with minify ---
        compressed_path = str(temp_output_dir / "compressed.dat")
        minify(traj_info, top_info, compressed_path, compress=True, ncpus=1)

        assert _is_zstd_traj(compressed_path), "minify output should be zstd-compressed"
        compressed_top, compressed_traj = describe(None, compressed_path)
        assert compressed_traj.nconfs == traj_info.nconfs, \
            "Compressed trajectory should have the same number of configurations"
        assert compressed_top.nbases == top_info.nbases, \
            "Compressed trajectory should have the same number of particles"

        # --- Step 2: run mean() on both trajectories with the same reference ---
        ref_conf = inbox(get_confs(top_info, traj_info, 0, 1)[0], center=True)

        mean_orig = mean(traj_info, top_info, ref_conf=ref_conf, ncpus=1)
        mean_comp = mean(compressed_traj, compressed_top, ref_conf=ref_conf, ncpus=1)

        assert mean_comp.positions.shape == (top_info.nbases, 3), \
            "Mean of compressed trajectory has wrong positions shape"
        assert mean_comp.a1s.shape == (top_info.nbases, 3), \
            "Mean of compressed trajectory has wrong a1s shape"
        assert mean_comp.a3s.shape == (top_info.nbases, 3), \
            "Mean of compressed trajectory has wrong a3s shape"

        a1_norms = np.linalg.norm(mean_comp.a1s, axis=1)
        a3_norms = np.linalg.norm(mean_comp.a3s, axis=1)
        np.testing.assert_allclose(a1_norms, np.ones(top_info.nbases), rtol=1e-5,
                                   err_msg="a1 vectors from compressed mean should be normalized")
        np.testing.assert_allclose(a3_norms, np.ones(top_info.nbases), rtol=1e-5,
                                   err_msg="a3 vectors from compressed mean should be normalized")

        # --- Step 3: positions and orientations should be numerically identical ---
        np.testing.assert_allclose(mean_comp.positions, mean_orig.positions, rtol=1e-5,
                                   err_msg="Mean positions differ between original and compressed")
        np.testing.assert_allclose(mean_comp.a1s, mean_orig.a1s, rtol=1e-5,
                                   err_msg="Mean a1 vectors differ between original and compressed")
        np.testing.assert_allclose(mean_comp.a3s, mean_orig.a3s, rtol=1e-5,
                                   err_msg="Mean a3 vectors differ between original and compressed")
