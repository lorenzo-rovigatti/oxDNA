"""
Tests for oxDNA_analysis_tools.distance module.

Tests cover:
- min_image() minimum image distance function
- vectorized_min_image() vectorized minimum image distance
- compute() helper function
- distance() API function
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

from oxDNA_analysis_tools.distance import (
    min_image,
    vectorized_min_image,
    compute,
    distance,
    cli_parser,
    ComputeContext,
    main
)
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs


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
# Unit Tests - min_image()
# =============================================================================

class TestMinImage:
    """Tests for the min_image() minimum image distance function."""

    def test_min_image_basic_distances(self):
        """Test basic distance calculations: same point, simple distance, and type."""
        box = 10.0

        # Same point should give zero
        p_same = np.array([5.0, 5.0, 5.0])
        assert min_image(p_same, p_same, box) == pytest.approx(0.0), "Same point should be zero"

        # Simple 3-4-5 triangle
        p1 = np.array([0.0, 0.0, 0.0])
        p2 = np.array([3.0, 4.0, 0.0])
        result = min_image(p1, p2, 100.0)  # Large box, no PBC effect
        assert result == pytest.approx(5.0), "Simple distance should be 5.0"

        # Should return float
        assert isinstance(result, float), "Should return Python float"

    def test_min_image_pbc_wrapping(self):
        """Test periodic boundary conditions in various scenarios."""
        box = 10.0

        # 1D wrap: points near opposite edges
        p1 = np.array([0.5, 0.5, 0.5])
        p2 = np.array([9.5, 0.5, 0.5])
        result_1d = min_image(p1, p2, box)
        assert result_1d == pytest.approx(1.0), "1D PBC wrap should give 1.0"

        # 3D wrap: all dimensions
        p3 = np.array([9.5, 9.5, 9.5])
        result_3d = min_image(p1, p3, box)
        assert result_3d == pytest.approx(np.sqrt(3)), "3D wrap should be sqrt(3)"

        # At box boundary (half-box)
        p4 = np.array([0.0, 0.0, 0.0])
        p5 = np.array([5.0, 0.0, 0.0])
        result_boundary = min_image(p4, p5, box)
        assert result_boundary == pytest.approx(5.0), "Half-box distance"

    def test_min_image_coords_outside_box(self):
        """Test handling of coordinates outside box and negative coordinates."""
        box = 10.0

        # Negative coordinates
        p1 = np.array([-0.5, 0.0, 0.0])
        p2 = np.array([0.5, 0.0, 0.0])
        assert min_image(p1, p2, box) == pytest.approx(1.0), "Negative coords should work"

        # Coordinates outside box (should wrap)
        p3 = np.array([15.0, 0.0, 0.0])  # Equivalent to 5.0
        p4 = np.array([25.0, 0.0, 0.0])  # Equivalent to 5.0
        assert min_image(p3, p4, box) == pytest.approx(0.0), "Outside coords should wrap"


# =============================================================================
# Unit Tests - vectorized_min_image()
# =============================================================================

class TestVectorizedMinImage:
    """Tests for vectorized_min_image() function."""

    def test_vectorized_basic_behavior(self):
        """Test vectorized function returns correct shapes and values."""
        box = 10.0

        # Single point each
        p1s_single = np.array([[5.0, 5.0, 5.0]])
        p2s_single = np.array([[5.0, 5.0, 5.0]])
        result_single = vectorized_min_image(p1s_single, p2s_single, box)
        assert result_single.shape == (1, 1), "Single point shape should be (1,1)"
        assert result_single[0, 0] == pytest.approx(0.0), "Same point should be zero"

        # Multiple points - check shape
        n_p1, n_p2 = 5, 3
        p1s = np.random.rand(n_p1, 3) * 10
        p2s = np.random.rand(n_p2, 3) * 10
        result = vectorized_min_image(p1s, p2s, box)
        assert result.shape == (n_p2, n_p1), f"Shape should be ({n_p2}, {n_p1})"

    def test_vectorized_matches_scalar(self):
        """Test that vectorized matches scalar min_image results."""
        p1s = np.array([
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0]
        ])
        p2s = np.array([
            [2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0]
        ])
        box = 20.0

        vectorized_result = vectorized_min_image(p1s, p2s, box)

        # Compare with scalar version
        for i, p2 in enumerate(p2s):
            for j, p1 in enumerate(p1s):
                scalar_result = min_image(p1, p2, box)
                assert vectorized_result[i, j] == pytest.approx(scalar_result, rel=1e-10), \
                    f"Mismatch at ({i},{j})"

    def test_vectorized_with_pbc(self):
        """Test vectorized function with PBC wrapping."""
        p1s = np.array([[0.5, 0.5, 0.5]])
        p2s = np.array([[9.5, 0.5, 0.5]])
        box = 10.0

        result = vectorized_min_image(p1s, p2s, box)
        assert result[0, 0] == pytest.approx(1.0), "PBC wrap should work"


# =============================================================================
# Unit Tests - compute()
# =============================================================================

class TestCompute:
    """Tests for the compute() helper function."""

    def test_compute_behavior(self, trajectory_info):
        """Test compute() returns correct array type, shape, and positive distances."""
        top_info, traj_info = trajectory_info

        n_pairs = 3
        p1s = [0, 1, 2]
        p2s = [3, 4, 5]

        ctx = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            p1s=p1s,
            p2s=p2s
        )

        chunk_size = 2
        result = compute(ctx, chunk_size=chunk_size, chunk_id=0)

        # Check type and shape
        assert isinstance(result, np.ndarray), "Should return numpy array"
        assert result.shape[0] == n_pairs, f"First dimension should be {n_pairs}"

        # Distances should be positive (scaled by 0.85 nm)
        assert np.all(result >= 0), "All distances should be non-negative"


# =============================================================================
# API Tests - distance() function
# =============================================================================

class TestDistanceFunction:
    """Tests for the distance() API function."""

    def test_distance_basic_behavior(self, trajectory_info):
        """Test distance() returns correct structure with positive values."""
        top_info, traj_info = trajectory_info

        p1ss = [[0]]
        p2ss = [[5]]

        result = distance([traj_info], [top_info], p1ss, p2ss, ncpus=1)

        # Check structure
        assert isinstance(result, list), "Should return list"
        assert len(result) == 1, "One trajectory"
        assert len(result[0]) == 1, "One pair"
        assert len(result[0][0]) == traj_info.nconfs, "One distance per conf"

        # All distances should be positive
        assert all(d >= 0 for d in result[0][0]), "Distances should be non-negative"

    def test_distance_multiple_pairs(self, trajectory_info):
        """Test distance calculation for multiple pairs."""
        top_info, traj_info = trajectory_info

        p1ss = [[0, 1, 2]]
        p2ss = [[5, 6, 7]]

        result = distance([traj_info], [top_info], p1ss, p2ss, ncpus=1)

        assert len(result[0]) == 3, "Should have three pairs"
        for pair_distances in result[0]:
            assert len(pair_distances) == traj_info.nconfs, "Each pair needs nconfs distances"

    def test_distance_same_particle_zero(self, trajectory_info):
        """Test that distance between same particle is zero."""
        top_info, traj_info = trajectory_info

        p1ss = [[0]]
        p2ss = [[0]]

        result = distance([traj_info], [top_info], p1ss, p2ss, ncpus=1)

        assert all(d == pytest.approx(0.0) for d in result[0][0]), "Same particle should be zero"

    def test_distance_symmetry(self, trajectory_info):
        """Test that distance(A, B) == distance(B, A)."""
        top_info, traj_info = trajectory_info

        result_forward = distance([traj_info], [top_info], [[0]], [[5]], ncpus=1)
        result_backward = distance([traj_info], [top_info], [[5]], [[0]], ncpus=1)

        np.testing.assert_allclose(result_forward[0][0], result_backward[0][0],
                                   err_msg="Distance should be symmetric")

    def test_distance_multiple_trajectories(self, trajectory_info):
        """Test distance calculation across multiple trajectories."""
        top_info, traj_info = trajectory_info

        # Use same trajectory twice
        p1ss = [[0], [0]]
        p2ss = [[5], [5]]

        result = distance(
            [traj_info, traj_info],
            [top_info, top_info],
            p1ss, p2ss,
            ncpus=1
        )

        assert len(result) == 2, "Should have two trajectory results"
        np.testing.assert_allclose(result[0][0], result[1][0],
                                   err_msg="Same trajectory should give same distances")


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_input_parsing(self):
        """Test input option parsing with various configurations."""
        parser = cli_parser()

        # No input gives None
        args_empty = parser.parse_args([])
        assert args_empty.input is None, "No input should be None"

        # Basic input
        args_basic = parser.parse_args(["-i", "traj.dat", "0", "5"])
        assert args_basic.input == [["traj.dat", "0", "5"]], "Basic input parsing"

        # Multiple pairs
        args_multi = parser.parse_args(["-i", "traj.dat", "0", "5", "1", "6"])
        assert args_multi.input == [["traj.dat", "0", "5", "1", "6"]], "Multiple pairs"

        # Multiple trajectories
        args_trajs = parser.parse_args([
            "-i", "traj1.dat", "0", "5",
            "-i", "traj2.dat", "1", "6"
        ])
        assert len(args_trajs.input) == 2, "Multiple trajectories"

    def test_parser_all_options(self):
        """Test parser with all options and defaults."""
        parser = cli_parser()

        # All options
        args = parser.parse_args([
            "-o", "output.png",
            "-f", "both",
            "-d", "data.json",
            "-n", "dist1", "dist2",
            "-p", "4",
            "-c",
            "-q",
            "-i", "traj.dat", "0", "5", "1", "6"
        ])

        assert args.output == ["output.png"], "Output option"
        assert args.format == ["both"], "Format option"
        assert args.data == ["data.json"], "Data option"
        assert args.names == [["dist1", "dist2"]], "Names option"
        assert args.parallel == [4], "Parallel option"
        assert args.cluster is True, "Cluster option"
        assert args.quiet is True, "Quiet option"

        # Defaults
        args_defaults = parser.parse_args(["-i", "traj.dat", "0", "5"])
        assert args_defaults.output is None
        assert args_defaults.format is None
        assert args_defaults.data is None
        assert args_defaults.names is None
        assert args_defaults.parallel is None
        assert args_defaults.cluster is False
        assert args_defaults.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_histogram_output(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates histogram output file."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "distance.png"

        test_args = [
            "distance.py",
            "-i", str(traj_copy), "0", "5",
            "-o", str(output_file),
            "-f", "histogram"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Histogram output should be created"

    def test_main_trajectory_output(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates trajectory plot output file."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "distance_traj.png"

        test_args = [
            "distance.py",
            "-i", str(traj_copy), "0", "5",
            "-o", str(output_file),
            "-f", "trajectory"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Trajectory plot should be created"

    def test_main_both_format(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates both histogram and trajectory plots."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "distance.png"

        test_args = [
            "distance.py",
            "-i", str(traj_copy), "0", "5",
            "-o", str(output_file),
            "-f", "both"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        hist_file = temp_output_dir / "distance_hist.png"
        traj_file = temp_output_dir / "distance_traj.png"
        assert hist_file.exists(), "Histogram should be created"
        assert traj_file.exists(), "Trajectory plot should be created"

    def test_main_json_data_output(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates valid JSON data file."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        data_file = temp_output_dir / "data.json"

        test_args = [
            "distance.py",
            "-i", str(traj_copy), "0", "5",
            "-o", str(temp_output_dir / "distance.png"),
            "-d", str(data_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert data_file.exists(), "JSON file should be created"

        # Verify JSON structure
        with open(data_file) as f:
            data = json.load(f)
        assert "0-5" in data, "Key should be particle pair"

    def test_main_default_output_and_multiple_pairs(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test default output filename and multiple particle pairs."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Default filename
        test_args_default = ["distance.py", "-i", str(traj_copy), "0", "5"]
        with patch.object(sys, 'argv', test_args_default):
            main()
        assert (temp_output_dir / "distance.png").exists(), "Default output should be created"

        # Multiple pairs
        output_multi = temp_output_dir / "multi_distance.png"
        test_args_multi = [
            "distance.py",
            "-i", str(traj_copy), "0", "5", "1", "6", "2", "7",
            "-o", str(output_multi)
        ]
        with patch.object(sys, 'argv', test_args_multi):
            main()
        assert output_multi.exists(), "Multiple pairs should work"

    def test_main_with_options(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() with parallel, quiet, and custom names options."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Quiet mode
        output_quiet = temp_output_dir / "quiet.png"
        test_args_quiet = [
            "distance.py", "-q",
            "-i", str(traj_copy), "0", "5",
            "-o", str(output_quiet)
        ]
        with patch.object(sys, 'argv', test_args_quiet):
            main()
        assert output_quiet.exists(), "Quiet mode should work"

        # Parallel
        output_parallel = temp_output_dir / "parallel.png"
        test_args_parallel = [
            "distance.py", "-p", "2",
            "-i", str(traj_copy), "0", "5",
            "-o", str(output_parallel)
        ]
        with patch.object(sys, 'argv', test_args_parallel):
            main()
        assert output_parallel.exists(), "Parallel should work"

        # Custom names
        data_file = temp_output_dir / "named_data.json"
        test_args_names = [
            "distance.py",
            "-i", str(traj_copy), "0", "5", "1", "6",
            "-n", "first_dist", "second_dist",
            "-o", str(temp_output_dir / "named.png"),
            "-d", str(data_file)
        ]
        with patch.object(sys, 'argv', test_args_names):
            main()
        assert data_file.exists(), "Custom names should work"

    def test_main_error_cases(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() error handling for invalid format and missing input."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Invalid format
        test_args_invalid = [
            "distance.py",
            "-i", str(traj_copy), "0", "5",
            "-f", "invalid_format"
        ]
        with patch.object(sys, 'argv', test_args_invalid):
            with pytest.raises(RuntimeError, match="Unrecognized graph format"):
                main()

        # Missing input
        test_args_no_input = ["distance.py"]
        with patch.object(sys, 'argv', test_args_no_input):
            with pytest.raises(SystemExit):
                main()


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_min_image_extreme_boxes(self):
        """Test min_image with very small and very large box sizes."""
        p1 = np.array([0.0, 0.0, 0.0])

        # Very small box
        p2_small = np.array([0.3, 0.0, 0.0])
        assert min_image(p1, p2_small, 1.0) == pytest.approx(0.3), "Small box"

        # Very large box
        p2_large = np.array([100.0, 0.0, 0.0])
        assert min_image(p1, p2_large, 1000000.0) == pytest.approx(100.0), "Large box"

    def test_distance_particle_spacing(self, trajectory_info):
        """Test distances between adjacent and distant particles."""
        top_info, traj_info = trajectory_info

        # Adjacent nucleotides should be close (~0.5-0.7 nm backbone)
        result_adjacent = distance([traj_info], [top_info], [[0]], [[1]], ncpus=1)
        avg_adj = np.mean(result_adjacent[0][0])
        assert 0.1 < avg_adj < 2.0, "Adjacent nucleotides should be close"

        # Distant particles should have non-zero distances
        max_particle = top_info.nbases - 1
        result_distant = distance([traj_info], [top_info], [[0]], [[max_particle]], ncpus=1)
        assert all(d > 0 for d in result_distant[0][0]), "Distant particles should have positive distance"

    def test_distance_variance_positive(self, trajectory_info):
        """Test that distance variance is non-negative across trajectory."""
        top_info, traj_info = trajectory_info

        result = distance([traj_info], [top_info], [[0]], [[10]], ncpus=1)

        variance = np.var(result[0][0])
        assert variance >= 0, "Variance should be non-negative"
