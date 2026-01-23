"""
Tests for oxDNA_analysis_tools.pca module.

Tests cover:
- align_positions() function
- pca() API function
- CLI argument parsing
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.pca import align_positions, pca, cli_parser, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, inbox


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
    """Get mean configuration for PCA."""
    top_info, _ = trajectory_info
    _, mean_info = describe(None, str(mean_conf_path))
    conf = get_confs(top_info, mean_info, 0, 1)[0]
    # Center the configuration
    cms = np.mean(conf.positions, axis=0)
    conf.positions -= cms
    return conf


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# align_positions() Tests
# =============================================================================

class TestAlignPositions:
    """Tests for the align_positions() SVD alignment function."""

    def test_align_positions_behavior(self):
        """Test align_positions returns correct type, preserves shape, and produces valid values."""
        # Test basic return type
        coords = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0]
        ])
        ref = coords - np.mean(coords, axis=0)

        result = align_positions(ref, coords.copy())

        assert isinstance(result, np.ndarray), "Should return numpy array"
        assert result.dtype in [np.float64, np.float32], "Should return float array"

        # Test shape preservation with larger array
        n_particles = 50
        ref_large = np.random.rand(n_particles, 3)
        ref_large = ref_large - np.mean(ref_large, axis=0)
        coords_large = np.random.rand(n_particles, 3)

        result_large = align_positions(ref_large, coords_large)
        assert result_large.shape == (n_particles, 3), "Shape should be preserved"

        # Test no NaN/Inf with offset coordinates
        coords_offset = np.random.rand(n_particles, 3) + np.array([10, 20, 30])
        result_offset = align_positions(ref_large, coords_offset)

        assert not np.any(np.isnan(result_offset)), "Result should not contain NaN"
        assert not np.any(np.isinf(result_offset)), "Result should not contain Inf"


# =============================================================================
# pca() API Tests
# =============================================================================

class TestPcaFunction:
    """Tests for the pca() API function."""

    def test_pca_behavior(self, trajectory_info, mean_conf):
        """Test pca() returns correct shapes, sorted descending eigenvalues, and non-negative values."""
        import matplotlib
        matplotlib.use('Agg')

        top_info, traj_info = trajectory_info

        coordinates, evalues, evectors = pca(traj_info, top_info, mean_conf, ncpus=1)

        n_dims = top_info.nbases * 3

        # Check shapes
        assert coordinates.shape == (traj_info.nconfs, n_dims), \
            f"Coordinates should be ({traj_info.nconfs}, {n_dims})"
        assert evalues.shape == (n_dims,), f"Eigenvalues should be ({n_dims},)"
        assert evectors.shape == (n_dims, n_dims), f"Eigenvectors should be ({n_dims}, {n_dims})"

        # Eigenvalues should be sorted descending (largest first)
        real_evalues = np.real(evalues)
        assert np.all(real_evalues[:-1] >= real_evalues[1:] - 1e-10), \
            "Eigenvalues should be sorted descending"

        # Real parts should be non-negative (small negative due to numerical errors ok)
        assert np.all(real_evalues >= -1e-10), "Eigenvalues should be non-negative"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_arguments(self):
        """Test that parser requires trajectory, meanfile, and outfile arguments."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-p", "4",
            "-c",
            "-n", "3",
            "-q",
            "trajectory.dat",
            "mean.dat",
            "output.json"
        ])

        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.meanfile == ["mean.dat"], "Meanfile not parsed"
        assert args.outfile == ["output.json"], "Outfile not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.cluster is True, "Cluster option not parsed"
        assert args.N == [3], "N option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["traj.dat", "mean.dat", "out.json"])
        assert args_defaults.quiet is False, "Quiet should default to False"
        assert args_defaults.cluster is False, "Cluster should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output_files(self, mini_traj_path, mean_conf_path, temp_output_dir, monkeypatch):
        """Test main() creates expected output files."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        mean_copy = temp_output_dir / "mean.dat"
        copy(mini_traj_path, traj_copy)
        copy(mean_conf_path, mean_copy)

        output_file = temp_output_dir / "pca_output.json"

        test_args = [
            "pca.py",
            str(traj_copy),
            str(mean_copy),
            str(output_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check that component file was created (default N=1)
        component_file = temp_output_dir / "pca_output0.json"
        assert component_file.exists(), "Component 0 JSON should be created"

        # Check scree plot was created
        assert (temp_output_dir / "scree.png").exists(), "Scree plot should be created"
