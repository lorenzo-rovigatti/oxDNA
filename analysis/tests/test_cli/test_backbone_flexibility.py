"""
Tests for oxDNA_analysis_tools.backbone_flexibility module.

Tests cover:
- rad2degree() utility function
- compute() helper function
- backbone_flexibility() API function
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

from oxDNA_analysis_tools.backbone_flexibility import (
    rad2degree,
    compute,
    backbone_flexibility,
    cli_parser,
    ComputeContext,
    main
)
from oxDNA_analysis_tools.UTILS.RyeReader import describe, strand_describe


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
def topology_path(test_resources):
    """Path to the topology file."""
    return test_resources / "rna_tile.top"


@pytest.fixture(scope="module")
def trajectory_info(topology_path, mini_traj_path):
    """Get topology and trajectory info."""
    top_info, traj_info = describe(str(topology_path), str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture(scope="module")
def system_info(topology_path):
    """Get system info from topology."""
    system, _ = strand_describe(str(topology_path))
    return system


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# Unit Tests - rad2degree()
# =============================================================================

class TestRad2Degree:
    """Tests for the rad2degree() utility function."""

    def test_rad2degree_conversions(self):
        """Test rad2degree converts various angles correctly."""
        # 0 radians -> 0 degrees
        assert rad2degree(0) == pytest.approx(0.0), "0 rad should be 0 deg"

        # pi radians -> 180 degrees
        assert rad2degree(np.pi) == pytest.approx(180.0), "pi rad should be 180 deg"

        # pi/2 radians -> 90 degrees
        assert rad2degree(np.pi / 2) == pytest.approx(90.0), "pi/2 rad should be 90 deg"

        # 2*pi radians -> 360 degrees
        assert rad2degree(2 * np.pi) == pytest.approx(360.0), "2*pi rad should be 360 deg"

        # negative angles
        assert rad2degree(-np.pi) == pytest.approx(-180.0), "-pi rad should be -180 deg"


# =============================================================================
# Unit Tests - compute()
# =============================================================================

class TestCompute:
    """Tests for the compute() helper function."""

    def test_compute_behavior(self, trajectory_info, system_info):
        """Test compute() returns correct tuple with proper shapes and valid angle ranges."""
        top_info, traj_info = trajectory_info

        ctx = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            system=system_info
        )

        result = compute(ctx, chunk_size=1, chunk_id=0)

        # Check return type
        assert isinstance(result, tuple), "Should return tuple"
        assert len(result) == 2, "Should return (torsions, dihedrals)"

        torsions, dihedrals = result

        # Check shapes
        n_strands = len(system_info.strands)
        expected_torsions = top_info.nbases - (2 * n_strands)
        expected_dihedrals = top_info.nbases - (3 * n_strands)

        assert len(torsions) == expected_torsions, f"Expected {expected_torsions} torsions"
        assert len(dihedrals) == expected_dihedrals, f"Expected {expected_dihedrals} dihedrals"

        # Check angle ranges (from arccos, should be 0-180)
        assert np.all(torsions >= 0), "Torsion angles should be non-negative"
        assert np.all(torsions <= 180), "Torsion angles should be <= 180"
        assert np.all(dihedrals >= 0), "Dihedral angles should be non-negative"
        assert np.all(dihedrals <= 180), "Dihedral angles should be <= 180"


# =============================================================================
# API Tests - backbone_flexibility() function
# =============================================================================

class TestBackboneFlexibilityFunction:
    """Tests for the backbone_flexibility() API function."""

    def test_backbone_flexibility_behavior(self, trajectory_info, system_info):
        """Test backbone_flexibility() returns correct arrays with valid values."""
        top_info, traj_info = trajectory_info

        torsions, dihedrals = backbone_flexibility(traj_info, top_info, system_info, ncpus=1)

        # Check return types
        assert isinstance(torsions, np.ndarray), "Torsions should be numpy array"
        assert isinstance(dihedrals, np.ndarray), "Dihedrals should be numpy array"

        # Check shapes
        n_strands = len(system_info.strands)
        expected_torsions = top_info.nbases - (2 * n_strands)
        expected_dihedrals = top_info.nbases - (3 * n_strands)

        assert len(torsions) == expected_torsions, "Torsion count mismatch"
        assert len(dihedrals) == expected_dihedrals, "Dihedral count mismatch"

        # Check mean angles are in valid range (averaged over trajectory)
        mean_torsion = np.mean(torsions)
        mean_dihedral = np.mean(dihedrals)

        assert 0 < mean_torsion < 180, "Mean torsion should be in valid range"
        assert 0 < mean_dihedral < 180, "Mean dihedral should be in valid range"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_behavior(self):
        """Test parser requires arguments and accepts all options with correct defaults."""
        parser = cli_parser()

        # Requires topology and trajectory
        with pytest.raises(SystemExit):
            parser.parse_args([])

        with pytest.raises(SystemExit):
            parser.parse_args(["top.top"])

        # All options
        args = parser.parse_args([
            "-p", "4",
            "-o", "plot.png",
            "-d", "data.json",
            "-q",
            "topology.top",
            "trajectory.dat"
        ])

        assert args.topology == ["topology.top"], "Topology not parsed"
        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.output == ["plot.png"], "Output option not parsed"
        assert args.data == ["data.json"], "Data option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Defaults
        args_defaults = parser.parse_args(["top.top", "traj.dat"])
        assert args_defaults.parallel is None
        assert args_defaults.output is None
        assert args_defaults.data is None
        assert args_defaults.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_outputs(self, mini_traj_path, topology_path, temp_output_dir, monkeypatch):
        """Test main() creates plot and data files with various options."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "top.top"
        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)

        # Test with explicit output and data file
        output_file = temp_output_dir / "flex.png"
        data_file = temp_output_dir / "data.json"

        test_args = [
            "backbone_flexibility.py",
            "-o", str(output_file),
            "-d", str(data_file),
            str(top_copy), str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output plot should be created"
        assert data_file.exists(), "Data file should be created"

        # Verify JSON structure
        with open(data_file) as f:
            data = json.load(f)

        assert "torsions" in data, "Data should have torsions"
        assert "dihedrals" in data, "Data should have dihedrals"
        assert isinstance(data["torsions"], list), "Torsions should be list"
        assert isinstance(data["dihedrals"], list), "Dihedrals should be list"

    def test_main_default_output(self, mini_traj_path, topology_path, temp_output_dir, monkeypatch):
        """Test main() uses default output filename."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "top.top"
        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)

        test_args = ["backbone_flexibility.py", str(top_copy), str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        default_output = temp_output_dir / "ramachandran.png"
        assert default_output.exists(), "Default output 'ramachandran.png' should be created"

    def test_main_with_parallel_and_quiet(self, mini_traj_path, topology_path, temp_output_dir, monkeypatch):
        """Test main() with parallel processing and quiet mode."""
        import matplotlib
        matplotlib.use('Agg')

        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "top.top"
        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)

        output_file = temp_output_dir / "parallel_flex.png"

        test_args = [
            "backbone_flexibility.py",
            "-p", "2",
            "-q",
            "-o", str(output_file),
            str(top_copy), str(traj_copy)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Parallel/quiet mode should create output"


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_backbone_flexibility_edge_cases(self, trajectory_info, system_info):
        """Test edge cases: single strand topology and high ncpus."""
        top_info, traj_info = trajectory_info

        # The test topology has 1 strand with 132 nucleotides
        assert len(system_info.strands) == 1, "Test assumes single strand"

        # Test single strand expected counts
        torsions, dihedrals = backbone_flexibility(traj_info, top_info, system_info, ncpus=1)
        assert len(torsions) == top_info.nbases - 2, "Single strand torsion count"
        assert len(dihedrals) == top_info.nbases - 3, "Single strand dihedral count"

        # Test high ncpus (more CPUs than configurations)
        torsions_high, dihedrals_high = backbone_flexibility(traj_info, top_info, system_info, ncpus=100)
        assert len(torsions_high) > 0, "Should produce torsions with high ncpus"
        assert len(dihedrals_high) > 0, "Should produce dihedrals with high ncpus"
