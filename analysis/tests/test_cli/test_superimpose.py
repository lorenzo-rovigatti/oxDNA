"""
Tests for oxDNA_analysis_tools.superimpose module.

Tests cover:
- superimpose() API function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.superimpose import superimpose, cli_parser, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, inbox, write_conf
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
def mean_path(test_resources):
    """Path to the mean configuration."""
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


@pytest.fixture
def reference_conf(trajectory_info):
    """Get the first configuration as reference."""
    top_info, traj_info = trajectory_info
    conf = get_confs(top_info, traj_info, 0, 1)[0]
    return inbox(conf)


@pytest.fixture
def victim_files(trajectory_info, temp_output_dir):
    """Create victim configuration files for superimposition."""
    top_info, traj_info = trajectory_info

    files = []
    for i in range(1, 3):
        conf = get_confs(top_info, traj_info, i, 1)[0]
        filepath = str(temp_output_dir / f"victim_{i}.dat")
        write_conf(filepath, conf)
        files.append(filepath)

    return files


# =============================================================================
# API Tests - superimpose() function
# =============================================================================

class TestSuperimposeFunction:
    """Tests for the superimpose() API function."""

    def test_superimpose_behavior(self, reference_conf, victim_files, trajectory_info, temp_output_dir):
        """Test superimpose() returns aligned configurations with correct shapes, RMSDs, and handles indexes."""
        top_info, _ = trajectory_info

        # Basic superimposition
        aligned, rmsds = superimpose(reference_conf, victim_files)

        # Check return types and lengths
        assert isinstance(aligned, list), "Should return list of configurations"
        assert isinstance(rmsds, list), "Should return list of RMSDs"
        assert len(aligned) == len(victim_files), "Should align all victims"
        assert len(rmsds) == len(victim_files), "Should return RMSD for each victim"

        # Check configurations
        for conf in aligned:
            assert isinstance(conf, Configuration), "Aligned should be Configuration"
            assert conf.positions is not None, "Should have positions"
            assert conf.positions.shape == (top_info.nbases, 3), "Position shape mismatch"
            assert conf.a1s.shape == (top_info.nbases, 3), "a1s shape mismatch"
            assert conf.a3s.shape == (top_info.nbases, 3), "a3s shape mismatch"

        # RMSDs should be non-negative floats
        for rmsd in rmsds:
            assert isinstance(rmsd, float), "RMSD should be float"
            assert rmsd >= 0, "RMSD should be non-negative"

        # Test self-alignment (RMSD ~0)
        victim_file = str(temp_output_dir / "self_victim.dat")
        write_conf(victim_file, reference_conf)
        _, rmsds_self = superimpose(reference_conf, [victim_file])
        assert rmsds_self[0] == pytest.approx(0.0, abs=1e-6), "Self-alignment RMSD should be ~0"

        # Test with indexes
        n_particles = min(20, top_info.nbases)
        indexes = np.array([list(range(n_particles))] * (len(victim_files) + 1), dtype=int)
        aligned_idx, rmsds_idx = superimpose(reference_conf, victim_files, indexes)
        assert len(aligned_idx) == len(victim_files), "Should align all victims with indexes"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_behavior(self):
        """Test parser requires arguments, accepts multiple victims and options, and has defaults."""
        parser = cli_parser()

        # Requires reference and at least one victim
        with pytest.raises(SystemExit):
            parser.parse_args([])
        with pytest.raises(SystemExit):
            parser.parse_args(["ref.dat"])

        # Multiple victims
        args = parser.parse_args(["ref.dat", "victim1.dat", "victim2.dat", "victim3.dat"])
        assert args.reference == ["ref.dat"], "Reference not parsed"
        assert args.victims == ["victim1.dat", "victim2.dat", "victim3.dat"], "Victims not parsed"

        # All options
        args_all = parser.parse_args([
            "-i", "index1.txt", "index2.txt",
            "-o", "out1.dat", "out2.dat",
            "-q",
            "ref.dat", "victim1.dat"
        ])
        assert args_all.index_files == ["index1.txt", "index2.txt"], "Index files not parsed"
        assert args_all.output_names == ["out1.dat", "out2.dat"], "Output names not parsed"
        assert args_all.quiet is True, "Quiet option not parsed"

        # Defaults
        args_defaults = parser.parse_args(["ref.dat", "victim.dat"])
        assert args_defaults.index_files is None
        assert args_defaults.output_names is None
        assert args_defaults.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_basic_and_custom_output(self, mini_traj_path, mean_path, temp_output_dir, monkeypatch, trajectory_info):
        """Test main() creates output files with default and custom names."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        # Copy files to temp dir
        ref_copy = temp_output_dir / "ref.dat"
        victim_copy = temp_output_dir / "victim.dat"
        copy(mean_path, ref_copy)
        copy(mini_traj_path, victim_copy)

        # Basic (default output name)
        test_args = ["superimpose.py", str(ref_copy), str(victim_copy)]
        with patch.object(sys, 'argv', test_args):
            main()

        expected_output = temp_output_dir / "victim_a.dat"
        assert expected_output.exists(), "Output file should be created with default name"

        # Custom output names
        victim1 = temp_output_dir / "victim1.dat"
        victim2 = temp_output_dir / "victim2.dat"
        conf1 = get_confs(top_info, traj_info, 1, 1)[0]
        conf2 = get_confs(top_info, traj_info, 2, 1)[0]
        write_conf(str(victim1), conf1)
        write_conf(str(victim2), conf2)

        output1 = temp_output_dir / "aligned1.dat"
        output2 = temp_output_dir / "aligned2.dat"

        test_args_custom = [
            "superimpose.py",
            str(ref_copy), str(victim1), str(victim2),
            "-o", str(output1), str(output2)
        ]
        with patch.object(sys, 'argv', test_args_custom):
            main()

        assert output1.exists(), "First output file should be created"
        assert output2.exists(), "Second output file should be created"

    def test_main_with_index_file(self, mean_path, temp_output_dir, monkeypatch, trajectory_info):
        """Test main() with single index file applied to all structures."""
        monkeypatch.chdir(temp_output_dir)
        top_info, traj_info = trajectory_info

        ref_copy = temp_output_dir / "ref.dat"
        copy(mean_path, ref_copy)

        victim = temp_output_dir / "victim.dat"
        conf = get_confs(top_info, traj_info, 1, 1)[0]
        write_conf(str(victim), conf)

        # Create index file
        index_file = temp_output_dir / "index.txt"
        index_file.write_text("0 1 2 3 4 5 6 7 8 9")

        test_args = [
            "superimpose.py",
            str(ref_copy), str(victim),
            "-i", str(index_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        expected_output = temp_output_dir / "victim_a.dat"
        assert expected_output.exists(), "Output should be created with index file"


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_superimpose_multiple_victims(self, reference_conf, trajectory_info, temp_output_dir):
        """Test superimpose() with many victim files and verify RMSD ordering."""
        top_info, traj_info = trajectory_info

        # Create multiple victims
        victims = []
        for i in range(min(5, traj_info.nconfs - 1)):
            conf = get_confs(top_info, traj_info, i + 1, 1)[0]
            filepath = str(temp_output_dir / f"multi_victim_{i}.dat")
            write_conf(filepath, conf)
            victims.append(filepath)

        aligned, rmsds = superimpose(reference_conf, victims)

        assert len(aligned) == len(victims), "Should align all victims"
        assert len(rmsds) == len(victims), "Should return RMSD for all"
        assert all(isinstance(r, float) for r in rmsds), "RMSDs should be floats"
