"""
Tests for oxDNA_analysis_tools.subset_trajectory module.

Tests cover:
- compute() helper function
- write_topologies() function
- subset() API function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.subset_trajectory import compute, write_topologies, subset, cli_parser, ComputeContext, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, strand_describe


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
# Unit Tests - compute()
# =============================================================================

class TestCompute:
    """Tests for the compute() helper function."""

    def test_compute_behavior(self, trajectory_info):
        """Test compute() returns list of non-empty strings for each index set."""
        top_info, traj_info = trajectory_info

        # Two sets of indexes
        indexes = [
            list(range(10)),  # First 10 particles
            list(range(10, 20))  # Next 10 particles
        ]

        ctx = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            indexes=indexes
        )

        result = compute(ctx, chunk_size=1, chunk_id=0)

        assert isinstance(result, list), "Should return list"
        assert len(result) == len(indexes), "Should have one string per index set"
        for i, r in zip(indexes, result):
            assert isinstance(r, str), "Each element should be string"
            n_lines = len(r.split('\n'))
            assert n_lines == len(i)+4, "Number of particles in output should match index count (3-line header plus newline at end)"

        # Test single index set
        ctx_single = ComputeContext(
            traj_info=traj_info,
            top_info=top_info,
            indexes=[list(range(20))]
        )
        result_single = compute(ctx_single, chunk_size=1, chunk_id=0)
        assert len(result_single) == 1, "Should have one result"


# =============================================================================
# Unit Tests - write_topologies()
# =============================================================================

class TestWriteTopologies:
    """Tests for the write_topologies() function."""

    def test_write_topologies_behavior(self, system_info, temp_output_dir):
        """Test write_topologies creates files with correct content."""
        indexes = [list(range(20)), list(range(20, 40))]
        outfiles = [str(temp_output_dir / "sub1"), str(temp_output_dir / "sub2")]

        # Use old_format since test topology is old-style
        old_format = system_info.strands[0].is_old()
        top_names = write_topologies(system_info, indexes, outfiles, old_format=old_format)

        assert len(top_names) == 2, "Should return two topology names"
        for name in top_names:
            assert Path(name).exists(), f"Topology file {name} should exist"

        # Test with single subset to verify content
        indexes_single = [list(range(10))]
        outfiles_single = [str(temp_output_dir / "subset")]
        top_names_single = write_topologies(system_info, indexes_single, outfiles_single, old_format=old_format)

        with open(top_names_single[0]) as f:
            lines = f.readlines()

        header = lines[0].strip().split()
        n_particles = int(header[0])
        assert n_particles == 10, "Topology should have 10 particles"


# =============================================================================
# API Tests - subset() function
# =============================================================================

class TestSubsetFunction:
    """Tests for the subset() API function."""

    def test_subset_behavior(self, trajectory_info, system_info, temp_output_dir):
        """Test subset() creates readable trajectory and topology files with correct particle counts."""
        top_info, traj_info = trajectory_info

        # Single subset
        n_particles = 25
        indexes = [list(range(n_particles))]
        outfiles = [str(temp_output_dir / "subset1")]

        subset(traj_info, top_info, system_info, indexes, outfiles, ncpus=1)

        # Check files exist
        dat_file = Path(str(temp_output_dir / "subset1.dat"))
        top_file = Path(str(temp_output_dir / "subset1.top"))
        assert dat_file.exists(), "Trajectory file should be created"
        assert top_file.exists(), "Topology file should be created"

        # Read and verify
        sub_top_info, sub_traj_info = describe(
            str(temp_output_dir / "subset1.top"),
            str(temp_output_dir / "subset1.dat")
        )
        assert sub_top_info.nbases == n_particles, "Subset should have correct particle count"
        assert sub_traj_info.nconfs == traj_info.nconfs, "Should have same number of configurations"

        # Multiple subsets
        indexes_multi = [
            list(range(0, 30)),
            list(range(30, 60)),
            list(range(60, 90))
        ]
        outfiles_multi = [
            str(temp_output_dir / "part1"),
            str(temp_output_dir / "part2"),
            str(temp_output_dir / "part3")
        ]

        subset(traj_info, top_info, system_info, indexes_multi, outfiles_multi, ncpus=1)

        for name in outfiles_multi:
            assert Path(f"{name}.dat").exists(), f"{name}.dat should exist"
            assert Path(f"{name}.top").exists(), f"{name}.top should exist"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_behavior(self):
        """Test parser requires arguments, handles -i option, and has correct defaults."""
        parser = cli_parser()

        # Requires trajectory and topology
        with pytest.raises(SystemExit):
            parser.parse_args([])
        with pytest.raises(SystemExit):
            parser.parse_args(["traj.dat"])

        # Index option
        args = parser.parse_args([
            "traj.dat", "top.top",
            "-i", "index1.txt", "out1",
            "-i", "index2.txt", "out2"
        ])

        assert args.trajectory == "traj.dat", "Trajectory not parsed"
        assert args.topology == "top.top", "Topology not parsed"
        assert len(args.index) == 2, "Should have two index/output pairs"
        assert args.index[0] == ["index1.txt", "out1"], "First pair not parsed"

        # All options
        args_all = parser.parse_args([
            "-p", "4",
            "-f",
            "-q",
            "traj.dat", "top.top",
            "-i", "index.txt", "output"
        ])

        assert args_all.parallel == 4, "Parallel option not parsed"
        assert args_all.old_format is True, "Old format option not parsed"
        assert args_all.quiet is True, "Quiet option not parsed"

        # Defaults
        args_defaults = parser.parse_args(["traj.dat", "top.top", "-i", "index.txt", "out"])
        assert args_defaults.parallel is None
        assert args_defaults.old_format is False
        assert args_defaults.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_basic_and_multiple_subsets(self, mini_traj_path, topology_path, temp_output_dir, monkeypatch):
        """Test main() creates output files for single and multiple subsets."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "top.top"
        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)

        # Single subset
        index_file = temp_output_dir / "index.txt"
        index_file.write_text(" ".join(str(i) for i in range(20)))

        test_args = [
            "subset_trajectory.py",
            str(traj_copy), str(top_copy),
            "-i", str(index_file), "subset"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert (temp_output_dir / "subset.dat").exists(), "Trajectory should be created"
        assert (temp_output_dir / "subset.top").exists(), "Topology should be created"

        # Multiple subsets
        index1 = temp_output_dir / "index1.txt"
        index2 = temp_output_dir / "index2.txt"
        index1.write_text(" ".join(str(i) for i in range(30)))
        index2.write_text(" ".join(str(i) for i in range(30, 60)))

        test_args_multi = [
            "subset_trajectory.py",
            str(traj_copy), str(top_copy),
            "-i", str(index1), "part1",
            "-i", str(index2), "part2"
        ]

        with patch.object(sys, 'argv', test_args_multi):
            main()

        assert (temp_output_dir / "part1.dat").exists(), "First subset should be created"
        assert (temp_output_dir / "part2.dat").exists(), "Second subset should be created"

    def test_main_invalid_index_file_raises(self, mini_traj_path, topology_path, temp_output_dir, monkeypatch):
        """Test that main() raises error for invalid index file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "top.top"
        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)

        # Create invalid index file
        index_file = temp_output_dir / "bad_index.txt"
        index_file.write_text("not integers")

        test_args = [
            "subset_trajectory.py",
            str(traj_copy), str(top_copy),
            "-i", str(index_file), "bad_output"
        ]

        with patch.object(sys, 'argv', test_args):
            with pytest.raises(RuntimeError, match="index file.*must be a space-seperated list"):
                main()


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_subset_edge_cases(self, trajectory_info, system_info, temp_output_dir):
        """Test subset() with small, non-contiguous indexes, and high ncpus."""
        top_info, traj_info = trajectory_info

        # Very small subset
        indexes_small = [list(range(5))]
        outfiles_small = [str(temp_output_dir / "small")]
        subset(traj_info, top_info, system_info, indexes_small, outfiles_small, ncpus=1)

        sub_top_small, _ = describe(
            str(temp_output_dir / "small.top"),
            str(temp_output_dir / "small.dat")
        )
        assert sub_top_small.nbases == 5, "Should have 5 particles"

        # Non-contiguous indexes
        indexes_sparse = [[0, 5, 10, 15, 20, 25, 30]]
        outfiles_sparse = [str(temp_output_dir / "sparse")]
        subset(traj_info, top_info, system_info, indexes_sparse, outfiles_sparse, ncpus=1)

        sub_top_sparse, _ = describe(
            str(temp_output_dir / "sparse.top"),
            str(temp_output_dir / "sparse.dat")
        )
        assert sub_top_sparse.nbases == 7, "Should have 7 particles"

        # High ncpus
        indexes_cpu = [list(range(20))]
        outfiles_cpu = [str(temp_output_dir / "highcpu")]
        subset(traj_info, top_info, system_info, indexes_cpu, outfiles_cpu, ncpus=100)
        assert (temp_output_dir / "highcpu.dat").exists(), "Should handle high ncpus"
