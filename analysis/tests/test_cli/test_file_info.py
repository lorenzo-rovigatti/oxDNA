"""
Tests for oxDNA_analysis_tools.file_info module.

Tests cover:
- file_info() API function
- print_info() function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch
from io import StringIO

import pytest

from oxDNA_analysis_tools.file_info import file_info, print_info, cli_parser, main
from oxDNA_analysis_tools.UTILS.RyeReader import describe


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
# API Tests - file_info() function
# =============================================================================

class TestFileInfoFunction:
    """Tests for the file_info() API function."""

    def test_file_info_behavior(self, mini_traj_path, trajectory_info):
        """Test file_info() returns dict with expected keys, correct values, multiple files, and valid times."""
        top_info, traj_info = trajectory_info

        # Single file
        result = file_info([str(mini_traj_path)])

        # Check all expected keys exist
        expected_keys = ['name', 'particles', 'n_confs', 'traj_size', 't_start', 't_end']
        for key in expected_keys:
            assert key in result, f"Result should have '{key}' key"

        # Check values are correct
        assert result['particles'][0] == top_info.nbases, "Particle count should match"
        assert result['n_confs'][0] == traj_info.nconfs, "Conf count should match"
        assert result['traj_size'][0] > 0, "File size should be positive"
        assert result['t_start'][0] <= result['t_end'][0], "Start time should be <= end time"

        # Multiple files
        result_multi = file_info([str(mini_traj_path), str(mini_traj_path)])
        assert len(result_multi['particles']) == 2, "Should have two particle counts"
        assert len(result_multi['n_confs']) == 2, "Should have two conf counts"
        assert result_multi['particles'][0] == result_multi['particles'][1], "Same file should have same particles"


# =============================================================================
# print_info() Tests
# =============================================================================

class TestPrintInfo:
    """Tests for the print_info() function."""

    def test_print_info_outputs_table(self, mini_traj_path, capsys):
        """Test that print_info() produces formatted output."""
        info = file_info([str(mini_traj_path)])
        labels = ["test_traj"]

        print_info(info, labels)

        captured = capsys.readouterr()
        # Should contain column headers and data
        assert "name" in captured.out, "Output should have 'name' header"
        assert "particles" in captured.out, "Output should have 'particles' header"
        assert "test_traj" in captured.out, "Output should have label"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_trajectories(self):
        """Test that parser requires at least one trajectory argument."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-l", "traj1", "traj2",
            "-q",
            "file1.dat",
            "file2.dat"
        ])

        assert args.trajectories == ["file1.dat", "file2.dat"], "Trajectories not parsed"
        assert args.labels == ["traj1", "traj2"], "Labels not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["file.dat"])
        assert args_defaults.labels is None, "Labels should default to None"
        assert args_defaults.quiet is False, "Quiet should default to False"

    def test_parser_multiple_trajectories(self):
        """Test parser accepts multiple trajectory files."""
        parser = cli_parser()
        args = parser.parse_args(["file1.dat", "file2.dat", "file3.dat"])

        assert len(args.trajectories) == 3, "Should accept multiple trajectories"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_single_trajectory(self, mini_traj_path, temp_output_dir, monkeypatch, capsys):
        """Test main() with single trajectory produces output."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        test_args = ["file_info.py", str(traj_copy)]

        with patch.object(sys, 'argv', test_args):
            main()

        captured = capsys.readouterr()
        # Should produce table output
        assert "particles" in captured.out, "Output should have particles column"

    def test_main_with_labels(self, mini_traj_path, temp_output_dir, monkeypatch, capsys):
        """Test main() with custom labels."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Note: trajectory must come after labels since -l consumes multiple args
        test_args = ["file_info.py", str(traj_copy), "-l", "custom_label"]

        with patch.object(sys, 'argv', test_args):
            main()

        captured = capsys.readouterr()
        assert "custom_label" in captured.out, "Output should have custom label"

    def test_main_label_count_mismatch_raises(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() raises error when label count doesn't match trajectory count."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        # Provide two labels but only one trajectory
        test_args = ["file_info.py", str(traj_copy), "-l", "label1", "label2"]

        with patch.object(sys, 'argv', test_args):
            with pytest.raises(RuntimeError, match="Number of trajectories does not match"):
                main()
