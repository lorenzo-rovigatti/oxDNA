"""
Tests for oxDNA_analysis_tools.decimate module.

Tests cover:
- decimate() API function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import pytest

from oxDNA_analysis_tools.decimate import decimate, cli_parser, main
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
# API Tests - decimate() function
# =============================================================================

class TestDecimateFunction:
    """Tests for the decimate() API function."""

    def test_decimate_behavior(self, mini_traj_path, temp_output_dir):
        """Test decimate with stride, start/stop, and particle preservation."""
        orig_top, orig_info = describe(None, str(mini_traj_path))

        # Test basic stride
        output_stride = temp_output_dir / "decimated.dat"
        decimate(
            traj=str(mini_traj_path),
            outfile=str(output_stride),
            start=0,
            stop=-1,
            stride=2
        )

        assert output_stride.exists(), "Output file should be created"
        _, dec_info = describe(None, str(output_stride))
        expected_confs = (orig_info.nconfs + 1) // 2
        assert dec_info.nconfs == expected_confs, \
            f"Expected {expected_confs} confs with stride=2, got {dec_info.nconfs}"

        # Test start/stop parameters
        output_range = temp_output_dir / "decimated_range.dat"
        start = 1
        stop = orig_info.nconfs - 1
        decimate(
            traj=str(mini_traj_path),
            outfile=str(output_range),
            start=start,
            stop=stop,
            stride=1
        )

        assert output_range.exists(), "Range output file should be created"
        _, range_info = describe(None, str(output_range))
        expected_range = stop - start
        assert range_info.nconfs == expected_range, \
            f"Expected {expected_range} confs with start/stop, got {range_info.nconfs}"

        # Test particle preservation
        output_particles = temp_output_dir / "decimated_particles.dat"
        decimate(
            traj=str(mini_traj_path),
            outfile=str(output_particles),
            start=0,
            stop=-1,
            stride=3
        )

        dec_top, _ = describe(None, str(output_particles))
        assert dec_top.nbases == orig_top.nbases, "Particle count should be preserved"


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
            "-s", "5",
            "-e", "100",
            "-d", "3",
            "-q",
            "input.dat",
            "output.dat"
        ])

        assert args.traj == "input.dat", "Trajectory not parsed"
        assert args.outfile == "output.dat", "Outfile not parsed"
        assert args.start == 5, "Start option not parsed"
        assert args.stop == 100, "Stop option not parsed"
        assert args.stride == 3, "Stride option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["input.dat", "output.dat"])
        assert args_defaults.start == 0, "Start should default to 0"
        assert args_defaults.stop == -1, "Stop should default to -1"
        assert args_defaults.stride == 10, "Stride should default to 10"
        assert args_defaults.quiet is False, "Quiet should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() creates output file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "decimated_main.dat"

        test_args = ["decimate.py", str(traj_copy), str(output_file)]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output file should be created"

    def test_main_with_options(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test main() with start, stop, and stride options."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        copy(mini_traj_path, traj_copy)

        output_file = temp_output_dir / "decimated_opts.dat"

        test_args = [
            "decimate.py",
            "-s", "1",
            "-d", "2",
            str(traj_copy),
            str(output_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output file should be created with options"
