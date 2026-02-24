"""
Tests for oxDNA_analysis_tools.generate_forces module.

Tests cover:
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from oxDNA_analysis_tools.generate_forces import cli_parser, main


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def test_resources():
    """Get the path to test resources directory."""
    return Path(__file__).parent.parent / "resources"


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_input_and_configuration(self):
        """Test that parser requires input file and configuration."""
        parser = cli_parser()

        with pytest.raises(SystemExit):
            parser.parse_args([])

        with pytest.raises(SystemExit):
            parser.parse_args(["input.txt"])

    def test_parser_all_options(self):
        """Test parser accepts all options."""
        parser = cli_parser()

        args = parser.parse_args([
            "-o", "forces.txt",
            "-f", "pairs.txt",
            "-s", "0.5",
            "-q",
            "input.txt",
            "conf.dat"
        ])

        assert args.inputfile == ["input.txt"], "Input file not parsed"
        assert args.configuration == ["conf.dat"], "Configuration not parsed"
        assert args.output == ["forces.txt"], "Output option not parsed"
        assert args.pairs == ["pairs.txt"], "Pairs option not parsed"
        assert args.stiff == [0.5], "Stiff option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

    def test_parser_defaults(self):
        """Test parser default values."""
        parser = cli_parser()

        args = parser.parse_args(["input.txt", "conf.dat"])

        assert args.output is None
        assert args.pairs is None
        assert args.stiff is None
        assert args.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_force_file(self, test_resources, temp_output_dir, monkeypatch):
        """Test main() creates forces output file."""
        monkeypatch.chdir(test_resources)

        output_file = temp_output_dir / "forces.txt"

        test_args = [
            "generate_forces.py",
            "-o", str(output_file),
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Forces file should be created"

    def test_main_default_output(self, test_resources, monkeypatch):
        """Test main() uses default output filename."""
        monkeypatch.chdir(test_resources)

        test_args = [
            "generate_forces.py",
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        default_output = test_resources / "forces.txt"
        assert default_output.exists(), "Default output 'forces.txt' should be created"
        # Clean up
        default_output.unlink()

    def test_main_creates_pairs_file(self, test_resources, temp_output_dir, monkeypatch):
        """Test main() creates pairs file when requested."""
        monkeypatch.chdir(test_resources)

        forces_file = temp_output_dir / "forces.txt"
        pairs_file = temp_output_dir / "pairs.txt"

        test_args = [
            "generate_forces.py",
            "-o", str(forces_file),
            "-f", str(pairs_file),
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert forces_file.exists(), "Forces file should be created"
        assert pairs_file.exists(), "Pairs file should be created"

    def test_main_custom_stiffness(self, test_resources, temp_output_dir, monkeypatch):
        """Test main() with custom stiffness value."""
        monkeypatch.chdir(test_resources)

        output_file = temp_output_dir / "forces.txt"

        test_args = [
            "generate_forces.py",
            "-o", str(output_file),
            "-s", "0.5",
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Forces file should be created"

        # Check that stiffness is set correctly
        with open(output_file) as f:
            content = f.read()
        assert "stiff = 0.5" in content, "Stiffness should be 0.5"

    def test_main_output_format(self, test_resources, temp_output_dir, monkeypatch):
        """Test that output file has correct mutual_trap format."""
        monkeypatch.chdir(test_resources)

        output_file = temp_output_dir / "forces.txt"

        test_args = [
            "generate_forces.py",
            "-o", str(output_file),
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        with open(output_file) as f:
            content = f.read()

        # Check for mutual_trap format elements
        assert "type = mutual_trap" in content, "Should contain mutual_trap type"
        assert "particle =" in content, "Should contain particle field"
        assert "ref_particle =" in content, "Should contain ref_particle field"
        assert "stiff =" in content, "Should contain stiff field"
        assert "r0 = 1.2" in content, "Should contain r0 = 1.2"
        assert "PBC=1" in content, "Should contain PBC=1"

    def test_main_quiet_mode(self, test_resources, temp_output_dir, monkeypatch):
        """Test main() with quiet mode."""
        monkeypatch.chdir(test_resources)

        output_file = temp_output_dir / "quiet_forces.txt"

        test_args = [
            "generate_forces.py",
            "-q",
            "-o", str(output_file),
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Quiet mode should still create output"

    def test_main_pairs_file_format(self, test_resources, temp_output_dir, monkeypatch):
        """Test that pairs file has correct format."""
        monkeypatch.chdir(test_resources)

        forces_file = temp_output_dir / "forces.txt"
        pairs_file = temp_output_dir / "pairs.txt"

        test_args = [
            "generate_forces.py",
            "-o", str(forces_file),
            "-f", str(pairs_file),
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        with open(pairs_file) as f:
            lines = f.readlines()

        # Each line should have two space-separated integers
        for line in lines:
            parts = line.strip().split()
            assert len(parts) == 2, "Each line should have 2 numbers"
            a, b = int(parts[0]), int(parts[1])
            assert a < b, "Smaller particle id should be first"


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_main_default_stiffness(self, test_resources, temp_output_dir, monkeypatch):
        """Test that default stiffness is 0.9."""
        monkeypatch.chdir(test_resources)

        output_file = temp_output_dir / "forces.txt"

        test_args = [
            "generate_forces.py",
            "-o", str(output_file),
            "input_rna",
            "minitraj.dat"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        with open(output_file) as f:
            content = f.read()

        assert "stiff = 0.9" in content, "Default stiffness should be 0.9"
