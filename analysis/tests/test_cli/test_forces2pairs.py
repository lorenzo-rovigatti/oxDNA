"""
Tests for oxDNA_analysis_tools.forces2pairs module.

Tests cover:
- forces2pairs() API function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from oxDNA_analysis_tools.forces2pairs import forces2pairs, cli_parser, main


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


@pytest.fixture
def sample_force_file(temp_output_dir):
    """Create a sample force file for testing."""
    force_content = """{
type = mutual_trap
particle = 0
stiff = 0.9
r0 = 1.2
ref_particle = 10
PBC=1
}
{
type = mutual_trap
particle = 10
stiff = 0.9
r0 = 1.2
ref_particle = 0
PBC=1
}
{
type = mutual_trap
particle = 1
stiff = 0.9
r0 = 1.2
ref_particle = 11
PBC=1
}
{
type = mutual_trap
particle = 11
stiff = 0.9
r0 = 1.2
ref_particle = 1
PBC=1
}
{
type = mutual_trap
particle = 5
stiff = 0.9
r0 = 1.2
ref_particle = 15
PBC=1
}
{
type = mutual_trap
particle = 15
stiff = 0.9
r0 = 1.2
ref_particle = 5
PBC=1
}
"""
    force_file = temp_output_dir / "forces.txt"
    force_file.write_text(force_content)
    return force_file


@pytest.fixture
def empty_force_file(temp_output_dir):
    """Create an empty force file."""
    force_file = temp_output_dir / "empty_forces.txt"
    force_file.write_text("")
    return force_file


# =============================================================================
# API Tests - forces2pairs() function
# =============================================================================

class TestForces2PairsFunction:
    """Tests for the forces2pairs() API function."""

    def test_forces2pairs_returns_list(self, sample_force_file):
        """Test forces2pairs() returns a list."""
        result = forces2pairs(str(sample_force_file))
        assert isinstance(result, list), "Should return list"

    def test_forces2pairs_returns_tuples(self, sample_force_file):
        """Test forces2pairs() returns list of tuples."""
        result = forces2pairs(str(sample_force_file))
        for item in result:
            assert isinstance(item, tuple), "Each item should be tuple"
            assert len(item) == 2, "Each tuple should have 2 elements"

    def test_forces2pairs_correct_pairs(self, sample_force_file):
        """Test forces2pairs() extracts correct pairs."""
        result = forces2pairs(str(sample_force_file))

        # Expected pairs (smaller id first)
        expected = [(0, 10), (1, 11), (5, 15)]

        assert len(result) == 3, "Should find 3 pairs"
        for pair in expected:
            assert pair in result, f"Pair {pair} should be found"

    def test_forces2pairs_empty_file(self, empty_force_file):
        """Test forces2pairs() with empty file."""
        result = forces2pairs(str(empty_force_file))
        assert result == [], "Empty file should return empty list"

    def test_forces2pairs_only_smaller_first(self, sample_force_file):
        """Test that pairs have smaller particle id first."""
        result = forces2pairs(str(sample_force_file))
        for a, b in result:
            assert a < b, f"Pair ({a}, {b}) should have smaller id first"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_force_file(self):
        """Test that parser requires force file."""
        parser = cli_parser()

        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options."""
        parser = cli_parser()

        args = parser.parse_args([
            "-o", "pairs.txt",
            "-q",
            "forces.txt"
        ])

        assert args.force_file == ["forces.txt"], "Force file not parsed"
        assert args.output == ["pairs.txt"], "Output option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

    def test_parser_defaults(self):
        """Test parser default values."""
        parser = cli_parser()

        args = parser.parse_args(["forces.txt"])

        assert args.force_file == ["forces.txt"]
        assert args.output is None
        assert args.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_output_file(self, sample_force_file, temp_output_dir, monkeypatch):
        """Test main() creates pairs output file."""
        monkeypatch.chdir(temp_output_dir)

        output_file = temp_output_dir / "pairs.txt"

        test_args = [
            "forces2pairs.py",
            "-o", str(output_file),
            str(sample_force_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output file should be created"

    def test_main_default_output(self, sample_force_file, temp_output_dir, monkeypatch):
        """Test main() uses default output filename."""
        monkeypatch.chdir(temp_output_dir)

        test_args = [
            "forces2pairs.py",
            str(sample_force_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        default_output = temp_output_dir / "pairs.txt"
        assert default_output.exists(), "Default output 'pairs.txt' should be created"

    def test_main_output_format(self, sample_force_file, temp_output_dir, monkeypatch):
        """Test that output file has correct format."""
        monkeypatch.chdir(temp_output_dir)

        output_file = temp_output_dir / "pairs.txt"

        test_args = [
            "forces2pairs.py",
            "-o", str(output_file),
            str(sample_force_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        with open(output_file) as f:
            lines = f.readlines()

        assert len(lines) == 3, "Should have 3 pairs"
        for line in lines:
            parts = line.strip().split()
            assert len(parts) == 2, "Each line should have 2 numbers"
            assert int(parts[0]) < int(parts[1]), "Smaller id should be first"

    def test_main_quiet_mode(self, sample_force_file, temp_output_dir, monkeypatch):
        """Test main() with quiet mode."""
        monkeypatch.chdir(temp_output_dir)

        output_file = temp_output_dir / "quiet_pairs.txt"

        test_args = [
            "forces2pairs.py",
            "-q",
            "-o", str(output_file),
            str(sample_force_file)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Quiet mode should still create output"


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_forces2pairs_single_pair(self, temp_output_dir):
        """Test forces2pairs with single pair."""
        force_content = """{
type = mutual_trap
particle = 5
stiff = 0.9
r0 = 1.2
ref_particle = 25
PBC=1
}
{
type = mutual_trap
particle = 25
stiff = 0.9
r0 = 1.2
ref_particle = 5
PBC=1
}
"""
        force_file = temp_output_dir / "single_pair.txt"
        force_file.write_text(force_content)

        result = forces2pairs(str(force_file))

        assert len(result) == 1, "Should find 1 pair"
        assert result[0] == (5, 25), "Pair should be (5, 25)"

    def test_forces2pairs_non_mutual_trap(self, temp_output_dir):
        """Test forces2pairs ignores non-mutual_trap forces."""
        force_content = """{
type = com
stiff = 1.0
r0 = 5
com_list = 0, 1, 2, 3
ref_list = 10, 11, 12, 13
}
{
type = mutual_trap
particle = 0
stiff = 0.9
r0 = 1.2
ref_particle = 10
PBC=1
}
{
type = mutual_trap
particle = 10
stiff = 0.9
r0 = 1.2
ref_particle = 0
PBC=1
}
"""
        force_file = temp_output_dir / "mixed_forces.txt"
        force_file.write_text(force_content)

        result = forces2pairs(str(force_file))

        # Should only find mutual_trap pairs
        assert len(result) == 1, "Should find 1 mutual_trap pair"

    def test_forces2pairs_large_particle_ids(self, temp_output_dir):
        """Test forces2pairs with large particle IDs."""
        force_content = """{
type = mutual_trap
particle = 10000
stiff = 0.9
r0 = 1.2
ref_particle = 20000
PBC=1
}
{
type = mutual_trap
particle = 20000
stiff = 0.9
r0 = 1.2
ref_particle = 10000
PBC=1
}
"""
        force_file = temp_output_dir / "large_ids.txt"
        force_file.write_text(force_content)

        result = forces2pairs(str(force_file))

        assert len(result) == 1, "Should find 1 pair"
        assert result[0] == (10000, 20000), "Pair should be (10000, 20000)"
