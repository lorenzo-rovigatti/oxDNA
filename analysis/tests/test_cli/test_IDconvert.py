"""
Tests for oxDNA_analysis_tools.IDconvert module.

Tests cover:
- IDconvert() API function: identity mapping, deletion, nicking, bijection
- Old-format topology support
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from oxDNA_analysis_tools.IDconvert import IDconvert, cli_parser, main


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def test_resources():
    """Get the path to test resources directory."""
    return Path(__file__).parent.parent / "resources"


@pytest.fixture(scope="module")
def rna_tile_top(test_resources):
    """Path to the old-format RNA topology in test resources."""
    return str(test_resources / "rna_tile.top")


@pytest.fixture
def two_strands_top(tmp_path):
    """A simple 2-strand new-format topology for identity-mapping tests."""
    top = tmp_path / "two_strands.top"
    top.write_text(
        "16 2 5->3\n"
        "ACGTACGT id=1 type=DNA circular=false\n"   # IDs 0-7
        "GCTAGCTA id=2 type=DNA circular=false\n"   # IDs 8-15
    )
    return str(top)


@pytest.fixture
def deletion_tops(tmp_path):
    """Reference and comparison topologies where the middle strand is deleted."""
    ref_top = tmp_path / "ref.top"
    ref_top.write_text(
        "48 3 5->3\n"
        "ACGTACGTACGTACGT id=1 type=DNA circular=false\n"   # IDs 0-15
        "TTTTTTTTTTTTTTTT id=2 type=DNA circular=false\n"   # IDs 16-31 (deleted)
        "GCTAGCTAGCTAGCTA id=3 type=DNA circular=false\n"   # IDs 32-47
    )
    cmp_top = tmp_path / "cmp.top"
    cmp_top.write_text(
        "32 2 5->3\n"
        "ACGTACGTACGTACGT id=1 type=DNA circular=false\n"   # IDs 0-15
        "GCTAGCTAGCTAGCTA id=2 type=DNA circular=false\n"   # IDs 16-31
    )
    return str(ref_top), str(cmp_top)


@pytest.fixture
def nicking_tops(tmp_path):
    """Reference with one 24-nt strand, comparison where it is nicked into 16+8 nt."""
    # Sequences are chosen so that the 8-nt suffix is unique within the 24-nt ref,
    # ensuring SequenceMatcher finds the correct matching block in phase 2.
    ref_top = tmp_path / "ref.top"
    ref_top.write_text(
        "24 1 5->3\n"
        "ACGTGAATCGAAATCCTGATCGAC id=1 type=DNA circular=false\n"   # IDs 0-23
    )
    cmp_top = tmp_path / "cmp.top"
    cmp_top.write_text(
        "24 2 5->3\n"
        "ACGTGAATCGAAATCC id=1 type=DNA circular=false\n"   # IDs 0-15  (5' fragment)
        "TGATCGAC id=2 type=DNA circular=false\n"           # IDs 16-23 (3' fragment)
    )
    return str(ref_top), str(cmp_top)


# =============================================================================
# API Tests - IDconvert() function
# =============================================================================

class TestIDconvertFunction:
    """Tests for the IDconvert() API function."""

    def test_returns_dict(self, two_strands_top):
        """IDconvert() returns a dict."""
        result = IDconvert(two_strands_top, two_strands_top)
        assert isinstance(result, dict)

    def test_identity_mapping(self, two_strands_top):
        """Same topology as ref and cmp: every nucleotide ID maps to itself."""
        id_map = IDconvert(two_strands_top, two_strands_top)
        assert len(id_map) == 16
        for ref_id, cmp_id in id_map.items():
            assert ref_id == cmp_id, f"ID {ref_id} should map to itself, got {cmp_id}"

    def test_old_format_identity(self, rna_tile_top):
        """Old-format topology mapped against itself gives a complete identity mapping."""
        id_map = IDconvert(rna_tile_top, rna_tile_top)
        assert len(id_map) == 132   # rna_tile.top has 132 nucleotides
        for ref_id, cmp_id in id_map.items():
            assert ref_id == cmp_id, f"ID {ref_id} should map to itself, got {cmp_id}"

    def test_deletion(self, deletion_tops):
        """Nucleotides belonging to a deleted strand are absent from the returned dict."""
        ref_top, cmp_top = deletion_tops
        id_map = IDconvert(ref_top, cmp_top)

        # Strand 1 (IDs 0-15): present in both → should be mapped
        for i in range(16):
            assert i in id_map, f"Nucleotide {i} (strand 1) should be mapped"

        # Strand 2 (IDs 16-31): deleted → should be absent
        for i in range(16, 32):
            assert i not in id_map, f"Nucleotide {i} (deleted strand) should be absent"

        # Strand 3 (IDs 32-47): present in both → should be mapped
        for i in range(32, 48):
            assert i in id_map, f"Nucleotide {i} (strand 3) should be mapped"

    def test_bijection(self, deletion_tops):
        """Each compare nucleotide ID appears at most once as a mapping target."""
        ref_top, cmp_top = deletion_tops
        id_map = IDconvert(ref_top, cmp_top)
        cmp_ids = list(id_map.values())
        assert len(cmp_ids) == len(set(cmp_ids)), \
            "Mapping must be injective — duplicate compare IDs found"

    def test_nicking(self, nicking_tops):
        """All nucleotides of a nicked strand are mapped across the nick site."""
        ref_top, cmp_top = nicking_tops
        id_map = IDconvert(ref_top, cmp_top)

        # All 24 ref nucleotides should be accounted for (nothing deleted)
        assert len(id_map) == 24, f"Expected 24 mapped nucleotides, got {len(id_map)}"

        # The nicking preserves sequence order so IDs map to themselves
        for ref_id, cmp_id in id_map.items():
            assert ref_id == cmp_id, \
                f"Nucleotide {ref_id} should map to {ref_id}, got {cmp_id}"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_requires_two_topology_files(self):
        """Parser exits when fewer than two positional arguments are given."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])
        with pytest.raises(SystemExit):
            parser.parse_args(["ref.top"])

    def test_positional_args(self):
        """Parser correctly assigns reference and compare positional arguments."""
        parser = cli_parser()
        args = parser.parse_args(["ref.top", "cmp.top"])
        assert args.reference == "ref.top"
        assert args.compare == "cmp.top"

    def test_index_file_optional(self):
        """Index file flag is optional and defaults to None."""
        parser = cli_parser()
        args = parser.parse_args(["ref.top", "cmp.top"])
        assert args.index_file is None

    def test_index_file_flag(self):
        """Parser accepts -i/--index flag with a file path."""
        parser = cli_parser()
        args = parser.parse_args(["ref.top", "cmp.top", "-i", "ids.txt"])
        assert args.index_file == ["ids.txt"]

    def test_quiet_flag(self):
        """Parser accepts -q/--quiet flag."""
        parser = cli_parser()
        args = parser.parse_args(["ref.top", "cmp.top", "-q"])
        assert args.quiet is True

    def test_quiet_default(self):
        """Quiet defaults to False."""
        parser = cli_parser()
        args = parser.parse_args(["ref.top", "cmp.top"])
        assert args.quiet is False


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_outputs_space_separated_ids(self, two_strands_top, capsys):
        """main() prints one space-separated token per reference nucleotide."""
        with patch.object(sys, 'argv',
                          ["IDconvert.py", two_strands_top, two_strands_top, "-q"]):
            main()
        tokens = capsys.readouterr().out.strip().split()
        assert len(tokens) == 16   # 16 nucleotides in two_strands_top

    def test_main_identity_values(self, two_strands_top, capsys):
        """main() outputs the correct translated IDs for an identity mapping."""
        with patch.object(sys, 'argv',
                          ["IDconvert.py", two_strands_top, two_strands_top, "-q"]):
            main()
        tokens = capsys.readouterr().out.strip().split()
        assert tokens == [str(i) for i in range(16)]

    def test_main_deleted_shown_as_DELETED(self, deletion_tops, capsys):
        """main() outputs 'DELETED' for nucleotides absent from the compare topology."""
        ref_top, cmp_top = deletion_tops
        with patch.object(sys, 'argv', ["IDconvert.py", ref_top, cmp_top, "-q"]):
            main()
        tokens = capsys.readouterr().out.strip().split()
        assert len(tokens) == 48                                # one token per ref nucleotide
        assert tokens[16:32] == ["DELETED"] * 16               # strand 2 was deleted

    def test_main_with_index_file(self, two_strands_top, tmp_path, capsys):
        """main() with -i restricts output to the queried IDs."""
        index_file = tmp_path / "ids.txt"
        index_file.write_text("0 1 2 3 4")
        with patch.object(sys, 'argv', [
            "IDconvert.py", two_strands_top, two_strands_top, "-q",
            "-i", str(index_file)
        ]):
            main()
        tokens = capsys.readouterr().out.strip().split()
        assert tokens == ["0", "1", "2", "3", "4"]

    def test_main_with_comma_separated_index_file(self, two_strands_top, tmp_path, capsys):
        """main() accepts comma-separated index files as well as space-separated."""
        index_file = tmp_path / "ids.txt"
        index_file.write_text("0,1,2,3,4")
        with patch.object(sys, 'argv', [
            "IDconvert.py", two_strands_top, two_strands_top, "-q",
            "-i", str(index_file)
        ]):
            main()
        tokens = capsys.readouterr().out.strip().split()
        assert tokens == ["0", "1", "2", "3", "4"]
