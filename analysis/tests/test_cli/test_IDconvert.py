"""
Tests for oxDNA_analysis_tools.IDconvert module.

Tests cover:
- IDconvert() API function: identity mapping, deletion, nicking, bijection
- Old-format topology support
- CLI argument parsing
- main() CLI entry point
- Circular strand handling
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
def greedy_steal_tops(tmp_path):
    """Topologies that expose the greedy set-based strand-matching bug.

    The reference has two *identical* strands (1 and 2) plus a third strand
    that differs by exactly one base at its 3' end.  The compare topology has
    only one copy of the all-A strand and one copy of the A+T strand.

    Reference:
        strand 1 (IDs  0-15): AAAAAAAAAAAAAAAA  ← correct match = cmp strand 1
        strand 2 (IDs 16-31): AAAAAAAAAAAAAAAA  ← duplicate; should be unmatched (deleted)
        strand 3 (IDs 32-47): AAAAAAAAAAAAAAAT  ← correct match = cmp strand 2

    Compare:
        strand 1 (IDs  0-15): AAAAAAAAAAAAAAAA
        strand 2 (IDs 16-31): AAAAAAAAAAAAAAAT

    The greedy phase-1 loop claims matches in ref-strand order.  After ref
    strand 1 takes cmp strand 1 (score=1.0), ref strand 2 (the duplicate)
    finds cmp strand 2 is the best *available* option (score≈0.94 ≥ 0.75) and
    claims it.  All compare strands are now exhausted before ref strand 3 is
    evaluated, so its nucleotides receive no match — even though cmp strand 2
    is a perfect sequence match for it.
    """
    ref_top = tmp_path / "ref.top"
    ref_top.write_text(
        "48 3 5->3\n"
        "AAAAAAAAAAAAAAAA id=1 type=DNA circular=false\n"  # IDs  0-15
        "AAAAAAAAAAAAAAAA id=2 type=DNA circular=false\n"  # IDs 16-31 (duplicate, deleted in cmp)
        "AAAAAAAAAAAAAAAT id=3 type=DNA circular=false\n"  # IDs 32-47
    )
    cmp_top = tmp_path / "cmp.top"
    cmp_top.write_text(
        "32 2 5->3\n"
        "AAAAAAAAAAAAAAAA id=1 type=DNA circular=false\n"  # IDs  0-15
        "AAAAAAAAAAAAAAAT id=2 type=DNA circular=false\n"  # IDs 16-31
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

    def test_greedy_set_misses_strand(self, greedy_steal_tops):
        """Duplicate ref strand steals the only viable cmp match for a distinct strand.

        ref strand 2 is an exact duplicate of ref strand 1 and has been deleted in the
        compare topology.  The greedy phase-1 loop processes ref strand 2 before ref
        strand 3 and assigns cmp strand 2 (score≈0.94, best available) to it, exhausting
        all compare strands.  ref strand 3 (IDs 32-47, sequence AAAAAAAAAAAAAAAT) is
        therefore left with no match, even though cmp strand 2 is a *perfect* match for
        it.  Phase 2 cannot recover because all cmp nucleotide IDs are already marked as
        used by the phase-1 assignments.

        All 16 nucleotides of ref strand 3 should be present in the returned mapping.
        """
        ref_top, cmp_top = greedy_steal_tops
        id_map = IDconvert(ref_top, cmp_top)

        # ref strand 3 (IDs 32-47, AAAAAAAAAAAAAAAT) is a perfect match for
        # cmp strand 2 (IDs 16-31, AAAAAAAAAAAAAAAT).  Every nucleotide must
        # appear in the mapping.
        for i in range(32, 48):
            assert i in id_map, (
                f"Nucleotide {i} (ref strand 3, 'AAAAAAAAAAAAAAAT') should map to "
                f"cmp strand 2 ('AAAAAAAAAAAAAAAT'), but the greedy 'used' set allowed "
                f"the duplicate ref strand 2 to steal that compare strand first."
            )

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

    def test_min_score_flag(self, two_strands_top, capsys):
        """--min-score flag is accepted and forwarded to IDconvert()."""
        with patch.object(sys, 'argv', [
            "IDconvert.py", two_strands_top, two_strands_top, "-q",
            "--min-score", "0.5"
        ]):
            main()
        tokens = capsys.readouterr().out.strip().split()
        assert len(tokens) == 16

    def test_min_block_flag(self, two_strands_top, capsys):
        """--min-block flag is accepted and forwarded to IDconvert()."""
        with patch.object(sys, 'argv', [
            "IDconvert.py", two_strands_top, two_strands_top, "-q",
            "--min-block", "4"
        ]):
            main()
        tokens = capsys.readouterr().out.strip().split()
        assert len(tokens) == 16


# =============================================================================
# Circular Strand Tests
# =============================================================================

@pytest.fixture
def circular_tops(tmp_path):
    """New-format topology with one circular strand (IDs 0-11, 12 nt)."""
    top = tmp_path / "circular.top"
    top.write_text(
        "12 1 5->3\n"
        "ACGTACGTACGT id=1 type=DNA circular=true\n"   # IDs 0-11
    )
    return str(top)


@pytest.fixture
def circular_nicked_tops(tmp_path):
    """
    Reference has a single 12-nt circular strand; comparison has it nicked
    into two 6-nt linear strands (simulating a nick introduced at position 6).
    """
    ref_top = tmp_path / "ref_circ.top"
    ref_top.write_text(
        "12 1 5->3\n"
        "ACGTACGTACGT id=1 type=DNA circular=true\n"   # IDs 0-11
    )
    cmp_top = tmp_path / "cmp_nicked.top"
    cmp_top.write_text(
        "12 2 5->3\n"
        "ACGTAC id=1 type=DNA circular=false\n"   # IDs 0-5  (5' half)
        "GTACGT id=2 type=DNA circular=false\n"   # IDs 6-11 (3' half)
    )
    return str(ref_top), str(cmp_top)


@pytest.fixture
def circular_rotated_tops(tmp_path):
    """
    Reference has a 12-nt circular strand; comparison has the same sequence
    rotated by 4 positions (simulating a different linearisation point).
    This documents the known limitation that a rotated circular strand may
    not achieve a full mapping via sequence similarity alone.
    """
    ref_top = tmp_path / "ref_rot.top"
    ref_top.write_text(
        "12 1 5->3\n"
        "ACGTACGTACGT id=1 type=DNA circular=true\n"    # IDs 0-11
    )
    cmp_top = tmp_path / "cmp_rot.top"
    cmp_top.write_text(
        "12 1 5->3\n"
        "ACGTACGTACGT id=1 type=DNA circular=false\n"   # IDs 0-11, same seq, linear
    )
    return str(ref_top), str(cmp_top)


@pytest.fixture
def circular_nicked_same_strand_tops(tmp_path):
    """
    Ref: 32-nt circular strand (IDs 0-31).
    Cmp: same strand nicked at position 8, stored as linear starting from ID 8
         (sequence rotated by 8 positions, IDs 0-31).
    Phase 1 score = 2*24/64 = 0.75 (exactly meets threshold).
    Wrap-around = 8 nt = min_block_size → phase 2 must recover them.
    """
    ref_top = tmp_path / "ref_circ32.top"
    ref_top.write_text(
        "32 1 5->3\n"
        "ACGTGAATCGAAATCCTGATACGTACGTGAAT id=1 type=DNA circular=true\n"
    )
    cmp_top = tmp_path / "cmp_nicked32.top"
    cmp_top.write_text(
        "32 1 5->3\n"
        "CGAAATCCTGATACGTACGTGAATACGTGAAT id=1 type=DNA circular=false\n"
    )
    return str(ref_top), str(cmp_top)


class TestCircularStrands:
    """Tests for circular-strand handling in IDconvert."""

    def test_circular_identity(self, circular_tops):
        """A circular strand mapped against itself gives a complete identity mapping."""
        id_map = IDconvert(circular_tops, circular_tops)
        assert len(id_map) == 12
        for ref_id, cmp_id in id_map.items():
            assert ref_id == cmp_id, f"ID {ref_id} should map to itself, got {cmp_id}"

    def test_circular_identity_bijection(self, circular_tops):
        """Identity mapping of a circular strand is injective (no duplicate targets)."""
        id_map = IDconvert(circular_tops, circular_tops)
        cmp_ids = list(id_map.values())
        assert len(cmp_ids) == len(set(cmp_ids))

    def test_nicked_circular_strand_recovered(self, circular_nicked_tops):
        """
        Phase 2 partially recovers nucleotides from a nicked circular strand.

        The reference has a 12-nt circular strand (ACGTACGTACGT); the
        comparison splits it into two 6-nt linear fragments (ACGTAC + GTACGT).
        Because the sequence is a repeating ACGT motif and the fragments are
        only 6 nt, SequenceMatcher cannot unambiguously resolve all 12
        positions — this is a known limitation for short circular strands with
        repetitive sequences.  The test documents the current behavior: at
        least 6 nucleotides are mapped, the result is injective, and all
        mapped targets lie within the 12-nt comparison topology.
        """
        ref_top, cmp_top = circular_nicked_tops
        id_map = IDconvert(ref_top, cmp_top, min_phase1_score=0.4, min_block_size=3)
        # At least the phase-1 matched fragment should appear in the mapping
        assert len(id_map) >= 6, (
            f"Expected at least 6 mapped nucleotides, got {len(id_map)}"
        )
        # Mapping must remain injective
        cmp_ids = list(id_map.values())
        assert len(cmp_ids) == len(set(cmp_ids)), "Mapping must be injective"
        # All mapped target IDs are within the 12-nt comparison topology
        assert all(0 <= cid < 12 for cid in cmp_ids), "All target IDs must be in range 0-11"

    def test_circular_nicked_same_strand(self, circular_nicked_same_strand_tops):
        """
        A circular strand nicked at a position other than its stored 5' end becomes a
        rotated linear strand.  Phase 2 must recover the wrap-around nucleotides that
        phase 1 misses due to SequenceMatcher's monotonicity constraint.
        All 32 nucleotides should be mapped, and the mapping must be injective.
        """
        ref_top, cmp_top = circular_nicked_same_strand_tops
        id_map = IDconvert(ref_top, cmp_top)
        assert len(id_map) == 32, (
            f"All 32 nucleotides should be mapped; got {len(id_map)}"
        )
        cmp_ids = list(id_map.values())
        assert len(cmp_ids) == len(set(cmp_ids)), "Mapping must be injective"

    def test_circular_vs_linear_rotation_documents_behavior(self, circular_rotated_tops):
        """
        Document current behavior when a circular strand is compared against a
        linear strand with the same sequence.

        A circular strand (ref) vs the same sequence as a linear strand (cmp)
        should achieve a high-scoring phase-1 match because the sequences are
        identical.  This test asserts that *at least* the overlapping bases are
        mapped, and records the known limitation that a purely rotated circular
        comparison may not map every base.
        """
        ref_top, cmp_top = circular_rotated_tops
        id_map = IDconvert(ref_top, cmp_top)
        # At minimum we expect some mapping — not an empty result
        assert len(id_map) > 0, "Expected at least some nucleotides to be mapped"
        # Mapping must be injective
        cmp_ids = list(id_map.values())
        assert len(cmp_ids) == len(set(cmp_ids)), "Mapping must be injective"
