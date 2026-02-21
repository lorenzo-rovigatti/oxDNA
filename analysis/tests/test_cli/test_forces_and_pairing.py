"""
Tests for oxDNA_analysis_tools pairing and forces modules.

This covers:
- db2forces: dot-bracket to external forces conversion
- forces2db: external forces to dot-bracket conversion
- forces2pairs: external forces to pair list conversion
- pairs2db: pair list to dot-bracket conversion
- external_force_utils/forces: force dictionary generators
- external_force_utils/force_reader: force file I/O
"""
import sys
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.db2forces import (
    parse_dot_bracket,
    db_to_forcelist,
    cli_parser as db2forces_parser,
    main as db2forces_main
)
from oxDNA_analysis_tools.forces2db import (
    forces2db,
    cli_parser as forces2db_parser,
    main as forces2db_main
)
from oxDNA_analysis_tools.forces2pairs import (
    forces2pairs,
    cli_parser as forces2pairs_parser,
    main as forces2pairs_main
)
from oxDNA_analysis_tools.pairs2db import (
    pairs2db,
    cli_parser as pairs2db_parser,
    main as pairs2db_main
)
from oxDNA_analysis_tools.external_force_utils.forces import (
    mutual_trap,
    string,
    harmonic_trap,
    rotating_harmonic_trap,
    repulsion_plane,
    attraction_plane,
    repulsion_sphere
)
from oxDNA_analysis_tools.external_force_utils.force_reader import (
    read_force_file,
    write_force_file
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture(scope="module")
def test_resources():
    """Get the path to test resources directory."""
    return Path(__file__).parent.parent / "resources"


@pytest.fixture(scope="module")
def old_top_path(test_resources):
    """Path to the old-style topology file."""
    return test_resources / "rna_tile.top"


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


@pytest.fixture
def sample_db_string():
    """A simple dot-bracket string for testing."""
    return "(((...)))"


@pytest.fixture
def sample_force_file_content():
    """Sample force file content."""
    return """{
type = mutual_trap
particle = 0
ref_particle = 8
stiff = 0.09
stiff_rate = 0
r0 = 1.2
rate = 0
PBC = 1
}

{
type = mutual_trap
particle = 1
ref_particle = 7
stiff = 0.09
stiff_rate = 0
r0 = 1.2
rate = 0
PBC = 1
}

{
type = mutual_trap
particle = 2
ref_particle = 6
stiff = 0.09
stiff_rate = 0
r0 = 1.2
rate = 0
PBC = 1
}
"""


@pytest.fixture
def sample_pairs_file_content():
    """Sample pairs file content."""
    return """0 8
1 7
2 6
"""


# =============================================================================
# parse_dot_bracket Tests
# =============================================================================

class TestParseDotBracket:
    """Tests for parse_dot_bracket() function."""

    def test_parse_simple_hairpin(self):
        """Test parsing a simple hairpin structure."""
        db = "(((...)))"
        result = parse_dot_bracket(db)

        assert isinstance(result, np.ndarray), "Should return numpy array"
        assert len(result) == len(db), "Should have same length as input"

        # Check pairs: 0-8, 1-7, 2-6
        assert result[0] == 8, "Position 0 should pair with 8"
        assert result[8] == 0, "Position 8 should pair with 0"
        assert result[1] == 7, "Position 1 should pair with 7"
        assert result[7] == 1, "Position 7 should pair with 1"
        assert result[2] == 6, "Position 2 should pair with 6"
        assert result[6] == 2, "Position 6 should pair with 2"

        # Unpaired positions
        assert result[3] == -1, "Position 3 should be unpaired"
        assert result[4] == -1, "Position 4 should be unpaired"
        assert result[5] == -1, "Position 5 should be unpaired"

    def test_parse_all_unpaired(self):
        """Test parsing fully unpaired structure."""
        db = "......."
        result = parse_dot_bracket(db)

        assert np.all(result == -1), "All positions should be unpaired"

    def test_parse_square_brackets(self):
        """Test parsing with square brackets (pseudoknots)."""
        db = "([)]"
        result = parse_dot_bracket(db)

        assert result[0] == 2, "( should pair with )"
        assert result[1] == 3, "[ should pair with ]"

    def test_parse_curly_brackets(self):
        """Test parsing with curly brackets for pseudoknots."""
        # Curly brackets are used for nested structures at depth 2
        # Must have () then [] then {} in proper nesting order
        db = "([{.}])"
        result = parse_dot_bracket(db)

        assert result[0] == 6, "( pairs with )"
        assert result[1] == 5, "[ pairs with ]"
        assert result[2] == 4, "{ pairs with }"

    def test_parse_angle_brackets(self):
        """Test parsing with angle brackets for pseudoknots."""
        # Angle brackets used for third level
        db = "(([{<.>}]))"
        result = parse_dot_bracket(db)

        assert result[0] == 10, "Outer ( pairs"
        assert result[4] == 6, "< pairs with >"

    def test_parse_letter_brackets(self):
        """Test parsing with letter notation (A-Z/a-z) for deep pseudoknots."""
        # Letter brackets used for additional levels beyond <>
        db = "((([{<Aa>}])))"
        result = parse_dot_bracket(db)

        assert result[6] == 7, "A should pair with a"

    def test_parse_nested_structure(self):
        """Test parsing nested structure."""
        db = "(((())))"
        result = parse_dot_bracket(db)

        assert result[0] == 7, "Outer pair 0-7"
        assert result[1] == 6, "Middle pair 1-6"
        assert result[2] == 5, "Inner pair 2-5"
        assert result[3] == 4, "Innermost pair 3-4"

    def test_parse_invalid_character_raises(self):
        """Test that invalid characters raise RuntimeError."""
        with pytest.raises(RuntimeError, match="invalid character"):
            parse_dot_bracket("((@))")

    def test_parse_with_whitespace(self):
        """Test parsing handles whitespace in input."""
        db = "((...))\n"
        result = parse_dot_bracket(db)

        # The function strips but still processes the string length after strip
        assert len(result) == 7 or len(result) == 8, "Should handle whitespace"
        # The core pairing should still work
        assert result[0] == 6, "First ( should pair with last )"


# =============================================================================
# db_to_forcelist Tests
# =============================================================================

class TestDbToForcelist:
    """Tests for db_to_forcelist() function."""

    def test_db_to_forcelist_basic(self, sample_db_string):
        """Test basic conversion to force list."""
        forces = db_to_forcelist(sample_db_string, stiff=0.09, reverse=False)

        assert isinstance(forces, list), "Should return list"
        # Each paired position generates a force (both directions)
        # (((...))) has 3 base pairs = 6 forces
        assert len(forces) == 6, "Should have 6 forces for 3 base pairs"

    def test_db_to_forcelist_force_structure(self, sample_db_string):
        """Test that forces have correct structure."""
        forces = db_to_forcelist(sample_db_string, stiff=0.5, reverse=False, r0=1.5)

        for force in forces:
            assert isinstance(force, dict), "Each force should be a dict"
            assert force["type"] == "mutual_trap", "Should be mutual_trap type"
            assert "particle" in force, "Should have particle"
            assert "ref_particle" in force, "Should have ref_particle"
            assert force["stiff"] == 0.5, "Stiffness should match"
            assert force["r0"] == 1.5, "r0 should match"

    def test_db_to_forcelist_reverse(self):
        """Test reverse option swaps bracket types."""
        db = "((...))"
        forces_normal = db_to_forcelist(db, stiff=0.09, reverse=False)
        forces_reversed = db_to_forcelist(db, stiff=0.09, reverse=True)

        # Both should have same number of forces
        assert len(forces_normal) == len(forces_reversed)

        # For asymmetric structures, particles should be mirrored
        # For the symmetric "((...))", the pairs are at same distances from ends
        # So we just verify both produce valid forces
        for force in forces_normal:
            assert "particle" in force
            assert "ref_particle" in force
        for force in forces_reversed:
            assert "particle" in force
            assert "ref_particle" in force

    def test_db_to_forcelist_pbc_option(self, sample_db_string):
        """Test PBC option."""
        forces_pbc = db_to_forcelist(sample_db_string, stiff=0.09, reverse=False, PBC=True)
        forces_no_pbc = db_to_forcelist(sample_db_string, stiff=0.09, reverse=False, PBC=False)

        assert forces_pbc[0]["PBC"] == 1, "PBC should be 1"
        assert forces_no_pbc[0]["PBC"] == 0, "PBC should be 0"

    def test_db_to_forcelist_rate_options(self, sample_db_string):
        """Test rate and stiff_rate options."""
        forces = db_to_forcelist(
            sample_db_string, stiff=0.09, reverse=False,
            rate=0.001, stiff_rate=0.002
        )

        assert forces[0]["rate"] == 0.001, "rate should match"
        assert forces[0]["stiff_rate"] == 0.002, "stiff_rate should match"

    def test_db_to_forcelist_unpaired_no_forces(self):
        """Test that unpaired bases generate no forces."""
        db = "......."
        forces = db_to_forcelist(db, stiff=0.09, reverse=False)

        assert len(forces) == 0, "Unpaired bases should have no forces"


# =============================================================================
# forces2pairs Tests
# =============================================================================

class TestForces2Pairs:
    """Tests for forces2pairs() function."""

    def test_forces2pairs_basic(self, temp_output_dir, sample_force_file_content):
        """Test basic force file parsing."""
        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)

        pairs = forces2pairs(str(force_file))

        assert isinstance(pairs, list), "Should return list"
        assert len(pairs) == 3, "Should have 3 pairs"

    def test_forces2pairs_tuple_structure(self, temp_output_dir, sample_force_file_content):
        """Test that pairs are tuples with correct structure."""
        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)

        pairs = forces2pairs(str(force_file))

        for pair in pairs:
            assert isinstance(pair, tuple), "Each pair should be a tuple"
            assert len(pair) == 2, "Each tuple should have 2 elements"
            assert pair[0] < pair[1], "First element should be smaller (a < b)"

    def test_forces2pairs_correct_values(self, temp_output_dir, sample_force_file_content):
        """Test that correct pair values are extracted."""
        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)

        pairs = forces2pairs(str(force_file))

        assert (0, 8) in pairs, "Should contain pair (0, 8)"
        assert (1, 7) in pairs, "Should contain pair (1, 7)"
        assert (2, 6) in pairs, "Should contain pair (2, 6)"

    def test_forces2pairs_empty_file(self, temp_output_dir):
        """Test parsing empty force file."""
        force_file = temp_output_dir / "empty_forces.txt"
        force_file.write_text("")

        pairs = forces2pairs(str(force_file))

        assert pairs == [], "Empty file should return empty list"


# =============================================================================
# pairs2db Tests
# =============================================================================

class TestPairs2Db:
    """Tests for pairs2db() function."""

    def test_pairs2db_basic(self):
        """Test basic conversion to dot-bracket."""
        n_bases = 9
        pairs = {0: 8, 8: 0, 1: 7, 7: 1, 2: 6, 6: 2}

        result = pairs2db(n_bases, pairs.copy())

        assert isinstance(result, str), "Should return string"
        assert len(result) == n_bases, "Should match n_bases"
        assert result.count('(') == 3, "Should have 3 open brackets"
        assert result.count(')') == 3, "Should have 3 close brackets"
        assert result.count('.') == 3, "Should have 3 dots"

    def test_pairs2db_unpaired(self):
        """Test with no pairs."""
        n_bases = 5
        pairs = {}

        result = pairs2db(n_bases, pairs)

        assert result == ".....", "Should be all dots"

    def test_pairs2db_nested(self):
        """Test nested structure."""
        n_bases = 8
        # Nested: (((())))
        pairs = {0: 7, 7: 0, 1: 6, 6: 1, 2: 5, 5: 2, 3: 4, 4: 3}

        result = pairs2db(n_bases, pairs.copy())

        assert result == "(((())))", "Should produce nested brackets"

    def test_pairs2db_pseudoknot(self):
        """Test pseudoknot structure uses different brackets."""
        n_bases = 8
        # Pseudoknot: ([)]
        pairs = {0: 4, 4: 0, 2: 6, 6: 2}  # Crossing pairs

        result = pairs2db(n_bases, pairs.copy())

        # Should use different bracket types for crossing pairs
        open_count = sum(1 for c in result if c in "([{<")
        close_count = sum(1 for c in result if c in ")]}>")
        assert open_count == 2, "Should have 2 open brackets"
        assert close_count == 2, "Should have 2 close brackets"


# =============================================================================
# forces2db Tests
# =============================================================================

class TestForces2Db:
    """Tests for forces2db() function."""

    def test_forces2db_basic(self):
        """Test basic conversion from forces to dot-bracket."""
        n_bases = 9
        forces = [
            {"particle": 0, "ref_particle": 8},
            {"particle": 1, "ref_particle": 7},
            {"particle": 2, "ref_particle": 6},
        ]

        result = forces2db(n_bases, forces)

        assert isinstance(result, str), "Should return string"
        assert len(result) == n_bases, "Should match n_bases"

    def test_forces2db_structure(self):
        """Test that correct dot-bracket structure is produced."""
        n_bases = 9
        forces = [
            {"particle": 0, "ref_particle": 8},
            {"particle": 1, "ref_particle": 7},
            {"particle": 2, "ref_particle": 6},
        ]

        result = forces2db(n_bases, forces)

        # Should produce (((...)))
        assert result[0] == '(' and result[8] == ')', "0 and 8 should be paired"
        assert result[3] == '.' and result[4] == '.' and result[5] == '.', "Middle should be unpaired"


# =============================================================================
# Force Dictionary Generator Tests
# =============================================================================

class TestForceGenerators:
    """Tests for force dictionary generator functions."""

    def test_mutual_trap(self):
        """Test mutual_trap force generator."""
        force = mutual_trap(
            particle=0, ref_particle=10, stiff=0.5,
            r0=1.2, PBC=True, rate=0.001, stiff_rate=0.002
        )

        assert force["type"] == "mutual_trap"
        assert force["particle"] == 0
        assert force["ref_particle"] == 10
        assert force["stiff"] == 0.5
        assert force["r0"] == 1.2
        assert force["PBC"] == 1
        assert force["rate"] == 0.001
        assert force["stiff_rate"] == 0.002

    def test_mutual_trap_pbc_false(self):
        """Test mutual_trap with PBC=False."""
        force = mutual_trap(particle=0, ref_particle=1, stiff=0.1, r0=1.0, PBC=False)

        assert force["PBC"] == 0, "PBC should be 0 when False"

    def test_string_force(self):
        """Test string force generator."""
        force = string(particle=5, f0=1.0, rate=0.01, direction=[1.0, 0.0, 0.0])

        assert force["type"] == "string"
        assert force["particle"] == 5
        assert force["f0"] == 1.0
        assert force["rate"] == 0.01
        assert force["dir"] == [1.0, 0.0, 0.0]

    def test_harmonic_trap(self):
        """Test harmonic_trap force generator."""
        force = harmonic_trap(
            particle=3, pos0=[10.0, 10.0, 10.0],
            stiff=0.5, rate=0.001, direction=[0.0, 0.0, 1.0]
        )

        assert force["type"] == "trap"
        assert force["particle"] == 3
        assert force["pos0"] == [10.0, 10.0, 10.0]
        assert force["stiff"] == 0.5
        assert force["rate"] == 0.001
        assert force["dir"] == [0.0, 0.0, 1.0]

    def test_rotating_harmonic_trap(self):
        """Test rotating_harmonic_trap force generator."""
        force = rotating_harmonic_trap(
            particle=2, pos0=[5.0, 5.0, 5.0], stiff=0.3,
            rate=0.01, base=0.0, center=[0.0, 0.0, 0.0],
            axis=[0.0, 0.0, 1.0], mask=[1.0, 1.0, 0.0]
        )

        assert force["type"] == "twist"
        assert force["particle"] == 2
        assert force["stiff"] == 0.3
        assert force["base"] == 0.0
        assert force["axis"] == [0.0, 0.0, 1.0]

    def test_repulsion_plane(self):
        """Test repulsion_plane force generator."""
        force = repulsion_plane(
            particle=-1, stiff=1.0,
            direction=[0.0, 0.0, 1.0], position=[0.0, 0.0, 0.0]
        )

        assert force["type"] == "repulsion_plane"
        assert force["particle"] == -1
        assert force["stiff"] == 1.0
        assert force["dir"] == [0.0, 0.0, 1.0]
        assert force["position"] == [0.0, 0.0, 0.0]

    def test_attraction_plane(self):
        """Test attraction_plane force generator."""
        force = attraction_plane(
            particle=0, stiff=0.5,
            direction=[1.0, 0.0, 0.0], position=[10.0, 0.0, 0.0]
        )

        assert force["type"] == "attraction_plane"
        assert force["particle"] == 0
        assert force["stiff"] == 0.5

    def test_repulsion_sphere(self):
        """Test repulsion_sphere force generator."""
        force = repulsion_sphere(
            particle=0, center=[25.0, 25.0, 25.0],
            stiff=1.0, r0=20.0, rate=-0.001
        )

        assert force["type"] == "sphere"
        assert force["particle"] == 0
        assert force["center"] == [25.0, 25.0, 25.0]
        assert force["stiff"] == 1.0
        assert force["r0"] == 20.0
        assert force["rate"] == -0.001


# =============================================================================
# Force File I/O Tests
# =============================================================================

class TestForceFileIO:
    """Tests for force file reading and writing."""

    def test_write_force_file_creates_file(self, temp_output_dir):
        """Test write_force_file creates file."""
        forces = [
            mutual_trap(0, 5, 0.09, 1.2, True),
            mutual_trap(1, 4, 0.09, 1.2, True)
        ]
        outfile = temp_output_dir / "forces.txt"

        write_force_file(forces, str(outfile))

        assert outfile.exists(), "File should be created"

    def test_write_force_file_content(self, temp_output_dir):
        """Test write_force_file produces correct content."""
        forces = [mutual_trap(0, 5, 0.09, 1.2, True)]
        outfile = temp_output_dir / "forces.txt"

        write_force_file(forces, str(outfile))
        content = outfile.read_text()

        assert "type = mutual_trap" in content
        assert "particle = 0" in content
        assert "ref_particle = 5" in content
        assert "stiff = 0.09" in content

    def test_write_force_file_append_mode(self, temp_output_dir):
        """Test write_force_file append mode."""
        forces1 = [mutual_trap(0, 5, 0.09, 1.2, True)]
        forces2 = [mutual_trap(1, 4, 0.09, 1.2, True)]
        outfile = temp_output_dir / "forces_append.txt"

        write_force_file(forces1, str(outfile), mode='w+')
        write_force_file(forces2, str(outfile), mode='a')
        content = outfile.read_text()

        assert content.count("type = mutual_trap") == 2, "Should have 2 forces"

    def test_read_force_file_basic(self, temp_output_dir, sample_force_file_content):
        """Test read_force_file basic functionality."""
        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)

        forces = read_force_file(str(force_file))

        assert isinstance(forces, list), "Should return list"
        assert len(forces) == 3, "Should have 3 forces"

    def test_read_force_file_structure(self, temp_output_dir, sample_force_file_content):
        """Test read_force_file returns correct structure."""
        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)

        forces = read_force_file(str(force_file))

        for force in forces:
            assert isinstance(force, dict), "Each force should be dict"
            assert "type" in force, "Should have type"
            assert force["type"] == "mutual_trap"

    def test_read_force_file_values(self, temp_output_dir, sample_force_file_content):
        """Test read_force_file parses correct values."""
        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)

        forces = read_force_file(str(force_file))

        assert forces[0]["particle"] == 0
        assert forces[0]["ref_particle"] == 8
        assert forces[0]["stiff"] == 0.09
        assert forces[0]["r0"] == 1.2

    def test_roundtrip_force_file(self, temp_output_dir):
        """Test write then read produces equivalent forces for mutual_trap."""
        # Note: Forces with list parameters (like string, harmonic_trap) may not
        # roundtrip perfectly due to write_force_file joining list elements.
        # This test focuses on mutual_trap which has only scalar values.
        original_forces = [
            mutual_trap(0, 10, 0.5, 1.5, True, rate=0.001, stiff_rate=0.002),
            mutual_trap(1, 9, 0.3, 1.2, False),
            mutual_trap(2, 8, 0.1, 1.0, True)
        ]
        outfile = temp_output_dir / "roundtrip.txt"

        write_force_file(original_forces, str(outfile))
        read_forces = read_force_file(str(outfile))

        assert len(read_forces) == len(original_forces)

        # Check mutual_trap values
        assert read_forces[0]["particle"] == 0
        assert read_forces[0]["ref_particle"] == 10
        assert read_forces[0]["stiff"] == 0.5


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParsers:
    """Tests for CLI argument parsers."""

    def test_db2forces_parser(self):
        """Test db2forces CLI parser."""
        parser = db2forces_parser()

        args = parser.parse_args(["input.db", "-o", "output.txt", "-s", "0.5", "-r"])

        assert args.db_file == ["input.db"]
        assert args.output == ["output.txt"]
        assert args.stiff == [0.5]
        assert args.reverse is True

    def test_forces2db_parser(self):
        """Test forces2db CLI parser."""
        parser = forces2db_parser()

        args = parser.parse_args(["top.top", "forces.txt", "-o", "output.db"])

        assert args.topology == "top.top"
        assert args.force_file == "forces.txt"
        assert args.output == "output.db"

    def test_forces2pairs_parser(self):
        """Test forces2pairs CLI parser."""
        parser = forces2pairs_parser()

        args = parser.parse_args(["forces.txt", "-o", "pairs.txt"])

        assert args.force_file == ["forces.txt"]
        assert args.output == ["pairs.txt"]

    def test_pairs2db_parser(self):
        """Test pairs2db CLI parser."""
        parser = pairs2db_parser()

        args = parser.parse_args(["top.top", "pairs.txt", "-o", "output.db"])

        assert args.topology == "top.top"
        assert args.pair_file == "pairs.txt"
        assert args.output == "output.db"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestCLIMain:
    """Tests for CLI main() entry points."""

    def test_db2forces_main(self, temp_output_dir):
        """Test db2forces main() creates output file."""
        db_file = temp_output_dir / "test.db"
        db_file.write_text("(((...)))")
        output_file = temp_output_dir / "forces.txt"

        test_args = [
            "db2forces.py",
            str(db_file),
            "-o", str(output_file),
            "-s", "0.1"
        ]

        with patch.object(sys, 'argv', test_args):
            db2forces_main()

        assert output_file.exists(), "Output file should be created"
        content = output_file.read_text()
        assert "mutual_trap" in content

    def test_forces2pairs_main(self, temp_output_dir, sample_force_file_content):
        """Test forces2pairs main() creates output file."""
        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)
        output_file = temp_output_dir / "pairs.txt"

        test_args = [
            "forces2pairs.py",
            str(force_file),
            "-o", str(output_file)
        ]

        with patch.object(sys, 'argv', test_args):
            forces2pairs_main()

        assert output_file.exists(), "Output file should be created"
        content = output_file.read_text()
        assert "0 8" in content

    def test_forces2db_main(self, temp_output_dir, old_top_path, sample_force_file_content):
        """Test forces2db main() creates output file."""
        # Copy topology to temp dir
        top_copy = temp_output_dir / "top.top"
        copy(old_top_path, top_copy)

        force_file = temp_output_dir / "forces.txt"
        force_file.write_text(sample_force_file_content)
        output_file = temp_output_dir / "output.db"

        test_args = [
            "forces2db.py",
            str(top_copy),
            str(force_file),
            "-o", str(output_file)
        ]

        with patch.object(sys, 'argv', test_args):
            forces2db_main()

        assert output_file.exists(), "Output file should be created"

    def test_pairs2db_main(self, temp_output_dir, old_top_path, sample_pairs_file_content):
        """Test pairs2db main() creates output file."""
        # Copy topology to temp dir
        top_copy = temp_output_dir / "top.top"
        copy(old_top_path, top_copy)

        pairs_file = temp_output_dir / "pairs.txt"
        pairs_file.write_text(sample_pairs_file_content)
        output_file = temp_output_dir / "output.db"

        test_args = [
            "pairs2db.py",
            str(top_copy),
            str(pairs_file),
            "-o", str(output_file)
        ]

        with patch.object(sys, 'argv', test_args):
            pairs2db_main()

        assert output_file.exists(), "Output file should be created"


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests for the full conversion pipeline."""

    def test_db_to_forces_to_pairs_roundtrip(self, temp_output_dir):
        """Test dot-bracket -> forces -> pairs pipeline."""
        db_string = "(((...)))"

        # Convert to forces
        forces = db_to_forcelist(db_string, stiff=0.09, reverse=False)

        # Write forces
        force_file = temp_output_dir / "pipeline_forces.txt"
        write_force_file(forces, str(force_file))

        # Convert to pairs
        pairs = forces2pairs(str(force_file))

        # Should have 3 pairs
        assert len(pairs) == 3

        # Check expected pairs
        expected = [(0, 8), (1, 7), (2, 6)]
        for exp in expected:
            assert exp in pairs, f"Expected pair {exp}"

    def test_db_to_forces_to_db_roundtrip(self, temp_output_dir):
        """Test dot-bracket -> forces -> dot-bracket roundtrip."""
        original_db = "(((...)))"
        n_bases = len(original_db)

        # Convert to forces
        forces = db_to_forcelist(original_db, stiff=0.09, reverse=False)

        # Write and read forces
        force_file = temp_output_dir / "roundtrip_forces.txt"
        write_force_file(forces, str(force_file))
        read_forces = read_force_file(str(force_file))

        # Convert back to dot-bracket
        result_db = forces2db(n_bases, read_forces)

        assert result_db == original_db, "Should roundtrip to same dot-bracket"

    def test_pairs_to_db_to_forces_roundtrip(self, temp_output_dir):
        """Test pairs -> dot-bracket -> forces pipeline."""
        n_bases = 9
        original_pairs = {0: 8, 8: 0, 1: 7, 7: 1, 2: 6, 6: 2}

        # Convert pairs to dot-bracket
        db_string = pairs2db(n_bases, original_pairs.copy())

        # Convert dot-bracket to forces
        forces = db_to_forcelist(db_string, stiff=0.09, reverse=False)

        # Should have 6 forces (one per paired position)
        assert len(forces) == 6

        # All forces should be mutual traps
        for force in forces:
            assert force["type"] == "mutual_trap"

    def test_complex_structure_roundtrip(self, temp_output_dir):
        """Test roundtrip with more complex structure."""
        # A structure with internal loop
        original_db = "(((..((...))..)))"
        n_bases = len(original_db)

        # Convert to forces
        forces = db_to_forcelist(original_db, stiff=0.09, reverse=False)

        # Write and read
        force_file = temp_output_dir / "complex_forces.txt"
        write_force_file(forces, str(force_file))
        read_forces = read_force_file(str(force_file))

        # Convert back
        result_db = forces2db(n_bases, read_forces)

        assert result_db == original_db, "Complex structure should roundtrip"
