"""
Tests for oxDNA_analysis_tools.UTILS.RyeReader module.

This is the core file reading/writing utility used by all analysis scripts.
Tests cover comprehensive type coverage and correctness for:
- Topology parsing (old and new formats)
- Trajectory parsing and indexing
- Configuration reading and writing
- Strand/System object creation
- Input file parameter extraction
- Periodic boundary condition utilities
"""
import json
import os
import pickle
import sys
from io import BytesIO
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.UTILS.RyeReader import (
    Chunker,
    linear_read,
    get_confs,
    get_top_info,
    get_top_info_from_traj,
    get_traj_info,
    describe,
    strand_describe,
    get_input_parameter,
    inbox,
    write_conf,
    conf_to_str,
    write_top,
    get_top_string,
)
from oxDNA_analysis_tools.UTILS.data_structures import (
    Chunk,
    ConfInfo,
    TrajInfo,
    TopInfo,
    Configuration,
    System,
    Strand,
    Monomer,
)


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
def old_top_path(test_resources):
    """Path to the old-style topology file."""
    return test_resources / "rna_tile.top"


@pytest.fixture(scope="module")
def input_file_path(test_resources):
    """Path to the input file."""
    return test_resources / "input_rna"


@pytest.fixture(scope="module")
def trajectory_info(mini_traj_path, old_top_path):
    """Get topology and trajectory info for the mini trajectory."""
    top_info, traj_info = describe(str(old_top_path), str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


@pytest.fixture
def sample_configuration():
    """Create a sample configuration for testing."""
    return Configuration(
        time=100,
        box=np.array([50.0, 50.0, 50.0]),
        energy=np.array([0.0, 0.0, 0.0]),
        positions=np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]),
        a1s=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        a3s=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
    )


@pytest.fixture
def new_style_top_content():
    """Content for a new-style topology file."""
    return "6 2 5->3\nACG type=DNA circular=False\nCGU type=RNA circular=False\n"


@pytest.fixture
def old_style_top_content():
    """Content for an old-style topology file."""
    return """6 2
1 A -1 1
1 C 0 2
1 G 1 -1
2 C -1 4
2 G 3 5
2 T 4 -1
"""


# =============================================================================
# Chunker Tests
# =============================================================================

class TestChunker:
    """Tests for the Chunker generator function."""

    def test_chunker_basic_behavior(self, temp_output_dir):
        """Test Chunker yields chunks with correct metadata."""
        test_file = temp_output_dir / "test_chunk.dat"
        content = b"0123456789" * 100  # 1000 bytes
        test_file.write_bytes(content)

        fsize = os.stat(test_file).st_size
        with open(test_file, 'rb') as f:
            chunks = list(Chunker(f, fsize, size=250))

        # Should have 4 chunks for 1000 bytes with 250 byte chunks
        assert len(chunks) == 4, "Should have 4 chunks"

        # Check first chunk
        assert chunks[0].offset == 0, "First chunk offset should be 0"
        assert chunks[0].is_last is False, "First chunk should not be last"
        assert chunks[0].file_size == fsize, "File size should be correct"

        # Check last chunk offset
        assert chunks[-1].offset == 750, "Last chunk offset should be 750"
        # Note: is_last is True only when chunk extends beyond file size
        # For exact-fit chunks, is_last may be False on last chunk

    def test_chunker_small_file(self, temp_output_dir):
        """Test Chunker with file smaller than chunk size."""
        test_file = temp_output_dir / "small.dat"
        content = b"small content"
        test_file.write_bytes(content)

        fsize = os.stat(test_file).st_size
        with open(test_file, 'rb') as f:
            chunks = list(Chunker(f, fsize, size=1000000))

        assert len(chunks) == 1, "Small file should yield one chunk"
        assert chunks[0].is_last is True, "Single chunk should be last"
        assert chunks[0].block == content, "Content should match"

    def test_chunker_returns_chunk_type(self, temp_output_dir):
        """Test that Chunker returns Chunk dataclass instances."""
        test_file = temp_output_dir / "type_test.dat"
        test_file.write_bytes(b"test content")

        fsize = os.stat(test_file).st_size
        with open(test_file, 'rb') as f:
            chunk = next(Chunker(f, fsize, size=100))

        assert isinstance(chunk, Chunk), "Should return Chunk dataclass"
        assert hasattr(chunk, 'block'), "Chunk should have block attribute"
        assert hasattr(chunk, 'offset'), "Chunk should have offset attribute"
        assert hasattr(chunk, 'is_last'), "Chunk should have is_last attribute"


# =============================================================================
# Topology Info Tests
# =============================================================================

class TestGetTopInfo:
    """Tests for get_top_info() function."""

    def test_get_top_info_old_format(self, old_top_path):
        """Test parsing old-style topology file."""
        top_info = get_top_info(str(old_top_path))

        assert isinstance(top_info, TopInfo), "Should return TopInfo"
        assert top_info.nbases == 132, "Should have 132 bases"
        assert top_info.path == str(old_top_path.absolute()), "Path should be absolute"

    def test_get_top_info_new_format(self, temp_output_dir, new_style_top_content):
        """Test parsing new-style topology file with 5->3 marker."""
        top_file = temp_output_dir / "new_style.top"
        top_file.write_text(new_style_top_content)

        top_info = get_top_info(str(top_file))

        assert isinstance(top_info, TopInfo), "Should return TopInfo"
        assert top_info.nbases == 6, "Should have 6 bases"

    def test_get_top_info_malformed_raises(self, temp_output_dir):
        """Test that malformed topology raises RuntimeError."""
        bad_top = temp_output_dir / "bad.top"
        bad_top.write_text("this is not a valid topology header with too many fields\n")

        with pytest.raises(RuntimeError, match="Malformed topology header"):
            get_top_info(str(bad_top))


class TestGetTopInfoFromTraj:
    """Tests for get_top_info_from_traj() function."""

    def test_get_top_info_from_traj(self, mini_traj_path):
        """Test extracting topology info from trajectory file."""
        top_info = get_top_info_from_traj(str(mini_traj_path))

        assert isinstance(top_info, TopInfo), "Should return TopInfo"
        assert top_info.nbases == 132, "Should detect 132 bases"
        assert top_info.path == "", "Path should be empty when read from traj"


# =============================================================================
# Trajectory Info Tests
# =============================================================================

class TestGetTrajInfo:
    """Tests for get_traj_info() function."""

    def test_get_traj_info_creates_index(self, mini_traj_path, temp_output_dir, monkeypatch):
        """Test that get_traj_info creates .pyidx file."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj_no_idx.dat"
        copy(mini_traj_path, traj_copy)

        # Remove any existing index
        idx_path = Path(str(traj_copy) + ".pyidx")
        if idx_path.exists():
            idx_path.unlink()

        traj_info = get_traj_info(str(traj_copy))

        assert isinstance(traj_info, TrajInfo), "Should return TrajInfo"
        assert traj_info.nconfs > 0, "Should have configurations"
        assert idx_path.exists(), "Should create .pyidx file"

    def test_get_traj_info_uses_existing_index(self, mini_traj_path):
        """Test that get_traj_info uses existing .pyidx file."""
        traj_info = get_traj_info(str(mini_traj_path))

        assert isinstance(traj_info, TrajInfo), "Should return TrajInfo"
        assert len(traj_info.idxs) == traj_info.nconfs, "Index count should match nconfs"

    def test_get_traj_info_regenerates_stale_index(self, mini_traj_path, temp_output_dir):
        """Test that stale index file is regenerated."""
        traj_copy = temp_output_dir / "traj_stale.dat"
        copy(mini_traj_path, traj_copy)

        # Create a fake stale index with wrong size
        idx_path = Path(str(traj_copy) + ".pyidx")
        stale_idx = [ConfInfo(0, 100, 0)]  # Wrong size
        with open(idx_path, 'wb') as f:
            pickle.dump(stale_idx, f)

        traj_info = get_traj_info(str(traj_copy))

        assert traj_info.nconfs > 1, "Should regenerate index with correct count"

    def test_get_traj_info_detects_velocities(self, mini_traj_path):
        """Test velocity detection in trajectory."""
        traj_info = get_traj_info(str(mini_traj_path))

        # incl_v should be boolean
        assert isinstance(traj_info.incl_v, bool), "incl_v should be boolean"

    def test_get_traj_info_absolute_path(self, mini_traj_path):
        """Test that path is stored as absolute."""
        traj_info = get_traj_info(str(mini_traj_path))

        assert os.path.isabs(traj_info.path), "Path should be absolute"

    def test_get_traj_info_idxs_are_confinfo(self, mini_traj_path):
        """Test that index entries are ConfInfo objects."""
        traj_info = get_traj_info(str(mini_traj_path))

        for idx in traj_info.idxs:
            assert isinstance(idx, ConfInfo), "Each index should be ConfInfo"
            assert idx.offset >= 0, "Offset should be non-negative"
            assert idx.size > 0, "Size should be positive"


# =============================================================================
# Describe Tests
# =============================================================================

class TestDescribe:
    """Tests for the describe() function."""

    def test_describe_with_topology(self, old_top_path, mini_traj_path):
        """Test describe() with both topology and trajectory files."""
        top_info, traj_info = describe(str(old_top_path), str(mini_traj_path))

        assert isinstance(top_info, TopInfo), "Should return TopInfo"
        assert isinstance(traj_info, TrajInfo), "Should return TrajInfo"
        assert top_info.nbases == 132, "Should have correct base count"
        assert traj_info.nconfs > 0, "Should have configurations"

    def test_describe_without_topology(self, mini_traj_path):
        """Test describe() with None topology - reads from trajectory."""
        top_info, traj_info = describe(None, str(mini_traj_path))

        assert isinstance(top_info, TopInfo), "Should return TopInfo"
        assert top_info.nbases == 132, "Should infer base count from traj"
        assert top_info.path == "", "Path should be empty"

    def test_describe_returns_tuple(self, old_top_path, mini_traj_path):
        """Test that describe returns correct tuple structure."""
        result = describe(str(old_top_path), str(mini_traj_path))

        assert isinstance(result, tuple), "Should return tuple"
        assert len(result) == 2, "Should have two elements"


# =============================================================================
# Get Confs Tests
# =============================================================================

class TestGetConfs:
    """Tests for the get_confs() function."""

    def test_get_confs_returns_configurations(self, trajectory_info):
        """Test get_confs returns list of Configuration objects."""
        top_info, traj_info = trajectory_info

        confs = get_confs(top_info, traj_info, start_conf=0, n_confs=2)

        assert isinstance(confs, list), "Should return list"
        assert len(confs) == 2, "Should return 2 configurations"
        for conf in confs:
            assert isinstance(conf, Configuration), "Each should be Configuration"

    def test_get_confs_positions_shape(self, trajectory_info):
        """Test that positions have correct shape."""
        top_info, traj_info = trajectory_info

        confs = get_confs(top_info, traj_info, start_conf=0, n_confs=1)

        assert confs[0].positions.shape == (top_info.nbases, 3), \
            f"Positions should be ({top_info.nbases}, 3)"

    def test_get_confs_orientations_shape(self, trajectory_info):
        """Test that orientation vectors have correct shape."""
        top_info, traj_info = trajectory_info

        confs = get_confs(top_info, traj_info, start_conf=0, n_confs=1)

        assert confs[0].a1s.shape == (top_info.nbases, 3), "a1s should be (nbases, 3)"
        assert confs[0].a3s.shape == (top_info.nbases, 3), "a3s should be (nbases, 3)"

    def test_get_confs_box_shape(self, trajectory_info):
        """Test that box vector has correct shape."""
        top_info, traj_info = trajectory_info

        confs = get_confs(top_info, traj_info, start_conf=0, n_confs=1)

        assert confs[0].box.shape == (3,), "Box should be shape (3,)"
        assert np.all(confs[0].box > 0), "Box dimensions should be positive"

    def test_get_confs_different_start(self, trajectory_info):
        """Test get_confs with different start positions."""
        top_info, traj_info = trajectory_info

        conf_0 = get_confs(top_info, traj_info, start_conf=0, n_confs=1)[0]
        conf_1 = get_confs(top_info, traj_info, start_conf=1, n_confs=1)[0]

        # Different configurations should have different times or positions
        assert conf_0.time != conf_1.time or not np.allclose(conf_0.positions, conf_1.positions), \
            "Different start should give different configurations"


# =============================================================================
# Linear Read Tests
# =============================================================================

class TestLinearRead:
    """Tests for the linear_read() function."""

    def test_linear_read_yields_configurations(self, trajectory_info):
        """Test that linear_read yields configuration lists."""
        top_info, traj_info = trajectory_info

        chunk_count = 0
        for chunk in linear_read(traj_info, top_info, chunk_size=2):
            chunk_count += 1
            assert isinstance(chunk, list), "Should yield lists"
            for conf in chunk:
                assert isinstance(conf, Configuration), "List should contain Configurations"

        assert chunk_count > 0, "Should yield at least one chunk"

    def test_linear_read_covers_all_confs(self, trajectory_info):
        """Test that linear_read covers all configurations."""
        top_info, traj_info = trajectory_info

        total_confs = 0
        for chunk in linear_read(traj_info, top_info, chunk_size=3):
            total_confs += len(chunk)

        assert total_confs == traj_info.nconfs, "Should read all configurations"


# =============================================================================
# Strand Describe Tests
# =============================================================================

class TestStrandDescribe:
    """Tests for the strand_describe() function."""

    def test_strand_describe_old_format(self, old_top_path):
        """Test strand_describe with old-style topology."""
        system, monomers = strand_describe(str(old_top_path))

        assert isinstance(system, System), "Should return System"
        assert isinstance(monomers, list), "Should return monomer list"
        assert len(monomers) == 132, "Should have 132 monomers"

    def test_strand_describe_system_structure(self, old_top_path):
        """Test that System has correct strand structure."""
        system, monomers = strand_describe(str(old_top_path))

        assert len(system.strands) >= 1, "Should have at least one strand"
        for strand in system:
            assert isinstance(strand, Strand), "System should contain Strands"

    def test_strand_describe_monomer_types(self, old_top_path):
        """Test that monomers have correct base types."""
        system, monomers = strand_describe(str(old_top_path))

        valid_bases = {'A', 'C', 'G', 'T', 'U'}
        for m in monomers:
            assert isinstance(m, Monomer), "Should be Monomer objects"
            assert m.btype in valid_bases, f"Base type {m.btype} should be valid"

    def test_strand_describe_new_format(self, temp_output_dir, new_style_top_content):
        """Test strand_describe with new-style topology."""
        top_file = temp_output_dir / "new_top.top"
        top_file.write_text(new_style_top_content)

        system, monomers = strand_describe(str(top_file))

        assert len(monomers) == 6, "Should have 6 monomers"
        assert len(system.strands) == 2, "Should have 2 strands"

    def test_strand_describe_strand_types(self, temp_output_dir, new_style_top_content):
        """Test that strand types are correctly parsed."""
        top_file = temp_output_dir / "typed_top.top"
        top_file.write_text(new_style_top_content)

        system, monomers = strand_describe(str(top_file))

        assert system[0].type == "DNA", "First strand should be DNA"
        assert system[1].type == "RNA", "Second strand should be RNA"

    def test_strand_describe_neighbor_links(self, old_top_path):
        """Test that n3/n5 neighbors are set correctly."""
        system, monomers = strand_describe(str(old_top_path))

        # In old topology format (3'-5' direction), first monomer in strand
        # is at 3' end, so it has n3 = -1 or None (no 3' neighbor)
        first_monomer = monomers[0]
        assert first_monomer.n3 == -1 or first_monomer.n3 is None, \
            "First monomer n3 should be -1 or None (3' end of strand)"

    def test_strand_describe_circular_strand(self, temp_output_dir):
        """Test parsing circular strand."""
        circular_top = temp_output_dir / "circular.top"
        circular_top.write_text("5 1 5->3\nACGTA type=DNA circular=True\n")

        system, monomers = strand_describe(str(circular_top))

        assert system[0].circular is True, "Strand should be circular"
        # First monomer should link to last
        assert monomers[0].n5 == 4, "First should link to last in circular"
        assert monomers[4].n3 == 0, "Last should link to first in circular"


# =============================================================================
# Get Input Parameter Tests
# =============================================================================

class TestGetInputParameter:
    """Tests for the get_input_parameter() function."""

    def test_get_input_parameter_basic(self, input_file_path):
        """Test basic parameter extraction."""
        result = get_input_parameter(str(input_file_path), "backend")

        assert result == "CPU", "Should extract 'CPU' for backend"

    def test_get_input_parameter_with_spaces(self, input_file_path):
        """Test parameter extraction handles whitespace."""
        result = get_input_parameter(str(input_file_path), "topology")

        assert result == "rna_tile.top", "Should extract topology value"

    def test_get_input_parameter_numeric(self, input_file_path):
        """Test extraction of numeric parameters."""
        result = get_input_parameter(str(input_file_path), "steps")

        assert result == "1000", "Should extract steps as string"

    def test_get_input_parameter_not_found_raises(self, input_file_path):
        """Test that missing parameter raises RuntimeError."""
        with pytest.raises(RuntimeError, match="Key.*not found"):
            get_input_parameter(str(input_file_path), "nonexistent_parameter")

    def test_get_input_parameter_with_comment(self, temp_output_dir):
        """Test parameter extraction ignores inline comments."""
        input_file = temp_output_dir / "input_comment"
        input_file.write_text("param = value # this is a comment\n")

        result = get_input_parameter(str(input_file), "param")

        assert result == "value", "Should ignore comment"


# =============================================================================
# Inbox Tests
# =============================================================================

class TestInbox:
    """Tests for the inbox() function."""

    def test_inbox_basic(self, sample_configuration):
        """Test inbox returns Configuration with positions in box."""
        result = inbox(sample_configuration)

        assert isinstance(result, Configuration), "Should return Configuration"
        # All positions should be within [0, box)
        for i in range(3):
            assert np.all(result.positions[:, i] >= 0), f"Dim {i} should be >= 0"
            assert np.all(result.positions[:, i] < result.box[i]), f"Dim {i} should be < box"

    def test_inbox_preserves_metadata(self, sample_configuration):
        """Test that inbox preserves time, box, energy."""
        result = inbox(sample_configuration)

        assert result.time == sample_configuration.time, "Time should be preserved"
        np.testing.assert_array_equal(result.box, sample_configuration.box), "Box should be preserved"
        np.testing.assert_array_equal(result.energy, sample_configuration.energy), "Energy should be preserved"

    def test_inbox_preserves_orientations(self, sample_configuration):
        """Test that inbox preserves orientation vectors."""
        result = inbox(sample_configuration)

        np.testing.assert_array_equal(result.a1s, sample_configuration.a1s), "a1s should be preserved"
        np.testing.assert_array_equal(result.a3s, sample_configuration.a3s), "a3s should be preserved"

    def test_inbox_no_center(self, sample_configuration):
        """Test inbox without centering wraps positions into [0, box)."""
        result = inbox(sample_configuration, center=False)

        # Positions should be wrapped into [0, box) in each dimension
        for i in range(3):
            assert np.all(result.positions[:, i] >= 0), \
                f"Dimension {i}: all positions should be >= 0"
            assert np.all(result.positions[:, i] < result.box[i]), \
                f"Dimension {i}: all positions should be < box ({result.box[i]})"

    def test_inbox_bc_centerpoint(self, sample_configuration):
        """Test inbox with 'bc' centerpoint (box center).

        This should produce identical wrapping behavior to center=False,
        as 'bc' is the default argument and centers the structure at box center.
        """
        result = inbox(sample_configuration, center=True, centerpoint='bc')

        # Positions should be wrapped into [0, box) in each dimension
        for i in range(3):
            assert np.all(result.positions[:, i] >= 0), \
                f"Dimension {i}: all positions should be >= 0"
            assert np.all(result.positions[:, i] < result.box[i]), \
                f"Dimension {i}: all positions should be < box ({result.box[i]})"

    def test_inbox_custom_centerpoint(self):
        """Test inbox with custom center point places center-of-mass at target.

        When a custom centerpoint is provided, the resulting center-of-mass
        of the configuration should be at that centerpoint.
        """
        # Create a simple configuration with known positions
        conf = Configuration(
            time=0,
            box=np.array([10.0, 10.0, 10.0]),
            energy=np.array([0.0, 0.0, 0.0]),
            positions=np.array([
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [7.0, 8.0, 9.0]
            ]),
            a1s=np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            a3s=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
        )

        # Use a custom centerpoint that is NOT box center
        # Box center is [5, 5, 5], so use [2, 2, 2] to test off-center centering
        custom_center = np.array([2.0, 2.0, 2.0])
        result = inbox(conf, center=True, centerpoint=custom_center)

        # Positions should still be in box
        for i in range(3):
            assert np.all(result.positions[:, i] >= 0), \
                f"Dimension {i}: all positions should be >= 0"
            assert np.all(result.positions[:, i] < result.box[i]), \
                f"Dimension {i}: all positions should be < box ({result.box[i]})"

        # The center of mass should now be at the custom centerpoint
        resulting_com = np.mean(result.positions, axis=0)
        np.testing.assert_allclose(
            resulting_com, custom_center, atol=0.5,
            err_msg=f"Center of mass should be at custom centerpoint {custom_center}, "
                    f"but got {resulting_com}"
        )

    def test_inbox_wraps_outside_positions(self):
        """Test that positions outside box are wrapped."""
        conf = Configuration(
            time=0,
            box=np.array([10.0, 10.0, 10.0]),
            energy=np.array([0.0, 0.0, 0.0]),
            positions=np.array([[15.0, -3.0, 25.0]]),  # All outside box
            a1s=np.array([[1.0, 0.0, 0.0]]),
            a3s=np.array([[0.0, 0.0, 1.0]])
        )

        result = inbox(conf, center=False)

        # All should be wrapped into [0, 10)
        for i in range(3):
            assert 0 <= result.positions[0, i] < 10, f"Dim {i} should be wrapped"


# =============================================================================
# Write Conf Tests
# =============================================================================

class TestWriteConf:
    """Tests for write_conf() and conf_to_str() functions."""

    def test_write_conf_creates_file(self, temp_output_dir, sample_configuration):
        """Test write_conf creates valid file."""
        outpath = temp_output_dir / "test_conf.dat"
        write_conf(str(outpath), sample_configuration)

        assert outpath.exists(), "File should be created"

    def test_write_conf_content_structure(self, temp_output_dir, sample_configuration):
        """Test write_conf creates correctly formatted content."""
        outpath = temp_output_dir / "test_conf_structure.dat"
        write_conf(str(outpath), sample_configuration)

        content = outpath.read_text()
        lines = content.strip().split('\n')

        assert lines[0].startswith('t = '), "First line should be time"
        assert lines[1].startswith('b = '), "Second line should be box"
        assert lines[2].startswith('E = '), "Third line should be energy"
        assert len(lines) == 3 + len(sample_configuration.positions), "Should have header + positions"

    def test_write_conf_append_mode(self, temp_output_dir, sample_configuration):
        """Test write_conf append mode."""
        outpath = temp_output_dir / "test_append.dat"
        write_conf(str(outpath), sample_configuration, append=False)
        write_conf(str(outpath), sample_configuration, append=True)

        content = outpath.read_text()
        # Should have two configurations
        assert content.count('t = ') == 2, "Should have two time headers"

    def test_write_conf_no_velocities(self, temp_output_dir, sample_configuration):
        """Test write_conf without velocities."""
        outpath = temp_output_dir / "test_no_vel.dat"
        write_conf(str(outpath), sample_configuration, include_vel=False)

        content = outpath.read_text()
        lines = content.strip().split('\n')
        # Without velocities, each particle line should have 9 values (pos + a1 + a3)
        particle_line = lines[3].split()
        assert len(particle_line) == 9, "Should have 9 values without velocities"

    def test_write_conf_with_velocities(self, temp_output_dir, sample_configuration):
        """Test write_conf with velocities."""
        outpath = temp_output_dir / "test_with_vel.dat"
        write_conf(str(outpath), sample_configuration, include_vel=True)

        content = outpath.read_text()
        lines = content.strip().split('\n')
        # With velocities, each particle line should have 15 values
        particle_line = lines[3].split()
        assert len(particle_line) == 15, "Should have 15 values with velocities"

    def test_conf_to_str_returns_string(self, sample_configuration):
        """Test conf_to_str returns string."""
        result = conf_to_str(sample_configuration)

        assert isinstance(result, str), "Should return string"
        assert 't = ' in result, "Should contain time header"

    def test_conf_to_str_matches_write_conf(self, temp_output_dir, sample_configuration):
        """Test conf_to_str produces same content as write_conf."""
        outpath = temp_output_dir / "compare.dat"
        write_conf(str(outpath), sample_configuration, include_vel=True)
        file_content = outpath.read_text()

        str_content = conf_to_str(sample_configuration, include_vel=True)

        assert str_content == file_content, "String and file content should match"

    def test_write_conf_roundtrip(self, temp_output_dir, trajectory_info):
        """Test write then read produces equivalent configuration."""
        top_info, traj_info = trajectory_info

        original = get_confs(top_info, traj_info, 0, 1)[0]
        outpath = temp_output_dir / "roundtrip.dat"
        write_conf(str(outpath), original, include_vel=True)

        # Read back
        rt_top, rt_traj = describe(None, str(outpath))
        roundtrip = get_confs(rt_top, rt_traj, 0, 1)[0]

        np.testing.assert_allclose(roundtrip.positions, original.positions, rtol=1e-5,
                                   err_msg="Positions should match after roundtrip")
        np.testing.assert_allclose(roundtrip.a1s, original.a1s, rtol=1e-5,
                                   err_msg="a1s should match after roundtrip")
        np.testing.assert_allclose(roundtrip.a3s, original.a3s, rtol=1e-5,
                                   err_msg="a3s should match after roundtrip")


# =============================================================================
# Write Top Tests
# =============================================================================

class TestWriteTop:
    """Tests for write_top() and get_top_string() functions."""

    def test_get_top_string_new_format(self, temp_output_dir, new_style_top_content):
        """Test get_top_string produces valid new-style topology."""
        top_file = temp_output_dir / "input.top"
        top_file.write_text(new_style_top_content)

        system, _ = strand_describe(str(top_file))
        result = get_top_string(system, old_format=False)

        assert '5->3' in result, "Should have 5->3 marker"
        assert 'type=DNA' in result, "Should have type info"

    def test_get_top_string_old_format(self, old_top_path):
        """Test get_top_string produces valid old-style topology."""
        system, _ = strand_describe(str(old_top_path))
        result = get_top_string(system, old_format=True)

        assert '5->3' not in result, "Should not have 5->3 marker"
        lines = result.strip().split('\n')
        assert len(lines) > 1, "Should have multiple lines"

    def test_write_top_creates_file(self, temp_output_dir, new_style_top_content):
        """Test write_top creates file."""
        top_file = temp_output_dir / "input.top"
        top_file.write_text(new_style_top_content)

        system, _ = strand_describe(str(top_file))

        outpath = temp_output_dir / "output.top"
        write_top(str(outpath), system, old_format=False)

        assert outpath.exists(), "Output file should be created"

    def test_write_top_roundtrip_new_format(self, temp_output_dir, new_style_top_content):
        """Test write then read produces equivalent system (new format)."""
        top_file = temp_output_dir / "input.top"
        top_file.write_text(new_style_top_content)

        system, monomers = strand_describe(str(top_file))

        outpath = temp_output_dir / "output.top"
        write_top(str(outpath), system, old_format=False)

        # Read back
        system2, monomers2 = strand_describe(str(outpath))

        assert len(system2.strands) == len(system.strands), "Should have same strand count"
        assert len(monomers2) == len(monomers), "Should have same monomer count"

    def test_write_top_preserves_sequences(self, temp_output_dir, new_style_top_content):
        """Test that topology write preserves sequences."""
        top_file = temp_output_dir / "seq_test.top"
        top_file.write_text(new_style_top_content)

        system, _ = strand_describe(str(top_file))
        original_seqs = [s.get_sequence() for s in system]

        outpath = temp_output_dir / "seq_out.top"
        write_top(str(outpath), system, old_format=False)

        system2, _ = strand_describe(str(outpath))
        roundtrip_seqs = [s.get_sequence() for s in system2]

        assert original_seqs == roundtrip_seqs, "Sequences should be preserved"


# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_empty_trajectory_raises(self, temp_output_dir):
        """Test that empty file raises error."""
        empty_traj = temp_output_dir / "empty.dat"
        empty_traj.write_text("")

        with pytest.raises(Exception):
            get_traj_info(str(empty_traj))

    def test_invalid_trajectory_raises(self, temp_output_dir):
        """Test that invalid trajectory format raises error."""
        bad_traj = temp_output_dir / "bad_traj.dat"
        bad_traj.write_text("this is not a valid trajectory\n")

        with pytest.raises(Exception):
            # Should fail when trying to parse particle lines
            top_info, traj_info = describe(None, str(bad_traj))
            get_confs(top_info, traj_info, 0, 1)

    def test_mismatched_base_count(self, temp_output_dir, trajectory_info):
        """Test handling of mismatched topology and trajectory.

        Note: The Cython reader doesn't validate base count match at runtime.
        This test verifies the topology parsing works independently.
        """
        top_info, traj_info = trajectory_info

        # Create topology with different base count
        wrong_top = temp_output_dir / "wrong.top"
        wrong_top.write_text("10 1\n" + "1 A -1 1\n" * 10)

        wrong_top_info = get_top_info(str(wrong_top))

        # Topology parsing should work independently
        assert wrong_top_info.nbases == 10, "Should parse declared base count"
        assert traj_info.nconfs > 0, "Original traj should have confs"

    def test_box_size_types(self, trajectory_info):
        """Test that box dimensions are correct types."""
        top_info, traj_info = trajectory_info
        conf = get_confs(top_info, traj_info, 0, 1)[0]

        assert conf.box.dtype in [np.float64, np.float32], "Box should be float type"

    def test_time_is_integer(self, trajectory_info):
        """Test that time is integer type."""
        top_info, traj_info = trajectory_info
        conf = get_confs(top_info, traj_info, 0, 1)[0]

        assert isinstance(conf.time, (int, np.integer)), "Time should be integer"


# =============================================================================
# Integration Tests
# =============================================================================

class TestIntegration:
    """Integration tests combining multiple RyeReader functions."""

    def test_full_read_write_cycle(self, temp_output_dir, trajectory_info):
        """Test complete read-modify-write-read cycle."""
        top_info, traj_info = trajectory_info

        # Read
        original_confs = [get_confs(top_info, traj_info, i, 1)[0]
                         for i in range(min(3, traj_info.nconfs))]

        # Write
        outpath = temp_output_dir / "full_cycle.dat"
        for i, conf in enumerate(original_confs):
            write_conf(str(outpath), conf, append=(i > 0))

        # Read back
        new_top_info, new_traj_info = describe(None, str(outpath))

        assert new_traj_info.nconfs == len(original_confs), "Should have same conf count"
        assert new_top_info.nbases == top_info.nbases, "Should have same base count"

    def test_strand_operations_on_parsed_system(self, old_top_path):
        """Test Strand object methods after parsing."""
        system, monomers = strand_describe(str(old_top_path))

        strand = system[0]

        # Test various Strand methods
        assert strand.get_length() == len(strand.monomers), "get_length should match"
        assert isinstance(strand.get_sequence(), str), "get_sequence should return string"
        assert isinstance(strand.is_circular(), bool), "is_circular should return bool"
        assert isinstance(strand.is_nucleic(), bool), "is_nucleic should return bool"
        assert strand.is_old() is True, "Old topology should be marked as old"

    def test_system_iteration(self, old_top_path):
        """Test System iteration and indexing."""
        system, _ = strand_describe(str(old_top_path))

        # Test iteration
        strand_count = 0
        for strand in system:
            strand_count += 1
            assert isinstance(strand, Strand)

        assert strand_count == len(system), "Iteration should cover all strands"

        # Test indexing
        first_strand = system[0]
        assert isinstance(first_strand, Strand), "Indexing should return Strand"
