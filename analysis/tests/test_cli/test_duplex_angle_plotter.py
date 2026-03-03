"""
Tests for oxDNA_analysis_tools.duplex_angle_plotter module.

Tests cover:
- rad2degree() function
- angle_between() function
- get_angle_between() function
- make_plots() function
- CLI argument parsing
- main() CLI entry point
"""
import sys
from pathlib import Path
from unittest.mock import patch
import json

import numpy as np
import pytest

from oxDNA_analysis_tools.duplex_angle_plotter import (
    rad2degree,
    angle_between,
    get_angle_between,
    make_plots,
    cli_parser,
    main
)


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


@pytest.fixture
def simple_angle_file(temp_output_dir):
    """
    Create a simple angle file with predictable data.

    Format: time, duplex_id, start1, end1, start2, end2, axisX, axisY, axisZ, helix_pos

    Note: Nucleotides are matched if they fall within EITHER strand range of a duplex.
    """
    angle_file = temp_output_dir / "test_angles.txt"

    # Header
    content = "time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n"

    # Timestep 0: Two duplexes with non-overlapping ranges
    # Duplex containing nucleotide 5 (in range 0-10) with axis [1, 0, 0]
    content += "0\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
    # Duplex containing nucleotide 55 (in range 50-60) with axis [0, 1, 0] - 90 degrees
    content += "0\t1\t50\t60\t150\t160\t0.0\t1.0\t0.0\t0.5\n"

    # Timestep 1: Same duplexes, different axes
    # Duplex with nucleotide 5: axis [1, 0, 0]
    content += "1\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
    # Duplex with nucleotide 55: axis [0.707, 0.707, 0] - 45 degrees
    content += "1\t1\t50\t60\t150\t160\t0.707\t0.707\t0.0\t0.5\n"

    # Timestep 2: One duplex missing (should produce NaN)
    # Only duplex with nucleotide 5
    content += "2\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"

    # Timestep 3: Both duplexes present again - needed to catch timestep 2
    content += "3\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
    content += "3\t1\t50\t60\t150\t160\t0.0\t0.0\t1.0\t0.5\n"

    with open(angle_file, 'w') as f:
        f.write(content)

    return angle_file


@pytest.fixture
def multi_angle_files(temp_output_dir):
    """Create multiple angle files for testing multiple trajectory handling."""
    files = []

    # File 1
    file1 = temp_output_dir / "angles1.txt"
    content1 = "time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n"
    content1 += "0\t0\t0\t5\t10\t15\t1.0\t0.0\t0.0\t0.5\n"
    content1 += "0\t1\t10\t15\t20\t25\t0.0\t1.0\t0.0\t0.5\n"
    with open(file1, 'w') as f:
        f.write(content1)
    files.append(file1)

    # File 2
    file2 = temp_output_dir / "angles2.txt"
    content2 = "time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n"
    content2 += "0\t0\t30\t35\t40\t45\t1.0\t0.0\t0.0\t0.5\n"
    content2 += "0\t1\t40\t45\t50\t55\t-1.0\t0.0\t0.0\t0.5\n"
    with open(file2, 'w') as f:
        f.write(content2)
    files.append(file2)

    return files


# =============================================================================
# Unit Tests - rad2degree()
# =============================================================================

class TestRad2Degree:
    """Tests for the rad2degree() function."""

    def test_rad2degree_known_values(self):
        """Test rad2degree with known conversion values."""
        # Execute - single call with multiple checks
        result_pi = rad2degree(np.pi)
        result_half_pi = rad2degree(np.pi / 2)
        result_zero = rad2degree(0)
        result_two_pi = rad2degree(2 * np.pi)

        # Progressive validation
        # 1. Correct type
        assert isinstance(result_pi, (float, np.floating)), "Result should be float"

        # 2. Values are reasonable (positive, finite)
        assert np.isfinite(result_pi), "Result should be finite"
        assert np.isfinite(result_half_pi), "Result should be finite"

        # 3. Correct content - known conversions
        assert np.isclose(result_pi, 180.0), "π radians should be 180 degrees"
        assert np.isclose(result_half_pi, 90.0), "π/2 radians should be 90 degrees"
        assert np.isclose(result_zero, 0.0), "0 radians should be 0 degrees"
        assert np.isclose(result_two_pi, 360.0), "2π radians should be 360 degrees"

    def test_rad2degree_array_input(self):
        """Test rad2degree with array input."""
        angles = np.array([0, np.pi/4, np.pi/2, np.pi])
        result = rad2degree(angles)

        # Check type and shape
        assert isinstance(result, np.ndarray), "Should work with arrays"
        assert result.shape == angles.shape, "Shape should be preserved"

        # Check values
        expected = np.array([0, 45, 90, 180])
        np.testing.assert_allclose(result, expected, rtol=1e-10)


# =============================================================================
# Unit Tests - angle_between()
# =============================================================================

class TestAngleBetween:
    """Tests for the angle_between() function."""

    def test_angle_between_orthogonal_vectors(self):
        """Test angle_between with orthogonal vectors."""
        # Setup - simple orthogonal vectors
        axis1 = np.array([1.0, 0.0, 0.0])
        axis2 = np.array([0.0, 1.0, 0.0])

        # Execute - single function call
        result = angle_between(axis1, axis2)

        # Progressive validation
        # 1. Correct type
        assert isinstance(result, (float, np.floating)), "Result should be float"

        # 2. Reasonable value (between 0 and π)
        assert 0 <= result <= np.pi, "Angle should be between 0 and π radians"
        assert np.isfinite(result), "Result should be finite"

        # 3. Correct content - orthogonal vectors should give π/2
        assert np.isclose(result, np.pi / 2, rtol=1e-10), "Orthogonal vectors should give π/2 radians"

    def test_angle_between_parallel_vectors(self):
        """Test angle_between with parallel vectors."""
        axis1 = np.array([1.0, 0.0, 0.0])
        axis2 = np.array([2.0, 0.0, 0.0])

        result = angle_between(axis1, axis2)

        # Type and bounds
        assert isinstance(result, (float, np.floating)), "Result should be float"
        assert 0 <= result <= np.pi, "Angle should be between 0 and π radians"

        # Parallel vectors should give 0
        assert np.isclose(result, 0.0, atol=1e-10), "Parallel vectors should give 0 radians"

    def test_angle_between_antiparallel_vectors(self):
        """Test angle_between with antiparallel vectors."""
        axis1 = np.array([1.0, 0.0, 0.0])
        axis2 = np.array([-1.0, 0.0, 0.0])

        result = angle_between(axis1, axis2)

        # Antiparallel vectors should give π
        assert np.isclose(result, np.pi, rtol=1e-10), "Antiparallel vectors should give π radians"

    def test_angle_between_45_degrees(self):
        """Test angle_between with 45-degree vectors."""
        axis1 = np.array([1.0, 0.0, 0.0])
        axis2 = np.array([1.0, 1.0, 0.0])  # Normalized will be [√2/2, √2/2, 0]

        result = angle_between(axis1, axis2)

        # 45 degrees = π/4 radians
        assert np.isclose(result, np.pi / 4, rtol=1e-10), "45-degree vectors should give π/4 radians"

    def test_angle_between_3d_vectors(self):
        """Test angle_between with general 3D vectors."""
        axis1 = np.array([1.0, 1.0, 1.0])
        axis2 = np.array([1.0, 0.0, 0.0])

        result = angle_between(axis1, axis2)

        # Type and bounds check
        assert isinstance(result, (float, np.floating)), "Result should be float"
        assert 0 <= result <= np.pi, "Angle should be between 0 and π radians"
        assert np.isfinite(result), "Result should be finite"

        # Calculate expected value
        dot_product = 1.0
        norm1 = np.sqrt(3.0)
        norm2 = 1.0
        expected = np.arccos(dot_product / (norm1 * norm2))
        assert np.isclose(result, expected, rtol=1e-10), "Should match calculated angle"

    def test_angle_between_nearly_parallel_vectors(self):
        """
        Tests that the function correctly handles floating-point precision issues
        when vectors are nearly parallel, which can cause dot_product/norms to be
        slightly greater than 1.0, leading to domain errors in arccos.
        """
        # Setup - vectors that are extremely close to parallel
        axis1 = np.array([1.0, 0.0, 0.0])
        axis2 = np.array([1.0, 1e-15, 0.0])  # Nearly parallel with tiny deviation

        # Execute - single function call
        result = angle_between(axis1, axis2)

        # Progressive validation
        # 1. Correct type
        assert isinstance(result, (float, np.floating)), "Result should be float"

        # 2. Result should be finite (not NaN due to arccos domain error)
        assert np.isfinite(result), "Result should be finite, not NaN from arccos domain error"

        # 3. Reasonable value (should be very close to 0 for nearly parallel vectors)
        assert 0 <= result <= np.pi, "Angle should be between 0 and π radians"
        assert result < 1e-10, "Nearly parallel vectors should give angle very close to 0"

        # Test with another precision edge case
        axis3 = np.array([1.0, 0.0, 0.0])
        axis4 = np.array([1.0, 1e-10, 1e-10])
        result2 = angle_between(axis3, axis4)

        assert np.isfinite(result2), "Should handle vectors with small perturbations"
        assert result2 < 1e-8, "Small perturbation should give small angle"


# =============================================================================
# API Tests - get_angle_between()
# =============================================================================

class TestGetAngleBetween:
    """Tests for the get_angle_between() function."""

    def test_get_angle_between_single_file(self, simple_angle_file):
        """Test get_angle_between with single file and particle pair."""
        # Setup
        files = [str(simple_angle_file)]
        p1s = [[5]]  # Search for nucleotide 5
        p2s = [[55]]  # Search for nucleotide 55
        invert_mask = [False]

        # Execute - single function call
        all_angles, means, medians, stdevs, representations = get_angle_between(
            files, p1s, p2s, invert_mask
        )

        # Progressive validation
        # 1. Correct types
        assert isinstance(all_angles, list), "all_angles should be list"
        assert isinstance(means, list), "means should be list"
        assert isinstance(medians, list), "medians should be list"
        assert isinstance(stdevs, list), "stdevs should be list"
        assert isinstance(representations, list), "representations should be list"

        # 2. Correct structure
        assert len(all_angles) == 1, "Should have one entry per file"
        assert len(all_angles[0]) == 1, "Should have one entry per particle pair"
        assert isinstance(all_angles[0][0], np.ndarray), "Angles should be numpy arrays"
        assert all_angles[0][0].shape[0] == 4, "Should have 4 timesteps"

        # 3. Reasonable values
        angles = all_angles[0][0]
        assert np.all((angles >= 0) | np.isnan(angles)), "Angles should be non-negative or NaN"
        assert np.all((angles <= 180) | np.isnan(angles)), "Angles should be <= 180 degrees or NaN"

        # 4. Correct content
        # Timestep 0: [1,0,0] vs [0,1,0] = 90 degrees
        assert np.isclose(angles[0], 90.0, rtol=1e-5), "First angle should be 90 degrees"
        # Timestep 1: [1,0,0] vs [0.707,0.707,0] = 45 degrees
        assert np.isclose(angles[1], 45.0, rtol=1e-2), "Second angle should be ~45 degrees"
        # Timestep 2: Missing second duplex, should be NaN
        assert np.isnan(angles[2]), "Third angle should be NaN (missing duplex)"
        # Timestep 3: [1,0,0] vs [0,0,1] = 90 degrees
        assert np.isclose(angles[3], 90.0, rtol=1e-5), "Fourth angle should be 90 degrees"

        # Check statistics
        assert len(means[0]) == 1, "Should have one mean per pair"
        assert len(medians[0]) == 1, "Should have one median per pair"
        assert len(stdevs[0]) == 1, "Should have one stdev per pair"
        assert len(representations[0]) == 1, "Should have one representation per pair"

        # Mean should be around (90 + 45 + 90) / 3 = 75 (ignoring NaN)
        assert np.isclose(means[0][0], 75.0, rtol=1e-5), "Mean should be between 70 and 80 degrees"
        # Representation should be 3/4 (3 valid out of 4 timesteps)
        assert np.isclose(representations[0][0], 3/4, rtol=1e-5), "Representation should be 3/4"

    def test_get_angle_between_with_invert(self, simple_angle_file):
        """Test get_angle_between with vector inversion."""
        files = [str(simple_angle_file)]
        p1s = [[5]]
        p2s = [[55]]
        invert_mask = [True]  # Invert first vector

        all_angles, _, _, _, _ = get_angle_between(files, p1s, p2s, invert_mask)

        # With inversion, [1,0,0] becomes [-1,0,0]
        # Timestep 0: [-1,0,0] vs [0,1,0] = still 90 degrees
        assert np.isclose(all_angles[0][0][0], 90.0, rtol=1e-5), "Inversion shouldn't affect orthogonal angle"

    def test_get_angle_between_multiple_pairs(self, temp_output_dir):
        """Test get_angle_between with multiple particle pairs in single file."""
        # Create file with data for multiple pairs
        angle_file = temp_output_dir / "multi_pair.txt"
        content = "time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n"
        # Timestep 0: Three duplexes
        content += "0\t0\t0\t5\t10\t15\t1.0\t0.0\t0.0\t0.5\n"  # Contains 3
        content += "0\t1\t10\t15\t20\t25\t0.0\t1.0\t0.0\t0.5\n"  # Contains 12
        content += "0\t2\t20\t25\t30\t35\t0.0\t0.0\t1.0\t0.5\n"  # Contains 22

        with open(angle_file, 'w') as f:
            f.write(content)

        files = [str(angle_file)]
        p1s = [[3, 12]]  # Two pairs
        p2s = [[12, 22]]
        invert_mask = [False, False]

        all_angles, means, medians, stdevs, representations = get_angle_between(
            files, p1s, p2s, invert_mask
        )

        # Structure checks
        assert len(all_angles[0]) == 2, "Should have two angle arrays"
        assert all_angles[0][0].shape[0] == 1, "Should have 1 timestep"
        assert all_angles[0][1].shape[0] == 1, "Should have 1 timestep"

        # Content checks
        # Pair 1: [1,0,0] vs [0,1,0] = 90 degrees
        assert np.isclose(all_angles[0][0][0], 90.0, rtol=1e-5)
        # Pair 2: [0,1,0] vs [0,0,1] = 90 degrees
        assert np.isclose(all_angles[0][1][0], 90.0, rtol=1e-5)

    def test_get_angle_between_multiple_files(self, multi_angle_files):
        """Test get_angle_between with multiple input files."""
        files = [str(f) for f in multi_angle_files]
        p1s = [[2], [32]]  # One pair per file
        p2s = [[12], [42]]
        invert_mask = [False, False]

        all_angles, means, medians, stdevs, representations = get_angle_between(
            files, p1s, p2s, invert_mask
        )

        # Structure
        assert len(all_angles) == 2, "Should have results for 2 files"
        assert len(means) == 2, "Should have means for 2 files"

        # Each file has one pair
        assert len(all_angles[0]) == 1, "File 1 should have 1 pair"
        assert len(all_angles[1]) == 1, "File 2 should have 1 pair"

    def test_get_angle_between_mismatched_inputs(self):
        """Test get_angle_between raises error with mismatched input lengths."""
        with pytest.raises(RuntimeError, match="Bad input arguments"):
            get_angle_between(
                files=["file1.txt", "file2.txt"],
                p1s=[[1]],  # Only one pair list
                p2s=[[2], [3]],
                invert_mask=[False]
            )

    def test_get_angle_between_mismatched_pairs(self, simple_angle_file):
        """Test get_angle_between raises error with mismatched pair lengths."""
        with pytest.raises(RuntimeError, match="Bad input arguments"):
            get_angle_between(
                files=[str(simple_angle_file)],
                p1s=[[1, 2]],  # Two p1 values
                p2s=[[3]],  # One p2 value
                invert_mask=[False]
            )

    def test_get_angle_between_exact_statistics(self, temp_output_dir):
        """
        Tests with known angle values (all exactly 90 degrees) to verify that
        mean is exactly 90, median is exactly 90, and stdev is exactly 0.
        """
        # Create angle file where all angles will be exactly 90 degrees
        angle_file = temp_output_dir / "exact_90.txt"
        content = "time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n"

        # Create multiple timesteps with orthogonal vectors (always 90 degrees)
        for t in range(10):
            # Duplex containing nucleotide 5: axis [1, 0, 0]
            content += f"{t}\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
            # Duplex containing nucleotide 55: axis [0, 1, 0] (orthogonal)
            content += f"{t}\t1\t50\t60\t150\t160\t0.0\t1.0\t0.0\t0.5\n"

        with open(angle_file, 'w') as f:
            f.write(content)

        files = [str(angle_file)]
        p1s = [[5]]
        p2s = [[55]]
        invert_mask = [False]

        # Execute - single function call
        all_angles, means, medians, stdevs, representations = get_angle_between(
            files, p1s, p2s, invert_mask
        )

        # Progressive validation
        # 1. Correct structure
        assert len(all_angles[0]) == 1, "Should have one angle array"
        assert all_angles[0][0].shape[0] == 10, "Should have 10 timesteps"

        # 2. All angles should be finite
        assert np.all(np.isfinite(all_angles[0][0])), "All angles should be finite"

        # 3. All angles should be exactly 90 degrees
        assert np.allclose(all_angles[0][0], 90.0, rtol=1e-10), "All angles should be exactly 90 degrees"

        # 4. Verify exact statistics
        # Mean should be exactly 90 degrees
        assert np.isclose(means[0][0], 90.0, rtol=1e-10), "Mean should be exactly 90 degrees"

        # Median should be exactly 90 degrees
        assert np.isclose(medians[0][0], 90.0, rtol=1e-10), "Median should be exactly 90 degrees"

        # Standard deviation should be exactly 0 (all values identical)
        assert np.isclose(stdevs[0][0], 0.0, atol=1e-10), "Standard deviation should be exactly 0"

        # Representation should be exactly 1.0 (all timesteps valid)
        assert representations[0][0] == 1.0, "Representation should be exactly 1.0"

        # Test alternate exact angle (all 0 degrees - parallel vectors)
        angle_file2 = temp_output_dir / "exact_0.txt"
        content2 = "time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n"
        for t in range(5):
            # Both duplexes with same axis [1, 0, 0]
            content2 += f"{t}\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
            content2 += f"{t}\t1\t50\t60\t150\t160\t1.0\t0.0\t0.0\t0.5\n"

        with open(angle_file2, 'w') as f:
            f.write(content2)

        all_angles2, means2, medians2, stdevs2, _ = get_angle_between(
            [str(angle_file2)], [[5]], [[55]], [False]
        )

        # All angles should be 0 degrees
        assert np.allclose(all_angles2[0][0], 0.0, atol=1e-10), "Parallel vectors give 0 degrees"
        assert np.isclose(means2[0][0], 0.0, atol=1e-10), "Mean should be exactly 0"
        assert np.isclose(stdevs2[0][0], 0.0, atol=1e-10), "Stdev should be exactly 0"

    def test_get_angle_between_malformed_input(self, temp_output_dir, capsys):
        """
        Tests that lines with non-numeric time fields are skipped with error messages,
        and processing continues with valid lines, producing correct results from good data.
        Note: The current implementation only handles time field parsing errors.
        Lines with missing fields after a valid time will cause IndexError.
        """
        # Create angle file with mix of good and malformed lines
        angle_file = temp_output_dir / "malformed.txt"
        content = "time\tduplex\tstart1\tend1\tstart2\tend2\taxisX\taxisY\taxisZ\thel_pos\n"

        # Timestep 0: Valid lines
        content += "0\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
        content += "0\t1\t50\t60\t150\t160\t0.0\t1.0\t0.0\t0.5\n"

        # Malformed line 1: Non-numeric time (this WILL be caught and skipped)
        content += "bad_time\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"

        # Timestep 1: Valid lines after malformed
        content += "1\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
        content += "1\t1\t50\t60\t150\t160\t0.707\t0.707\t0.0\t0.5\n"

        # Malformed line 2: Another non-numeric time
        content += "another_bad_time\t1\t5\t15\t105\t115\t0.0\t1.0\t0.0\t0.5\n"

        # Timestep 2: Valid lines
        content += "2\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"
        content += "2\t1\t50\t60\t150\t160\t0.0\t0.0\t1.0\t0.5\n"

        # Malformed line 3: Empty string as time
        content += "\t0\t0\t10\t100\t110\t1.0\t0.0\t0.0\t0.5\n"

        # Timestep 3: Valid lines (should still process)
        content += "3\t0\t0\t10\t100\t110\t0.0\t1.0\t0.0\t0.5\n"
        content += "3\t1\t50\t60\t150\t160\t1.0\t0.0\t0.0\t0.5\n"

        with open(angle_file, 'w') as f:
            f.write(content)

        files = [str(angle_file)]
        p1s = [[5]]
        p2s = [[55]]
        invert_mask = [False]

        # Execute - single function call
        all_angles, means, medians, stdevs, representations = get_angle_between(
            files, p1s, p2s, invert_mask
        )

        # Capture stderr/stdout to check for error messages
        captured = capsys.readouterr()

        # Progressive validation
        # 1. Should have successfully computed angles despite malformed lines
        assert len(all_angles[0]) == 1, "Should have one angle array"
        assert all_angles[0][0].shape[0] == 4, "Should have 4 timesteps (0, 1, 2, 3)"

        # 2. Check that error messages were printed for malformed lines
        assert "ERROR" in captured.out or "error" in captured.out.lower(), \
            "Should print error messages for malformed lines"
        assert "skip" in captured.out.lower(), "Should mention skipping bad lines"

        # 3. Valid timesteps should have correct angles
        angles = all_angles[0][0]

        # Timestep 0: [1,0,0] vs [0,1,0] = 90 degrees
        assert np.isclose(angles[0], 90.0, rtol=1e-5), "Timestep 0 should be 90 degrees"

        # Timestep 1: [1,0,0] vs [0.707,0.707,0] = 45 degrees
        assert np.isclose(angles[1], 45.0, rtol=1e-2), "Timestep 1 should be ~45 degrees"

        # Timestep 2: [1,0,0] vs [0,0,1] = 90 degrees
        assert np.isclose(angles[2], 90.0, rtol=1e-5), "Timestep 2 should be 90 degrees"

        # Timestep 3: [0,1,0] vs [1,0,0] = 90 degrees
        assert np.isclose(angles[3], 90.0, rtol=1e-5), "Timestep 3 should be 90 degrees"

        # 4. All angles should be finite (malformed lines didn't corrupt processing)
        assert np.all(np.isfinite(angles)), "All computed angles should be finite"

        # 5. Statistics should be reasonable
        assert 40 < means[0][0] < 95, "Mean should be reasonable despite malformed lines"
        assert representations[0][0] == 1.0, "All 4 timesteps should have valid angles"

        # 6. Verify multiple error outputs were generated (3 malformed lines)
        assert len(captured.out) > 100, "Should have printed error information for all bad lines"
        # Count the number of "ERROR" occurrences
        error_count = captured.out.count("ERROR")
        assert error_count >= 3, f"Should have at least 3 error messages, found {error_count}"


# =============================================================================
# Tests - make_plots()
# =============================================================================

class TestMakePlots:
    """Tests for the make_plots() function."""

    def test_make_plots_histogram(self, temp_output_dir):
        """Test make_plots creates histogram output."""
        # Setup - synthetic angle data
        all_angles = [[np.array([45.0, 90.0, 135.0])]]
        names = ["Test Angle"]
        outfile = str(temp_output_dir / "test_hist.png")

        # Execute
        result = make_plots(all_angles, names, outfile, hist=True, line=False)

        # Validation
        # 1. Function returns None (no return statement)
        assert result is None, "make_plots should return None"

        # 2. File is created
        assert Path(outfile).exists(), "Histogram file should be created"

        # 3. File has reasonable size (not empty)
        assert Path(outfile).stat().st_size > 0, "Output file should not be empty"

    def test_make_plots_trajectory(self, temp_output_dir):
        """Test make_plots creates trajectory line plot."""
        all_angles = [[np.array([45.0, 60.0, 75.0, 90.0])]]
        names = ["Trajectory"]
        outfile = str(temp_output_dir / "test_traj.png")

        result = make_plots(all_angles, names, outfile, hist=False, line=True)

        assert result is None, "make_plots should return None"
        assert Path(outfile).exists(), "Trajectory file should be created"
        assert Path(outfile).stat().st_size > 0, "Output file should not be empty"

    def test_make_plots_both(self, temp_output_dir):
        """Test make_plots creates both histogram and trajectory plots."""
        all_angles = [[np.array([30.0, 60.0, 90.0])]]
        names = ["Both Plots"]
        outfile = str(temp_output_dir / "test_both.png")

        result = make_plots(all_angles, names, outfile, hist=True, line=True)

        assert result is None, "make_plots should return None"

        # Should create both files
        hist_file = temp_output_dir / "test_both_hist.png"
        traj_file = temp_output_dir / "test_both_traj.png"

        assert hist_file.exists(), "Histogram file should be created"
        assert traj_file.exists(), "Trajectory file should be created"
        assert hist_file.stat().st_size > 0, "Histogram file should not be empty"
        assert traj_file.stat().st_size > 0, "Trajectory file should not be empty"

    def test_make_plots_multiple_series(self, temp_output_dir):
        """Test make_plots with multiple data series."""
        all_angles = [
            [np.array([45.0, 90.0])],
            [np.array([60.0, 120.0])]
        ]
        names = ["Series 1", "Series 2"]
        outfile = str(temp_output_dir / "multi_series.png")

        result = make_plots(all_angles, names, outfile, hist=True, line=False)

        assert result is None, "make_plots should return None"
        assert Path(outfile).exists(), "Output file should be created"

    def test_make_plots_with_nans(self, temp_output_dir):
        """Test make_plots handles NaN values correctly."""
        all_angles = [[np.array([45.0, np.nan, 90.0, np.nan, 135.0])]]
        names = ["With NaNs"]
        outfile = str(temp_output_dir / "with_nans.png")

        # Should not raise error
        result = make_plots(all_angles, names, outfile, hist=True, line=True)

        assert result is None, "make_plots should handle NaNs"
        assert Path(temp_output_dir / "with_nans_hist.png").exists(), "Should create histogram"
        assert Path(temp_output_dir / "with_nans_traj.png").exists(), "Should create trajectory"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_basic_arguments(self):
        """Test parser accepts basic required arguments."""
        parser = cli_parser()

        args = parser.parse_args([
            "-i", "angles.txt", "1", "2",
            "-o", "output.png"
        ])

        assert args.input == [["angles.txt", "1", "2"]], "Input should be parsed"
        assert args.output == "output.png", "Output should be parsed"

    def test_parser_multiple_inputs(self):
        """Test parser accepts multiple -i flags."""
        parser = cli_parser()

        args = parser.parse_args([
            "-i", "file1.txt", "1", "2",
            "-i", "file2.txt", "3", "4",
            "-o", "output.png"
        ])

        assert len(args.input) == 2, "Should have two input groups"
        assert args.input[0] == ["file1.txt", "1", "2"], "First input parsed"
        assert args.input[1] == ["file2.txt", "3", "4"], "Second input parsed"

    def test_parser_format_options(self):
        """Test parser accepts format options."""
        parser = cli_parser()

        args_hist = parser.parse_args([
            "-i", "file.txt", "1", "2",
            "-f", "histogram"
        ])
        assert args_hist.format == "histogram", "Histogram format should be parsed"

        args_traj = parser.parse_args([
            "-i", "file.txt", "1", "2",
            "-f", "trajectory"
        ])
        assert args_traj.format == "trajectory", "Trajectory format should be parsed"

        args_both = parser.parse_args([
            "-i", "file.txt", "1", "2",
            "-f", "both"
        ])
        assert args_both.format == "both", "Both format should be parsed"

    def test_parser_invert_mask(self):
        """Test parser accepts invert mask."""
        parser = cli_parser()

        args = parser.parse_args([
            "-i", "file.txt", "1", "2", "3", "4",
            "-v", "1", "0"
        ])

        assert args.invert_mask == ["1", "0"], "Invert mask should be parsed"

    def test_parser_data_output(self):
        """Test parser accepts data output option."""
        parser = cli_parser()

        args = parser.parse_args([
            "-i", "file.txt", "1", "2",
            "-d", "data.json"
        ])

        assert args.data == "data.json", "Data output should be parsed"

    def test_parser_names(self):
        """Test parser accepts custom names."""
        parser = cli_parser()

        args = parser.parse_args([
            "-i", "file.txt", "1", "2",
            "-n", "Name1", "Name2"
        ])

        assert args.names == ["Name1", "Name2"], "Names should be parsed"

    def test_parser_quiet_flag(self):
        """Test parser accepts quiet flag."""
        parser = cli_parser()

        args_quiet = parser.parse_args([
            "-i", "file.txt", "1", "2",
            "-q"
        ])
        assert args_quiet.quiet is True, "Quiet flag should be True"

        args_default = parser.parse_args([
            "-i", "file.txt", "1", "2"
        ])
        assert args_default.quiet is False, "Quiet should default to False"

    def test_parser_defaults(self):
        """Test parser has correct default values."""
        parser = cli_parser()

        args = parser.parse_args(["-i", "file.txt", "1", "2"])

        assert args.output is None, "Output should default to None"
        assert args.format is None, "Format should default to None"
        assert args.data is None, "Data should default to None"
        assert args.names is None, "Names should default to None"
        assert args.invert_mask is None, "Invert mask should default to None"
        assert args.quiet is False, "Quiet should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_histogram(self, simple_angle_file, temp_output_dir):
        """Test main() creates histogram output with default settings."""
        output_file = temp_output_dir / "main_output.png"

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-o", str(output_file),
            "-q"  # Quiet mode
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check output exists
        assert output_file.exists(), "Output file should be created"
        assert output_file.stat().st_size > 0, "Output file should not be empty"

    def test_main_creates_both_plots(self, simple_angle_file, temp_output_dir):
        """Test main() creates both histogram and trajectory plots."""
        output_file = temp_output_dir / "both_plots.png"

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-o", str(output_file),
            "-f", "both",
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        hist_file = temp_output_dir / "both_plots_hist.png"
        traj_file = temp_output_dir / "both_plots_traj.png"

        assert hist_file.exists(), "Histogram file should be created"
        assert traj_file.exists(), "Trajectory file should be created"

    def test_main_creates_trajectory_plot(self, simple_angle_file, temp_output_dir):
        """Test main() creates trajectory plot."""
        output_file = temp_output_dir / "trajectory.png"

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-o", str(output_file),
            "-f", "trajectory",
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Trajectory file should be created"

    def test_main_creates_data_json(self, simple_angle_file, temp_output_dir):
        """Test main() creates JSON data output."""
        output_file = temp_output_dir / "plot.png"
        data_file = temp_output_dir / "angles_data"

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-o", str(output_file),
            "-d", str(data_file),
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        json_file = temp_output_dir / "angles_data.json"
        assert json_file.exists(), "JSON data file should be created"

        # Validate JSON structure
        with open(json_file, 'r') as f:
            data = json.load(f)

        assert isinstance(data, dict), "JSON should be a dictionary"
        assert "5-55" in data, "Should contain angle pair key"
        assert isinstance(data["5-55"], list), "Values should be lists"

    def test_main_with_custom_names(self, simple_angle_file, temp_output_dir):
        """Test main() with custom names for data series."""
        output_file = temp_output_dir / "named_output.png"
        data_file = temp_output_dir / "named_data"

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-o", str(output_file),
            "-d", str(data_file),
            "-n", "Custom Name",
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check JSON uses custom name
        json_file = temp_output_dir / "named_data.json"
        with open(json_file, 'r') as f:
            data = json.load(f)

        # Note: JSON uses particle IDs as keys, not custom names
        # Custom names are for plot legends
        assert len(data) > 0, "JSON should contain data"

    def test_main_with_invert_mask(self, simple_angle_file, temp_output_dir):
        """Test main() with invert mask."""
        output_file = temp_output_dir / "inverted.png"

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-o", str(output_file),
            "-v", "1",
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output should be created with invert mask"

    def test_main_multiple_inputs(self, multi_angle_files, temp_output_dir):
        """Test main() with multiple input files."""
        output_file = temp_output_dir / "multi_input.png"

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(multi_angle_files[0]), "2", "12",
            "-i", str(multi_angle_files[1]), "32", "42",
            "-o", str(output_file),
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert output_file.exists(), "Output should be created for multiple inputs"

    def test_main_default_output_name(self, simple_angle_file, temp_output_dir, monkeypatch):
        """Test main() uses default output filename when not specified."""
        monkeypatch.chdir(temp_output_dir)

        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        default_file = temp_output_dir / "angle.png"
        assert default_file.exists(), "Should create default 'angle.png' file"
        # Clean up
        default_file.unlink()

    def test_main_invalid_format_raises_error(self, simple_angle_file, temp_output_dir):
        """Test main() raises error with invalid format."""
        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55",
            "-f", "invalid_format",
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            with pytest.raises(RuntimeError, match="Unrecognized graph format"):
                main()

    def test_main_mismatched_invert_mask_raises_error(self, simple_angle_file):
        """Test main() raises error when invert mask length doesn't match pairs."""
        test_args = [
            "duplex_angle_plotter.py",
            "-i", str(simple_angle_file), "5", "55", "6", "26",  # 2 pairs
            "-v", "1",  # Only 1 mask value
            "-q"
        ]

        with patch.object(sys, 'argv', test_args):
            with pytest.raises(RuntimeError, match="invert mask"):
                main()
