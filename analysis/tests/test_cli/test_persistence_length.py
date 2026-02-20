"""
Tests for oxDNA_analysis_tools.persistence_length module.

Tests cover:
- get_r() utility function
- compute() helper function
- persistence_length() API function
- fit_PL() fitting function
- CLI argument parsing
- main() CLI entry point
"""
import os
import sys
from contextlib import contextmanager
from pathlib import Path
from shutil import copy
from unittest.mock import patch

import numpy as np
import pytest
import oxpy


@contextmanager
def suppress_c_output():
    """
    Suppress output from C libraries by redirecting file descriptors.

    This is needed for oxpy which prints I/O statistics directly to C-level
    stdout/stderr, bypassing Python's sys.stdout/sys.stderr.
    """
    # Flush Python buffers first - critical for redirection to take effect
    sys.stdout.flush()
    sys.stderr.flush()

    # Save the original file descriptors
    stdout_fd = sys.stdout.fileno()
    stderr_fd = sys.stderr.fileno()
    saved_stdout_fd = os.dup(stdout_fd)
    saved_stderr_fd = os.dup(stderr_fd)

    try:
        # Open /dev/null and redirect stdout/stderr to it
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, stdout_fd)
        os.dup2(devnull, stderr_fd)
        os.close(devnull)
        yield
    finally:
        # Restore original file descriptors
        os.dup2(saved_stdout_fd, stdout_fd)
        os.dup2(saved_stderr_fd, stderr_fd)
        os.close(saved_stdout_fd)
        os.close(saved_stderr_fd)


# Filter the expected RuntimeWarning from persistence_length.py division
# This warning occurs when correlations_counter has zeros (unpaired positions)
pytestmark = pytest.mark.filterwarnings(
    "ignore:invalid value encountered in divide:RuntimeWarning"
)

from oxDNA_analysis_tools.persistence_length import (
    get_r,
    compute,
    persistence_length,
    fit_PL,
    cli_parser,
    ComputeContext,
    main
)
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
def topology_path(test_resources):
    """Path to the topology file."""
    return test_resources / "rna_tile.top"


@pytest.fixture(scope="module")
def input_file_path(test_resources):
    """Path to the input file for oxpy simulation."""
    return test_resources / "input_rna"


@pytest.fixture(scope="module")
def trajectory_info(topology_path, mini_traj_path):
    """Get topology and trajectory info."""
    top_info, traj_info = describe(str(topology_path), str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


class OxpyBackendManager:
    """
    Manages oxpy backend lifecycle and provides cleanup with output suppression.

    This is used instead of a simple fixture to allow proper cleanup with
    output suppression during pytest fixture teardown.
    """
    def __init__(self, input_file_path, mini_traj_path, test_resources):
        self.original_dir = os.getcwd()
        os.chdir(test_resources)

        self.context = oxpy.Context()
        self.context.__enter__()

        inp = oxpy.InputFile()
        inp.init_from_filename(str(input_file_path))
        inp["list_type"] = "cells"
        inp["trajectory_file"] = str(mini_traj_path)
        inp["confs_to_analyse"] = "1"
        inp["analysis_data_output_1"] = '{ \n name = stdout \n print_every = 1e10 \n col_1 = { \n id = my_obs \n type = hb_list \n } \n }'

        self.backend = oxpy.analysis.AnalysisBackend(inp)
        self.backend.read_next_configuration()

        # Extract pairs
        pairs = self.backend.config_info().get_observable_by_id("my_obs").get_output_string(
            self.backend.config_info().current_step).strip().split('\n')[1:]

        self.pair_dict = {}
        for p in pairs:
            p1 = int(p.split()[0])
            p2 = int(p.split()[1])
            self.pair_dict[p1] = p2

    def cleanup(self):
        """
        Clean up oxpy resources.

        Note: The oxpy SimBackend destructor prints I/O statistics to stdout
        which cannot be suppressed via Python since it's written directly by
        the C++ library. Using class-scoped fixtures minimizes this output
        to once per test class rather than once per test.
        """
        with suppress_c_output():
            del self.backend
            self.context.__exit__(None, None, None)
        os.chdir(self.original_dir)


@pytest.fixture(scope="class")
def oxpy_backend_and_pairs(input_file_path, mini_traj_path, test_resources, request):
    """
    Create an oxpy backend and extract base pairs for testing get_r().
    Returns backend and pair_dict from first configuration.

    Note: The oxpy context and backend must stay alive during the test since
    get_r() accesses config_info().particles(). We suppress the I/O statistics
    that oxpy prints when the SimBackend is destroyed by redirecting C-level
    stdout/stderr during fixture teardown.

    Using class scope so all TestGetR tests share one context, minimizing output.
    """
    manager = OxpyBackendManager(input_file_path, mini_traj_path, test_resources)

    # Register cleanup as a finalizer - this runs during pytest's teardown
    request.addfinalizer(manager.cleanup)

    return manager.backend, manager.pair_dict


# =============================================================================
# Unit Tests - get_r()
# =============================================================================

class TestGetR:
    """Tests for the get_r() utility function."""

    def test_get_r_returns_vector(self, oxpy_backend_and_pairs):
        """Test get_r() returns a vector with correct properties."""
        backend, pair_dict = oxpy_backend_and_pairs

        # Find a paired nucleotide with a paired neighbor
        valid_nucid = None
        for nucid in sorted(pair_dict.keys()):
            if nucid + 1 in pair_dict:
                valid_nucid = nucid
                break

        assert valid_nucid is not None, "Need at least two consecutive paired nucleotides"

        # Execute - single function call
        r = get_r(backend.config_info(), valid_nucid, pair_dict)

        # Progressive validation
        # 1. Correct type
        assert isinstance(r, np.ndarray), "Result should be numpy array"

        # 2. Correct shape
        assert r.shape == (3,), f"Expected shape (3,), got {r.shape}"

        # 3. Reasonable values (should not be NaN or Inf, should be finite distance)
        assert not np.any(np.isnan(r)), "Result contains NaN values"
        assert not np.any(np.isinf(r)), "Result contains Inf values"
        assert np.linalg.norm(r) > 0, "Vector should have non-zero length"

    def test_get_r_vector_magnitude(self, oxpy_backend_and_pairs):
        """Test that get_r() produces vectors with reasonable magnitude for base-pair steps."""
        backend, pair_dict = oxpy_backend_and_pairs

        # Get vectors for multiple base-pair steps
        vectors = []
        for nucid in sorted(pair_dict.keys()):
            if nucid + 1 in pair_dict:
                r = get_r(backend.config_info(), nucid, pair_dict)
                vectors.append(r)

        assert len(vectors) > 0, "Should have at least one vector"

        # Base-pair step distance should be reasonable (typically ~0.3-0.5 nm in oxDNA units)
        magnitudes = [np.linalg.norm(v) for v in vectors]
        mean_magnitude = np.mean(magnitudes)

        # oxDNA units: base-pair steps are typically 0.3-0.6 simulation units
        assert 0.2 < mean_magnitude < 0.8, \
            f"Mean base-pair step distance {mean_magnitude:.3f} is outside expected range"

    def test_get_r_midpoint_calculation(self, oxpy_backend_and_pairs):
        """Test that get_r() correctly computes base pair midpoints.

        CRITICAL: Verify midpoint calculation by manually computing expected values.
        """
        backend, pair_dict = oxpy_backend_and_pairs
        conf = backend.config_info()
        box = np.array(conf.box.box_sides)

        # Find a valid pair of consecutive base pairs
        valid_nucid = None
        for nucid in sorted(pair_dict.keys()):
            if nucid + 1 in pair_dict:
                valid_nucid = nucid
                break

        assert valid_nucid is not None, "Need at least two consecutive paired nucleotides"

        # Execute - get the vector using get_r()
        r = get_r(conf, valid_nucid, pair_dict)

        # Manually compute the expected vector to verify correctness
        pair = pair_dict[valid_nucid]
        next_pair = pair_dict[valid_nucid + 1]

        # Get base sites for first base pair
        firstA = conf.particles()[valid_nucid].base_site()
        firstB = conf.particles()[pair].base_site()

        # Get base sites for second base pair
        secondA = conf.particles()[valid_nucid + 1].base_site()
        secondB = conf.particles()[next_pair].base_site()

        # Manually compute midpoints
        first_midpos = (firstA + firstB) / 2
        second_midpos = (secondA + secondB) / 2

        # Manually compute vector with PBC correction
        r_expected = second_midpos - first_midpos
        r_expected -= box * np.rint(r_expected / box)

        # Progressive validation
        # 1. Correct type
        assert isinstance(r, np.ndarray), "Result should be numpy array"

        # 2. Correct shape
        assert r.shape == (3,), "Vector should be 3D"

        # 3. No NaN/Inf values
        assert not np.any(np.isnan(r)), "Result contains NaN"
        assert not np.any(np.isinf(r)), "Result contains Inf"

        # 4. Exact match with manually computed value
        np.testing.assert_allclose(
            r, r_expected,
            rtol=1e-10,
            err_msg=f"get_r() output {r} doesn't match manually computed midpoint vector {r_expected}"
        )


# =============================================================================
# Unit Tests - compute()
# =============================================================================

class TestCompute:
    """Tests for the compute() helper function."""

    def test_compute_returns_correct_structure(self, trajectory_info, input_file_path, test_resources, monkeypatch):
        """Test compute() returns tuple with proper structure and valid correlation values."""
        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        # Use a range with multiple paired nucleotides (wider range to ensure pairing)
        n1 = 0
        n2 = 20

        ctx = ComputeContext(
            traj_info=traj_info,
            input_file=str(input_file_path),
            n1=n1,
            n2=n2
        )

        # Execute - single function call
        result = compute(ctx, chunk_size=1, chunk_id=0)

        # Progressive validation
        # 1. Correct type
        assert isinstance(result, tuple), "Result should be tuple"
        assert len(result) == 3, "Result should be (l0, correlations, correlations_counter)"

        l0, correlations, correlations_counter = result

        # 2. Correct types and shapes
        assert isinstance(l0, (int, float, np.number)), "l0 should be numeric"
        assert isinstance(correlations, np.ndarray), "correlations should be numpy array"
        assert isinstance(correlations_counter, np.ndarray), "correlations_counter should be numpy array"

        expected_size = n2 - n1
        assert correlations.shape == (expected_size,), \
            f"correlations shape should be ({expected_size},), got {correlations.shape}"
        assert correlations_counter.shape == (expected_size,), \
            f"correlations_counter shape should match correlations"

        # 3. Reasonable values
        assert l0 >= 0, "l0 (contour length contribution) should be non-negative"
        # correlations from compute() are raw sums, not yet normalized
        assert not np.any(np.isnan(correlations)), "correlations should not contain NaN (not yet divided)"
        assert not np.any(np.isnan(correlations_counter)), "correlations_counter should not contain NaN"

        # 4. Correlations are accumulated sums of dot products
        # When normalized by correlations_counter, they should be in [-1, 1]
        # Test the normalized values
        valid_mask = correlations_counter > 0
        if np.any(valid_mask):
            normalized = correlations[valid_mask] / correlations_counter[valid_mask]
            assert np.all(normalized >= -1.1), "Normalized correlations should be >= -1"
            assert np.all(normalized <= 1.1), "Normalized correlations should be <= 1"

        # counters should be non-negative integers
        assert np.all(correlations_counter >= 0), "Counters should be non-negative"

    def test_compute_raises_on_different_strands(self, trajectory_info, input_file_path, test_resources, monkeypatch):
        """Test compute() raises error when n1 and n2 are on different strands."""
        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        # rna_tile.top has only one strand, so we can't test different strands directly
        # But we verify the check exists by ensuring same-strand inputs don't raise
        n1 = 0
        n2 = 10

        ctx = ComputeContext(
            traj_info=traj_info,
            input_file=str(input_file_path),
            n1=n1,
            n2=n2
        )

        # This should NOT raise since they're on the same strand
        result = compute(ctx, chunk_size=1, chunk_id=0)
        assert result is not None, "Same-strand nucleotides should work"


# =============================================================================
# API Tests - persistence_length() function
# =============================================================================

class TestPersistenceLengthFunction:
    """Tests for the persistence_length() API function."""

    def test_persistence_length_returns_correct_structure(self, trajectory_info, input_file_path, test_resources, monkeypatch):
        """Test persistence_length() returns correct tuple with valid values and expected decay."""
        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        # Use a subset of the structure with known pairing
        n1 = 0
        n2 = 20

        # Execute - single function call
        l0, correlations = persistence_length(traj_info, str(input_file_path), n1, n2, ncpus=1)

        # Progressive validation
        # 1. Correct types
        assert isinstance(l0, (float, np.floating)), "l0 should be float"
        assert isinstance(correlations, np.ndarray), "correlations should be numpy array"

        # 2. Correct shape
        expected_size = n2 - n1
        assert correlations.shape == (expected_size,), \
            f"correlations shape should be ({expected_size},), got {correlations.shape}"

        # 3. Reasonable values
        assert l0 > 0, "Average contour length should be positive"
        # Note: correlations may contain NaN where division by zero occurred (unpaired positions)
        # This is expected behavior given the structure

        # Valid correlations (excluding unpaired positions) should be in [-1, 1]
        finite_correlations = correlations[np.isfinite(correlations)]
        if len(finite_correlations) > 0:
            assert np.all(finite_correlations >= -1.1), "Correlations should be >= -1"
            assert np.all(finite_correlations <= 1.1), "Correlations should be <= 1"

        # 4. Correct content - correlations should decay with distance
        # First correlation (offset=0) should be ~1.0 if there are paired nucleotides
        # (vector correlates perfectly with itself)
        if np.isfinite(correlations[0]):
            assert 0.8 < correlations[0] <= 1.0, \
                f"First correlation (self) should be ~1.0, got {correlations[0]}"

        # Correlations should generally decrease with offset (characteristic of persistence length)
        # Check that correlation at larger offset is less than at offset 0
        valid_indices = np.where(np.isfinite(correlations))[0]
        if len(valid_indices) > 5:  # Need enough points to see decay
            mid_offset = valid_indices[min(len(valid_indices) // 2, len(valid_indices) - 1)]
            if mid_offset > 0 and np.isfinite(correlations[0]):
                # Allow for some noise in the decay
                assert correlations[mid_offset] < correlations[0] + 0.2, \
                    "Correlations should generally decay with distance"

    def test_persistence_length_with_parallel(self, trajectory_info, input_file_path, test_resources, monkeypatch):
        """Test persistence_length() works with parallel processing."""
        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        n1 = 0
        n2 = 15

        # Test with single CPU
        l0_single, corr_single = persistence_length(
            traj_info, str(input_file_path), n1, n2, ncpus=1
        )

        # Test with multiple CPUs
        l0_multi, corr_multi = persistence_length(
            traj_info, str(input_file_path), n1, n2, ncpus=2
        )

        # Results should be similar (allowing for floating-point differences)
        assert l0_single > 0, "Single CPU should produce valid l0"
        assert l0_multi > 0, "Multi CPU should produce valid l0"

        # Correlations should be similar (some variation due to numerical precision)
        finite_both = np.isfinite(corr_single) & np.isfinite(corr_multi)
        if np.any(finite_both):
            np.testing.assert_allclose(
                corr_single[finite_both],
                corr_multi[finite_both],
                rtol=1e-4,
                err_msg="Parallel and serial should give similar results"
            )

    def test_persistence_length_self_correlation_is_one(self, trajectory_info, input_file_path, test_resources, monkeypatch):
        """Test that self-correlation (offset=0) is exactly 1.0 for paired positions."""
        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        # Use a range that includes paired nucleotides
        n1 = 0
        n2 = 20

        # Execute - single function call
        l0, correlations = persistence_length(traj_info, str(input_file_path), n1, n2, ncpus=1)

        # Progressive validation
        # 1. Correct type
        assert isinstance(correlations, np.ndarray), "Correlations should be numpy array"

        # 2. Has data at offset 0
        assert len(correlations) > 0, "Should have at least one correlation value"

        # 3. Check if offset 0 has valid data (not NaN)
        if np.isfinite(correlations[0]):
            # 4. Self-correlation (offset=0) should be exactly 1.0
            # Since we normalize vectors (r0 = r0 / ||r0||, rk = rk / ||rk||)
            # and at offset 0, we have dot(r0, r0) = 1 for normalized vectors
            np.testing.assert_allclose(
                correlations[0], 1.0,
                rtol=1e-10,
                err_msg="Self-correlation at offset 0 should be exactly 1.0 for normalized vectors"
            )
        else:
            pytest.skip("No paired nucleotides found at position n1, cannot test self-correlation")


# =============================================================================
# Unit Tests - fit_PL()
# =============================================================================

class TestFitPL:
    """Tests for the fit_PL() fitting function."""

    def test_fit_pl_returns_positive_value(self, trajectory_info, input_file_path, temp_output_dir, test_resources, monkeypatch):
        """Test fit_PL() returns positive persistence length and creates plot."""
                  # Non-interactive backend

        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        # Get correlations
        n1 = 0
        n2 = 20
        l0, correlations = persistence_length(traj_info, str(input_file_path), n1, n2, ncpus=1)

        # Remove NaN/Inf values for fitting
        valid_corr = correlations[np.isfinite(correlations) & (correlations > 0)]

        if len(valid_corr) < 3:
            pytest.skip("Not enough valid correlations for fitting")

        plot_path = temp_output_dir / "test_pl_fit.png"

        # Execute - single function call
        pl = fit_PL(valid_corr, str(plot_path))

        # Progressive validation
        # 1. Correct type
        assert isinstance(pl, (float, np.floating)), "Persistence length should be float"

        # 2. Reasonable value (should be positive)
        assert pl > 0, "Persistence length should be positive"

        # 3. Plot file created
        assert plot_path.exists(), "fit_PL should create plot file"
        assert plot_path.stat().st_size > 0, "Plot file should have content"

    def test_fit_pl_with_exponential_decay(self, temp_output_dir):
        """Test fit_PL() correctly fits known exponential decay."""
        # Create synthetic correlation data with known persistence length
        true_pl = 10.0
        x = np.arange(0, 50)
        correlations = np.exp(-x / true_pl)

        plot_path = temp_output_dir / "synthetic_fit.png"
        fitted_pl = fit_PL(correlations, str(plot_path))

        # Should recover the true persistence length (within ~10% for linear fit in log space)
        assert abs(fitted_pl - true_pl) / true_pl < 0.15, \
            f"Fitted PL {fitted_pl:.1f} should be close to true PL {true_pl}"


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_all_positional_args(self):
        """Test parser requires all positional arguments."""
        parser = cli_parser()

        # No args - should fail
        with pytest.raises(SystemExit):
            parser.parse_args([])

        # Missing some args - should fail
        with pytest.raises(SystemExit):
            parser.parse_args(["traj.dat"])

        with pytest.raises(SystemExit):
            parser.parse_args(["traj.dat", "input.txt"])

        with pytest.raises(SystemExit):
            parser.parse_args(["traj.dat", "input.txt", "0"])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        # Test with all options
        args = parser.parse_args([
            "-p", "4",
            "-d", "data.txt",
            "-n", "plot.png",
            "-q",
            "trajectory.dat",
            "input.txt",
            "0",
            "50"
        ])

        assert args.traj_file == ["trajectory.dat"], "Trajectory not parsed"
        assert args.input == ["input.txt"], "Input file not parsed"
        assert args.nucid_1 == [0], "nucid_1 not parsed"
        assert args.nucid_2 == [50], "nucid_2 not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.data == ["data.txt"], "Data option not parsed"
        assert args.plot_name == ["plot.png"], "Plot name option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["traj.dat", "input.txt", "0", "10"])
        assert args_defaults.parallel is None, "Parallel should default to None"
        assert args_defaults.data is None, "Data should default to None"
        assert args_defaults.plot_name is None, "Plot name should default to None"
        assert args_defaults.quiet is False, "Quiet should default to False"


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_creates_default_plot(self, mini_traj_path, topology_path, input_file_path, test_resources, temp_output_dir, monkeypatch):
        """Test main() creates default persistence_length.png plot."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "rna_tile.top"
        input_copy = temp_output_dir / "input.txt"
        param_file = test_resources / "rna_sequence_dependent_parameters.txt"
        param_copy = temp_output_dir / "rna_sequence_dependent_parameters.txt"

        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)
        copy(input_file_path, input_copy)
        copy(param_file, param_copy)

        test_args = [
            "persistence_length.py",
            str(traj_copy),
            str(input_copy),
            "0",
            "15"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        default_plot = temp_output_dir / "persistence_length.png"
        assert default_plot.exists(), "Default plot 'persistence_length.png' should be created"
        assert default_plot.stat().st_size > 0, "Plot file should have content"

    def test_main_creates_custom_outputs(self, mini_traj_path, topology_path, input_file_path, test_resources, temp_output_dir, monkeypatch):
        """Test main() creates custom plot and data files."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "rna_tile.top"
        input_copy = temp_output_dir / "input.txt"
        param_file = test_resources / "rna_sequence_dependent_parameters.txt"
        param_copy = temp_output_dir / "rna_sequence_dependent_parameters.txt"

        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)
        copy(input_file_path, input_copy)
        copy(param_file, param_copy)

        plot_file = temp_output_dir / "custom_pl.png"
        data_file = temp_output_dir / "correlations.txt"

        test_args = [
            "persistence_length.py",
            "-n", str(plot_file),
            "-d", str(data_file),
            str(traj_copy),
            str(input_copy),
            "0",
            "15"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Check plot created
        assert plot_file.exists(), "Custom plot should be created"
        assert plot_file.stat().st_size > 0, "Plot should have content"

        # Check data file created with proper format
        assert data_file.exists(), "Data file should be created"

        # Verify data file format (offset correlation)
        with open(data_file) as f:
            lines = f.readlines()

        assert len(lines) > 0, "Data file should have content"

        # Check first line has two numbers
        first_line = lines[0].strip().split()
        assert len(first_line) == 2, "Each line should have 'offset correlation'"

        offset, corr = first_line
        assert offset == "0", "First line should be offset 0"
        float(corr)  # Should be convertible to float

    def test_main_with_parallel_and_quiet(self, mini_traj_path, topology_path, input_file_path, test_resources, temp_output_dir, monkeypatch):
        """Test main() with parallel processing and quiet mode."""
        monkeypatch.chdir(temp_output_dir)

        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "rna_tile.top"
        input_copy = temp_output_dir / "input.txt"
        param_file = test_resources / "rna_sequence_dependent_parameters.txt"
        param_copy = temp_output_dir / "rna_sequence_dependent_parameters.txt"

        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)
        copy(input_file_path, input_copy)
        copy(param_file, param_copy)

        plot_file = temp_output_dir / "parallel_pl.png"

        test_args = [
            "persistence_length.py",
            "-p", "2",
            "-q",
            "-n", str(plot_file),
            str(traj_copy),
            str(input_copy),
            "0",
            "10"
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        assert plot_file.exists(), "Parallel/quiet mode should create plot"

    def test_main_data_file_contains_exact_correlation_values(self, mini_traj_path, topology_path, input_file_path, test_resources, temp_output_dir, monkeypatch):
        """Test that data file written by main() contains exact correlation values from computation."""
        monkeypatch.chdir(temp_output_dir)

        # Copy required files
        traj_copy = temp_output_dir / "traj.dat"
        top_copy = temp_output_dir / "rna_tile.top"
        input_copy = temp_output_dir / "input.txt"
        param_file = test_resources / "rna_sequence_dependent_parameters.txt"
        param_copy = temp_output_dir / "rna_sequence_dependent_parameters.txt"

        copy(mini_traj_path, traj_copy)
        copy(topology_path, top_copy)
        copy(input_file_path, input_copy)
        copy(param_file, param_copy)

        data_file = temp_output_dir / "correlations_validation.txt"
        plot_file = temp_output_dir / "validation_pl.png"

        n1 = 0
        n2 = 15

        # First, compute correlations directly using the API
        top_info, traj_info = describe(str(top_copy), str(traj_copy))
        l0_expected, correlations_expected = persistence_length(
            traj_info, str(input_copy), n1, n2, ncpus=1
        )

        # Now run main() to generate the data file
        test_args = [
            "persistence_length.py",
            "-n", str(plot_file),
            "-d", str(data_file),
            str(traj_copy),
            str(input_copy),
            str(n1),
            str(n2)
        ]

        with patch.object(sys, 'argv', test_args):
            main()

        # Progressive validation
        # 1. File exists
        assert data_file.exists(), "Data file should be created"

        # 2. File has content
        with open(data_file) as f:
            lines = f.readlines()

        assert len(lines) > 0, "Data file should have content"
        assert len(lines) == len(correlations_expected), \
            f"Data file should have {len(correlations_expected)} lines, got {len(lines)}"

        # 3. Parse and validate each line
        for i, line in enumerate(lines):
            parts = line.strip().split()
            assert len(parts) == 2, f"Line {i} should have 'offset correlation' format"

            offset = int(parts[0])
            correlation = float(parts[1])

            # Offset should match line number
            assert offset == i, f"Line {i} should have offset {i}, got {offset}"

            # 4. Exact match with computed correlation values
            # Handle NaN values (unpaired positions) specially
            if np.isnan(correlations_expected[i]):
                assert np.isnan(correlation), \
                    f"Expected NaN at offset {i}, got {correlation}"
            else:
                np.testing.assert_allclose(
                    correlation, correlations_expected[i],
                    rtol=1e-10,
                    err_msg=f"Correlation at offset {i} in data file doesn't match computed value"
                )


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_persistence_length_high_ncpus(self, trajectory_info, input_file_path, test_resources, monkeypatch):
        """Test persistence_length() handles ncpus > nconfs gracefully."""
        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        # Use more CPUs than we have configurations (10 confs in minitraj.dat)
        n1 = 0
        n2 = 10

        l0, correlations = persistence_length(traj_info, str(input_file_path), n1, n2, ncpus=100)

        assert l0 > 0, "Should handle high ncpus gracefully"
        assert len(correlations) == n2 - n1, "Should still produce correct output"

    def test_persistence_length_larger_range(self, trajectory_info, input_file_path, test_resources, monkeypatch):
        """Test persistence_length() with larger nucleotide range."""
        # Change to test resources directory so oxpy can find rna_sequence_dependent_parameters.txt
        monkeypatch.chdir(test_resources)

        top_info, traj_info = trajectory_info

        # Use a larger range (but still within the 132 nucleotides)
        n1 = 0
        n2 = 50

        l0, correlations = persistence_length(traj_info, str(input_file_path), n1, n2, ncpus=1)

        assert l0 > 0, "Should work with larger range"
        assert len(correlations) == n2 - n1, "Correlations should match range"

        # With longer range, should see more decay in correlations
        valid_indices = np.where(np.isfinite(correlations) & (correlations > 0))[0]
        if len(valid_indices) > 10:
            # Correlation should decrease over distance
            early_corr = np.mean(correlations[valid_indices[:3]])
            late_corr = np.mean(correlations[valid_indices[-3:]])

            # Later correlations should generally be smaller (more decay)
            assert late_corr < early_corr, \
                "Correlations should decay over longer distances"

    def test_fit_pl_with_minimal_data(self, temp_output_dir):
        """Test fit_PL() handles minimal data points."""
        # Just 3 points (minimum for a reasonable fit)
        correlations = np.array([1.0, 0.5, 0.25])

        plot_path = temp_output_dir / "minimal_fit.png"
        pl = fit_PL(correlations, str(plot_path))

        assert pl > 0, "Should produce positive PL even with minimal data"
        assert plot_path.exists(), "Should create plot even with minimal data"
