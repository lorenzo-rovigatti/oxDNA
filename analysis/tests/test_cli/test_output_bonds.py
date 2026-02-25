"""
Tests for oxDNA_analysis_tools.output_bonds module.

Tests cover:
- output_bonds() API function
- parse_header() helper function
- _resolve_fields() helper function
- CLI argument parsing
- main() CLI entry point

Note: These tests require oxpy and a valid input file.
"""
import json
import sys
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from oxDNA_analysis_tools.UTILS.RyeReader import describe

# Skip all tests if oxpy is not available
pytest.importorskip("oxpy")

from oxDNA_analysis_tools.output_bonds import (
    output_bonds,
    parse_header,
    _resolve_fields,
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


@pytest.fixture(scope="module")
def mini_traj_path(test_resources):
    """Path to the mini trajectory file."""
    return test_resources / "minitraj.dat"


@pytest.fixture(scope="module")
def input_file_path(test_resources):
    """Path to the input file."""
    return test_resources / "input_rna"


@pytest.fixture(scope="module")
def trajectory_info(mini_traj_path, test_resources):
    """Get topology and trajectory info for the mini trajectory."""
    top_file = test_resources / "rna_tile.top"
    top_info, traj_info = describe(str(top_file), str(mini_traj_path))
    return top_info, traj_info


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary directory for output files."""
    return tmp_path


# =============================================================================
# Helper Function Tests
# =============================================================================

class TestParseHeader:
    """Tests for the parse_header() helper function."""

    def test_parse_header_extracts_potentials(self):
        """Test parse_header extracts potential names from header string."""
        header = "#id1 id2 fene bexc stack nexc hb cr_stack cx_stack total, t = 1000000"
        result = parse_header(header)

        assert isinstance(result, list), "Should return list"
        assert "fene" in result, "Should contain 'fene'"
        assert "total" in result, "Should contain 'total'"
        # id1 and id2 should be excluded
        assert "id1" not in result, "Should not contain 'id1'"
        assert "id2" not in result, "Should not contain 'id2'"


class TestResolveFields:
    """Tests for the _resolve_fields() helper function."""

    def test_resolve_fields_none_returns_all(self):
        """None field_names returns all fields unchanged."""
        pot_names = ['fene', 'hb', 'total']
        indices, active = _resolve_fields(pot_names, None)
        assert indices is None
        assert active == pot_names

    def test_resolve_fields_all_keyword(self):
        """'all' keyword returns all fields."""
        pot_names = ['fene', 'hb', 'total']
        indices, active = _resolve_fields(pot_names, ['all'])
        assert indices is None
        assert active == pot_names

    def test_resolve_fields_subset(self):
        """Requesting a subset returns correct indices and names."""
        pot_names = ['fene', 'bexc', 'stack', 'nexc', 'hb', 'total']
        indices, active = _resolve_fields(pot_names, ['fene', 'hb'])
        assert indices == [0, 4]
        assert active == ['fene', 'hb']

    def test_resolve_fields_dh_alias(self):
        """'dh' resolves to Debye-Huckel field."""
        pot_names = ['fene', 'hb', 'Debye-Huckel', 'total']
        indices, active = _resolve_fields(pot_names, ['dh'])
        assert indices == [2]
        assert active == ['Debye-Huckel']

    def test_resolve_fields_unknown_field_warns(self):
        """Unknown field names are skipped; valid ones are kept."""
        pot_names = ['fene', 'hb', 'total']
        indices, active = _resolve_fields(pot_names, ['fene', 'nonexistent'])
        assert indices == [0]
        assert active == ['fene']

    def test_resolve_fields_no_valid_fields_returns_all(self):
        """If no valid fields given, fall back to all fields."""
        pot_names = ['fene', 'hb', 'total']
        indices, active = _resolve_fields(pot_names, ['nonexistent'])
        assert indices is None
        assert active == pot_names


# =============================================================================
# API Tests - output_bonds() function
# =============================================================================

class TestOutputBondsFunction:
    """Tests for the output_bonds() API function."""

    def test_output_bonds_visualize_returns_energy_array(
        self, trajectory_info, input_file_path, test_resources, monkeypatch
    ):
        """output_bonds(visualize=True) returns mean per-nucleotide energy array and field names."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info

        energies, pot_names = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=1, ncpus=1
        )

        assert isinstance(energies, np.ndarray), "Energies should be a numpy array"
        assert isinstance(pot_names, list), "Potential names should be a list"
        assert energies.shape == (top_info.nbases, len(pot_names)), \
            "Energies shape should be (nbases, n_potentials)"
        assert all(isinstance(n, str) for n in pot_names), "Potential names should be strings"
        # Returned array is already the mean; values should be finite
        assert np.all(np.isfinite(energies)), "Mean per-nucleotide energies should be finite"

    def test_output_bonds_data_file_output(
        self, trajectory_info, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """output_bonds(data_file=...) writes pair data incrementally and returns None energies."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info
        out = temp_output_dir / "bonds.txt"

        energies, pot_names = output_bonds(
            traj_info, top_info, str(input_file_path),
            data_file=str(out)
        )

        assert energies is None, "Energies should be None when visualize=False"
        assert out.exists(), "Data file should be created"
        content = out.read_text()
        assert content.startswith('# p1 p2'), "File should start with field header"
        # Should contain at least one pair line (two integer IDs + energies)
        pair_lines = [l for l in content.splitlines() if l and not l.startswith('#')]
        assert len(pair_lines) > 0, "Should contain pair data lines"
        parts = pair_lines[0].split()
        assert len(parts) == 2 + len(pot_names), \
            "Each pair line should have p, q, then one value per field"

    def test_output_bonds_field_filtering(
        self, trajectory_info, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """Field filtering reduces the number of energy columns in output."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info

        # All fields (parallel path via visualize=True)
        energies_all, pot_names_all = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=1, ncpus=1
        )

        # Filter to the first two fields only
        requested = pot_names_all[:2]
        energies_filtered, pot_names_filtered = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=1, ncpus=1,
            fields=requested
        )

        assert pot_names_filtered == requested, "Filtered pot_names should match request"
        assert energies_filtered.shape == (top_info.nbases, 2), \
            "Filtered energies should have 2 columns"

        # Verify values match the corresponding columns in the full array
        for fi, name in enumerate(requested):
            ai = pot_names_all.index(name)
            np.testing.assert_array_almost_equal(
                energies_filtered[:, fi], energies_all[:, ai],
                err_msg=f"Filtered column '{name}' should match full column"
            )

    def test_output_bonds_conversion_factor(
        self, trajectory_info, input_file_path, test_resources, monkeypatch
    ):
        """Conversion factor scales all energy values proportionally."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info

        e1, _ = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=1, ncpus=1
        )
        e2, _ = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, conversion_factor=2, ncpus=1
        )

        mean1 = np.mean(np.abs(e1))
        mean2 = np.mean(np.abs(e2))
        assert mean1 > 0, "Non-zero energies expected"
        ratio = mean2 / mean1
        assert 1.5 < ratio < 2.5, f"Conversion factor should scale energies, got ratio {ratio}"

    def test_output_bonds_visualize_plus_data_file(
        self, trajectory_info, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """visualize=True with data_file uses serial path: both energies and file are produced."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info
        out = temp_output_dir / "bonds_v.txt"

        energies, pot_names = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, data_file=str(out)
        )

        assert isinstance(energies, np.ndarray), "Should return energy array when visualize=True"
        assert energies.shape == (top_info.nbases, len(pot_names))
        assert out.exists(), "Data file should also be created"
        content = out.read_text()
        assert content.startswith('# p1 p2'), "Data file should have field header"

    def test_output_bonds_parallel_with_serial_mode_raises(
        self, trajectory_info, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """Combining -p with any serial-mode flag raises RuntimeError."""
        monkeypatch.chdir(test_resources)
        top_info, traj_info = trajectory_info

        # -t with -p should raise
        with pytest.raises(RuntimeError, match="incompatible with -p"):
            output_bonds(
                traj_info, top_info, str(input_file_path),
                traj_view=str(temp_output_dir / "traj.json"), ncpus=2
            )

        # -d with -p should raise
        with pytest.raises(RuntimeError, match="incompatible with -p"):
            output_bonds(
                traj_info, top_info, str(input_file_path),
                data_file=str(temp_output_dir / "bonds.txt"), ncpus=2
            )

        # -v with -p alone should NOT raise (parallel path)
        energies, _ = output_bonds(
            traj_info, top_info, str(input_file_path),
            visualize=True, ncpus=2
        )
        assert isinstance(energies, np.ndarray)


# =============================================================================
# CLI Parser Tests
# =============================================================================

class TestCLIParser:
    """Tests for the CLI argument parser."""

    def test_parser_requires_arguments(self):
        """Test that parser requires inputfile and trajectory arguments."""
        parser = cli_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_parser_all_options(self):
        """Test parser accepts all options and sets correct defaults."""
        parser = cli_parser()

        args = parser.parse_args([
            "-v", "output.json",
            "-p", "4",
            "-u", "pNnm",
            "-q",
            "input.txt",
            "trajectory.dat"
        ])

        assert args.inputfile == ["input.txt"], "Inputfile not parsed"
        assert args.trajectory == ["trajectory.dat"], "Trajectory not parsed"
        assert args.outfile == ["output.json"], "Output option not parsed"
        assert args.parallel == [4], "Parallel option not parsed"
        assert args.units == ["pNnm"], "Units option not parsed"
        assert args.quiet is True, "Quiet option not parsed"

        # Test defaults
        args_defaults = parser.parse_args(["input.txt", "traj.dat"])
        assert args_defaults.quiet is False, "Quiet should default to False"
        assert args_defaults.outfile is None, "Outfile should default to None"
        assert args_defaults.traj_view is None, "traj_view should default to None"
        assert args_defaults.fields is None, "fields should default to None"
        assert args_defaults.data_file is None, "data_file should default to None"
        assert args_defaults.force_print is False, "force_print should default to False"

    def test_parser_new_flags(self):
        """Test that new -t, -f, -d, --force_print flags are parsed correctly."""
        parser = cli_parser()
        args = parser.parse_args([
            "-t", "traj_out.json",
            "-f", "fene", "hb",
            "-d", "bonds.txt",
            "--force_print",
            "input.txt",
            "trajectory.dat"
        ])
        assert args.traj_view == ["traj_out.json"]
        assert args.fields == ["fene", "hb"]
        assert args.data_file == ["bonds.txt"]
        assert args.force_print is True


# =============================================================================
# CLI main() Tests
# =============================================================================

class TestMain:
    """Tests for the main() CLI entry point."""

    def test_main_visualize_creates_files(
        self, mini_traj_path, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """main() with -v creates oxView overlay JSON files with correct structure."""
        monkeypatch.chdir(test_resources)
        output_base = temp_output_dir / "energies.json"

        with patch.object(sys, 'argv', [
            "output_bonds.py",
            "-v", str(output_base),
            str(input_file_path),
            str(mini_traj_path)
        ]):
            main()

        json_files = list(temp_output_dir.glob("energies_*.json"))
        assert len(json_files) > 0, "Should create at least one JSON file per potential"

        # Verify JSON structure: {"field_name (units)": [per_nucleotide_values]}
        with open(json_files[0]) as f:
            data = json.load(f)
        assert len(data) == 1, "Each file should have exactly one top-level key"
        values = list(data.values())[0]
        assert isinstance(values, list), "Values should be a list"

    def test_main_traj_view_creates_files(
        self, mini_traj_path, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """main() with -t creates per-frame oxView trajectory JSON files."""
        monkeypatch.chdir(test_resources)
        output_base = temp_output_dir / "traj_overlay.json"
        _, traj_info = describe(None, str(mini_traj_path))

        with patch.object(sys, 'argv', [
            "output_bonds.py",
            "-t", str(output_base),
            str(input_file_path),
            str(mini_traj_path)
        ]):
            main()

        json_files = list(temp_output_dir.glob("traj_overlay_*.json"))
        assert len(json_files) > 0, "Should create at least one trajectory JSON file"

        # Verify structure: {"field (units)": [[frame0_values], [frame1_values], ...]}
        with open(json_files[0]) as f:
            data = json.load(f)
        assert len(data) == 1, "Each file should have exactly one top-level key"
        frames = list(data.values())[0]
        assert isinstance(frames, list), "Value should be a list of frames"
        assert len(frames) == traj_info.nconfs, "Should have one entry per configuration"
        assert isinstance(frames[0], list), "Each frame entry should be a list of per-nucleotide values"

    def test_main_data_file_flag(
        self, mini_traj_path, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """main() with -d writes bond data to the specified file."""
        monkeypatch.chdir(test_resources)
        data_out = temp_output_dir / "bonds.txt"

        with patch.object(sys, 'argv', [
            "output_bonds.py",
            "-d", str(data_out),
            str(input_file_path),
            str(mini_traj_path)
        ]):
            main()

        assert data_out.exists(), "Data file should be created"
        content = data_out.read_text()
        assert content.startswith('# p1 p2'), "Data file should start with field header"
        pair_lines = [l for l in content.splitlines() if l and not l.startswith('#')]
        assert len(pair_lines) > 0, "Data file should contain pair data"

    def test_main_fields_filter(
        self, mini_traj_path, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """main() with -f only writes the requested energy fields."""
        monkeypatch.chdir(test_resources)
        data_all = temp_output_dir / "bonds_all.txt"
        data_filtered = temp_output_dir / "bonds_filtered.txt"

        # Run with all fields
        with patch.object(sys, 'argv', [
            "output_bonds.py", "-d", str(data_all),
            str(input_file_path), str(mini_traj_path)
        ]):
            main()

        # Header format: "# p1 p2 FIELD1 FIELD2 ... (units)"
        # split() gives ['#', 'p1', 'p2', 'FIELD1', ...]; index 3 is the first field
        all_header = data_all.read_text().splitlines()[0]
        first_field = all_header.split()[3]

        # Run with only the first field
        with patch.object(sys, 'argv', [
            "output_bonds.py",
            "-f", first_field,
            "-d", str(data_filtered),
            str(input_file_path), str(mini_traj_path)
        ]):
            main()

        filtered_header = data_filtered.read_text().splitlines()[0]
        assert first_field in filtered_header, "Filtered header should mention the requested field"
        # Pair lines in the filtered file should have fewer columns (p q field) vs all-fields
        all_pair_lines = [l for l in data_all.read_text().splitlines() if l and not l.startswith('#')]
        filt_pair_lines = [l for l in data_filtered.read_text().splitlines() if l and not l.startswith('#')]
        assert len(all_pair_lines[0].split()) > len(filt_pair_lines[0].split()), \
            "Filtered output should have fewer columns than full output"

    def test_main_with_units_option(
        self, mini_traj_path, input_file_path, test_resources, temp_output_dir, monkeypatch
    ):
        """main() with different -u options produces correctly scaled energies."""
        monkeypatch.chdir(test_resources)

        output_pn = temp_output_dir / "energies_pn.json"
        with patch.object(sys, 'argv', [
            "output_bonds.py",
            "-v", str(output_pn),
            "-u", "pNnm",
            str(input_file_path),
            str(mini_traj_path)
        ]):
            main()

        output_ox = temp_output_dir / "energies_ox.json"
        with patch.object(sys, 'argv', [
            "output_bonds.py",
            "-v", str(output_ox),
            "-u", "oxDNA",
            str(input_file_path),
            str(mini_traj_path)
        ]):
            main()

        pn_files = list(temp_output_dir.glob("energies_pn_*.json"))
        ox_files = list(temp_output_dir.glob("energies_ox_*.json"))
        assert len(pn_files) > 0, "pNnm units should create files"
        assert len(ox_files) > 0, "oxDNA units should create files"

        # Compare the same field between runs; use 'total' for reliability.
        pn_by_field = {f.stem.split('energies_pn_', 1)[1]: f for f in pn_files}
        ox_by_field = {f.stem.split('energies_ox_', 1)[1]: f for f in ox_files}
        common_fields = set(pn_by_field) & set(ox_by_field)
        assert len(common_fields) > 0, "Both runs should produce the same fields"

        check_field = 'total' if 'total' in common_fields else sorted(common_fields)[0]
        with open(pn_by_field[check_field]) as f:
            pn_data = list(json.load(f).values())[0]
        with open(ox_by_field[check_field]) as f:
            ox_data = list(json.load(f).values())[0]

        pn_mean = np.mean(np.abs(pn_data))
        ox_mean = np.mean(np.abs(ox_data))
        assert ox_mean > 0, f"'{check_field}' should have non-zero energies in oxDNA units"
        ratio = pn_mean / ox_mean
        assert 30 < ratio < 55, f"pNnm/oxDNA ratio should be ~41.42, got {ratio}"
