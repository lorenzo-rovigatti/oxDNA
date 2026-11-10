# oxDNA Analysis Tools (oat) instructions

This file provides guidance to coding agents working on the `analysis/` directory of the oxDNA repository. The instructions in the top-level `AGENTS.md` also apply.

## Project Overview

oxDNA Analysis Tools (`oat`) is a Python package for structural analysis of oxDNA coarse-grained DNA/RNA simulation trajectories. It provides both a CLI (`oat <command>`) and importable Python library.

**Documentation:** https://lorenzo-rovigatti.github.io/oxDNA/oat/index.html built from ../docs/source/oat/

## Build & Installation

```bash
# Install the whole package (Python modules + compiled Cython extensions), from the analysis directory
python -m pip install . --no-build-isolation   # --no-build-isolation avoids numpy header errors

# Build only the Cython utilities (cython_utils/) in place, without reinstalling the package
python -m build

# Verify installation
oat config
```

Use `python -m build` when iterating on `.pyx` files; use `pip install` to install or update the full pipeline (including the `oat` entry point). Either way, `.pyx` changes are not picked up until rebuilt.

## Running Tests

```bash
# Run pytest test suite
pytest tests/test_cli/

# Run a single test file
pytest tests/test_cli/test_mean.py

# Check code coverage 
coverage run --source=oxDNA_analysis_tools -m pytest tests/test_cli/
coverage combine
coverage report

# Run legacy bash test suite (from tests/ directory)
cd tests && ./test.sh
```

## Architecture

### CLI Entry Point
`src/oxDNA_analysis_tools/entry.py` - The `oat` command dynamically imports and executes analysis scripts. Running `oat` with no arguments lists all available commands.

### Analysis Script Pattern
Each analysis script (e.g., `mean.py`, `align.py`, `pca.py`) follows this structure:
1. `ComputeContext` namedtuple - holds parameters for parallel execution
2. `compute(ctx, chunk_size, chunk_id)` - parallelized function for processing trajectory chunks
3. Main function (same name as file, e.g., `mean()`) - importable API for the analysis
4. `cli_parser()` - returns argparse parser (used for auto-generating docs)
5. `main()` - CLI entry point, must take no arguments

### Core Modules (UTILS/)
- **RyeReader.py** - File I/O: `describe()` parses metadata, `get_confs()` reads configurations, `write_conf()` writes output
- **data_structures.py** - Core types: `Configuration`, `TopInfo`, `TrajInfo`
- **oat_multiprocesser.py** - Parallelization framework using chunk-based processing
- **geom.py** - Geometry/alignment utilities
- **oxview.py** - Generates JSON overlays for oxView visualization

### Cython Extensions
`cython_utils/get_confs.pyx` contains performance-critical trajectory reading code. Changes require rebuilding.

### Compressed (zstd) Trajectories
`RyeReader` transparently reads zstd-compressed trajectories in addition to plain text ones; `describe()`/`get_confs()` dispatch on `_is_zstd_traj()`, and `write_conf(..., compress=True, level=3)` writes them.

File layout: `[8-byte header][zstd frame 0][zstd frame 1]...`, one independently compressed frame per configuration. Header = magic `OXD\x01` + version byte + compression level byte + 2 reserved bytes.

- `_index_compressed()` builds the frame index by seeking: `_zstd_frame_size_from_file()` parses each zstd frame/block header to get the frame size without reading the payload, so indexing large files stays streaming-safe. Don't replace this with `zstandard.stream_reader(read_across_frames=False)` — it does not reliably stop at frame boundaries.
- `_get_confs_compressed()` decompresses only the requested frames and hands the bytes to `cget_confs_zst()` in `get_confs.pyx`, which reuses the same `parse_conf` as the text reader.
- zstd Frame Header Descriptor bit layout (easy to get wrong): bits 7–6 = Frame_Content_Size_Flag, bit 5 = Single_Segment_Flag, bit 2 = Content_Checksum_Flag, bits 1–0 = Dict_ID_Flag.

### Documenting new tools
Every new command must be documented in both places, in alphabetical order with the existing entries:
- the CLI documentation, `docs/source/oat/cli.md`: an `argparse` directive that renders the command's `cli_parser()`:
  ```
  .. argparse::
      :filename: ../analysis/src/oxDNA_analysis_tools/<script>.py
      :func: cli_parser
      :prog: oat <script>
  ```
- the API documentation, `docs/source/oat/api.md`: an `autosummary` and an `autofunction` entry for the main function and for every other public function of the module, whose docstrings must follow the existing format (description, `Parameters`, `Returns`).

Also add the command to the `opts` list in `oat-completion.sh`. New options of existing commands only need their argparse `help` text and docstrings to be updated. Build the documentation (see the top-level `AGENTS.md`) to check that the new entries render.

## Key Conventions
- Trajectories are processed in chunks for memory efficiency and parallelization
- The `-p` flag on most commands sets CPU count for parallel processing
- Index files are space-separated particle IDs (generated by oxView's "Download Selected Base List")
- Output JSON files are overlays for oxView visualization
- There are two possible topology files which readers must handle. Old-style oxDNA files are written 3'-5' (opposite of traditional 5'-3'), while 'new-style' topologies are written 5'-3'. The correct topology reader (from RyeReader) must always be used when handling files.
- When writing tests, prefer end-to-end tests where the output of a single function call is put through increasingly rigorous consistency and correctness checking rather than calling the function multiple times for each check.
- There will be times where multiple function calls are necessary, for example, checking that a function works with both old-style and new-style topologies.
- Whenever possible when writing tests, verify correctness of output content, not just shape and data type of output.

## Dependencies
- Python >=3.11 but must work with the latest python
- numpy, scikit-learn, matplotlib, zstandard
- oxpy (oxDNA Python bindings, >=3.2.2)
