# Tests of the CUDA implementation of the NA interaction

This directory checks that `interaction_type = NA` (oxDNA2 + oxRNA2 + DNA/RNA hybrid, `CUDADRHInteraction`) gives the same results on the GPU as on the CPU, with the neighbour-matrix kernel and with the edge-based one (`use_edge = true`).

Everything here needs an NVIDIA GPU and an oxDNA binary compiled with `-DCUDA=ON`. The same binary is used for the CPU and for the CUDA runs.

## Automated tests: `test_cuda_na.py`

```bash
# plain python
python3 test_cuda_na.py /path/to/oxDNA [-v]
# or pytest
OXDNA_BIN=/path/to/oxDNA pytest test_cuda_na.py
```

`OXDNA_CUDA_PRECISION` (`float`, `mixed`, the default, or `double`) selects the precision of the CUDA backend; `double` requires a binary built with `-DCUDA_DOUBLE=ON` and skips the edge tests, since `use_edge` is not available in double precision. If the binary cannot use a GPU the tests are skipped. The whole suite takes a few seconds.

The tests run on four small systems, which were relaxed and sampled with the CPU code: a DNA duplex (`dna_duplex`), an RNA duplex (`rna_duplex`), a DNA:RNA hybrid with a single-stranded overhang (`hybrid_overhang`) and a 179-nucleotide system with DNA and RNA duplexes, overhangs, hairpins and single strands (`mixed`). Each directory contains the topology (`init.top`) and the configuration (`init.dat`).

| Test | What it checks |
|---|---|
| `test_energy_cpu_vs_gpu` | the potential energy computed on the GPU (`CUDA_print_energy`) equals the CPU one, for the neighbour-matrix kernel, for `use_edge = true` and for `use_edge = true` with `edge_n_forces = 4` |
| `test_nve_trajectories` | 300-step NVE runs started from the same configuration on CPU and GPU end up at the same positions and total energy (this tests forces and torques), in the same three modes |
| `test_edge_forces_match_matrix_forces` | the per-nucleotide forces, torques and energies dumped (`CUDA_dump_file`) with the edge-based kernels equal those of the neighbour-matrix kernel up to round-off |
| `test_external_force` | a CUDA external force (trap) acts on NA nucleotides as it does on the CPU |
| `test_unsupported_options_fail` | options that are not supported (the stress tensor) make oxDNA exit with an error |

The tolerances are explained in the comments at the top of the file. The CPU code uses interpolated meshes for some angular modulations, whereas the GPU code uses the analytic functions, so energies differ by about 1e-5 and forces by up to about 5e-3 of the largest force even in double precision.

## Validation scripts: `validation/`

These are the tools that were used to validate the implementation more thoroughly. They are not run automatically (`parity.py` needs the Python bindings, `oxpy`, for the CPU reference) and they read their systems from a data directory (`NA_VALIDATION_DATA`, default `./data`) containing, for each system, `init.top`, `init.dat` and, optionally, `traj.dat` (a trajectory of a CPU run). Results go to `NA_VALIDATION_RUNS` (default `./runs`). The four directories next to this file can be used as a data directory:

```bash
export NA_VALIDATION_DATA=$PWD NA_VALIDATION_RUNS=/some/scratch/dir
```

* `parity.py`: compares per-nucleotide forces, body-frame torques and energies computed by CUDA (dumped with `CUDA_dump_file`, after the first force evaluation) and by the CPU code (`cpu_ref.py`, through `oxpy`, summing over all pairs), and prints the largest relative deviations. `--set "key = value"` adds an option to both inputs, so `--set "use_edge = true"` checks the edge-based kernels; `--seqdep` uses the sequence-dependent parameters, `--hybrid-only` restricts the comparison to the DNA/RNA terms (`CUDA_NA_dump_hybrid_only`), and `--terms` also compares the hydrogen-bonding energy.
  ```bash
  python3 validation/parity.py --bin /path/to/oxDNA --precision mixed --oxpy /path/containing/oxpy \
      --systems dna_duplex,rna_duplex,hybrid_overhang,mixed
  ```
* `nve_run.py`: energy conservation versus time step in an NVE run, on CPU or CUDA. The potential energy is computed on the CPU from the (GPU) positions, so forces that are not the gradient of the energy show up as a drift of the total energy.
  ```bash
  python3 validation/nve_run.py --system mixed --backend CUDA --precision double --dts 0.001,0.002 \
      --bin /path/to/oxDNA --tag test --out nve_test.json
  ```
* `thermo_run.py`: thermostatted run on CPU or CUDA, which writes averages (with block errors) of the kinetic and potential energy, of the energy of each interaction term and of the number of hydrogen bonds to a json file.
  ```bash
  mkdir -p res/thermo
  python3 validation/thermo_run.py --system mixed --T 50C --steps 2000000 --seed 1 --backend CUDA --bin /path/to/oxDNA \
      --tag cuda_mixed --out res/thermo/cuda_mixed_s1.json
  ```
* `thermo_summary.py`: run from the directory that contains `res/thermo/`, it averages the runs of `thermo_run.py` over seeds and gives the deviation of each CUDA precision from the CPU in units of the standard error.

The `CUDA_dump_file`, `CUDA_dump_call` and `CUDA_NA_dump_hybrid_only` options used by `parity.py` are meant for debugging and validation only.

## Reference results

Obtained on an NVIDIA GeForce RTX 3060 Laptop GPU (compute capability 8.6) with a Release build of this branch.

* `pytest test_cuda_na.py`: 5 passed (28 subtests) with `OXDNA_CUDA_PRECISION=mixed` and with `float`. With a `-DCUDA_DOUBLE=ON` build and `OXDNA_CUDA_PRECISION=double`, 4 passed and 1 skipped (the `use_edge` test, which is not available in double precision).
* `validation/parity.py` (four systems, mixed precision): the per-nucleotide deviations from the CPU code are the same with and without `use_edge = true`, and are due to the interpolated meshes used by the CPU code (see above).
* Sanity check on a 4 833-nucleotide system: 300 000-step thermostatted runs stay finite, and the potential energy per nucleotide agrees within the noise between the neighbour-matrix kernel and the edge-based one (also with `CUDA_sort_every = 10`).
* Accuracy against the CPU code (mixed precision, 300-step NVE runs, `dt = 0.001`, same results with the neighbour-matrix kernel and with `use_edge = true`):

| System | Relative difference of the potential energy (step 1) | Largest position difference after 300 steps (simulation units) | Total energy per nucleotide after 300 steps (CPU / GPU) |
|---|---:|---:|---|
| `dna_duplex` | 8e-7 | 1e-4 | -0.897328 / -0.897349 |
| `rna_duplex` | 7e-7 | 2e-6 | -1.091487 / -1.091486 |
| `hybrid_overhang` | 7e-7 | 3e-6 | -1.032886 / -1.032887 |
| `mixed` | 8e-7 | 1e-5 | -0.928045 / -0.928044 |

* Time per step on one CPU core (single-threaded oxDNA) and on the GPU, mixed precision, thermostatted (`john`) systems obtained by tiling the `mixed` system (1, 8, 27, 125 and 343 copies), measured as the wall-clock difference between a long and a short run, so that start-up costs are excluded:

| Nucleotides | CPU | GPU, `use_edge = false` | GPU, `use_edge = true` | CPU / GPU (edge) |
|---:|---:|---:|---:|---:|
| 179 | 0.38 ms | 0.20 ms | 0.071 ms | 5 |
| 1 432 | 3.0 ms | 0.20 ms | 0.070 ms | 43 |
| 4 833 | 10 ms | 0.22 ms | 0.11 ms | 91 |
| 22 375 | 55 ms | 0.74 ms | 0.26 ms | 210 |
| 61 397 | ~190 ms | 2.1 ms | 0.6 ms | ~320 |

The CPU timing at 61 397 nucleotides is noisy (188-262 ms in three runs, since the CPU was shared with other jobs), and the GPU one varied between 0.54 and 0.66 ms. These numbers are a reference for regressions, not a benchmark: they come from a single GPU and a single CPU.
