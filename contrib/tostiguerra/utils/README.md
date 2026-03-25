# ANNaMo utilities

This directory contains the `annamo` command-line tool, which automates the
setup and launch of ANNaMo simulations starting from a single JSON input file.

## Requirements

- Python 3.8+
- oxDNA compiled from source (see the main oxDNA repository)
- `oxDNA` and `confGenerator` available in PATH (found in `build/bin/` after compilation)

## Usage
```bash
# Generate all input files without launching the simulation
python annamo.py prepare system.json

# Generate all input files and launch the simulation
python annamo.py run system.json
```

After `prepare`, you can freely edit `input.dat` before launching manually:
```bash
oxDNA input.dat
```

## Input JSON

Only `material` and `strands` are required. All other fields have sensible defaults.
```json
{
  "material":             "DNA",
  "strands":              [["GAA", "GTG", "ACA"], ["TGT", "CAC", "TTC"]],
  "temperature":          30,
  "steps":                2e9,
  "box_size":             30,
  "swap":                 true,
  "seed":                 104123,
  "print_conf_interval":  1e6,
  "print_energy_every":   1e3,
  "oxdna_overrides":      {}
}
```

| Field | Default | Description |
|---|---|---|
| `material` | — | `"DNA"` or `"RNA"` |
| `strands` | — | list of strands, each a list of bead sequences |
| `temperature` | `30` | degrees C |
| `steps` | `2e9` | MD steps |
| `box_size` | `30` | simulation box side (internal units) |
| `swap` | `true` | `true` → λ=1 (bond-swapping); `false` → λ=10 (no swapping) |
| `seed` | random | RNG seed |
| `print_conf_interval` | `1e6` | steps between saved configurations |
| `print_energy_every` | `1e3` | steps between energy output |
| `oxdna_overrides` | `{}` | raw key=value pairs appended to `input.dat` |

## Bead design

Each bead represents a nucleotide sequence of ideally 3 nt (range 2–4).
Strand division should follow native contacts as described in the ANNaMo paper
(Tosti Guerra et al., *J. Chem. Phys.* 2024, DOI: 10.1063/5.0202829): beads are
chosen so that native base pairs fall at bead boundaries rather than within a
single bead.

## Output files

| File | Description |
|---|---|
| `topology.top` | bead topology |
| `dHdS_matrix.dat` | nearest-neighbour interaction matrix |
| `input.dat` | oxDNA input file |
| `init_conf.dat` | random initial configuration |
| `trajectory.dat` | simulation trajectory (`run` only) |
| `energy.dat` | energy output (`run` only) |
