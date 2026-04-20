# ANNaMo utilities

This directory contains command-line tools to prepare and run ANNaMo simulations.

In addition to `annamo.py`, which generates all input files starting from a bead-level JSON description, this directory also includes a bead-division utility to automatically convert nucleotide-level systems into ANNaMo bead representations.

## Requirements

- Python 3.8+
- oxDNA compiled from source (see the main oxDNA repository)
- `oxDNA` and `confGenerator` available in PATH (found in `build/bin/` after compilation)

## Usage

### ANNaMo simulation setup

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

## Automatic bead division

For systems specified at nucleotide resolution (e.g. from native base-pair lists), the `bead_division.py` script automatically groups nucleotides into ANNaMo beads.

Each bead represents a sequence of 2–4 nucleotides, with a preference for 3 nt, following the ANNaMo coarse-graining scheme.

### Basic usage

```bash
# Bead division from pairing file only
python bead_division.py pairs.dat

# Strand-aware bead division using an oxDNA topology file
python bead_division.py pairs.dat -t topology.top

# Generate ANNaMo JSON and native bead bonds
python bead_division.py pairs.dat -t topology.top \
  --json-out system.json \
  --bonds-out native_bonds.txt
```

## JSON template support

When generating a JSON file for ANNaMo, it is often useful to reuse simulation parameters
(e.g. temperature, steps, box size, etc.) from an existing configuration.

This can be done using:

```bash
python bead_division.py pairs.dat -t topology.top \
  --json-out system.json \
  --json-template template.json
```

The template file should be a valid ANNaMo JSON input. All fields from the template are preserved,
except for:

- `strands`, which is always replaced by the newly generated bead strands.

This allows you to:
- keep simulation parameters consistent across systems,
- avoid rewriting boilerplate JSON,
- separate geometry (bead division) from simulation setup.

## Pairing file format

The pairing file describes native base pairs at nucleotide resolution.

Each line can be:

```
i
```

for an unbound nucleotide, or:

```
i j
```

for a native base pair between nucleotides i and j.

The file supports multiple formats:

- complete (all nucleotides explicitly listed),
- sparse (unbound nucleotides omitted),
- symmetric (`i j` and `j i` both present),
- non-symmetric (each pair listed only once).

The script automatically:
- fills missing nucleotides as unbound,
- completes missing reverse pairs (`j i`),
- checks consistency if both directions are provided.

Contradictory definitions raise an error.

## Topology support

The script supports both:

- old oxDNA topology format (`N N_strands`)
- new oxDNA topology format (`N N_strands 5->3`)

When a topology file is provided:

- bead division becomes strand-aware: no bead spans multiple strands,
- strand sequences are reconstructed,
- output JSON is directly compatible with `annamo.py`.

Without a topology file:

- bead division is based only on pairing information,
- strand information is not available,
- only native bead bonds can be generated reliably.

## Bead design rules

Beads are constructed according to the ANNaMo coarse-graining scheme
(Tosti Guerra et al., *J. Chem. Phys.* 2024, DOI: 10.1063/5.0202829):

- bead size is ideally 3 nt, with 2 nt and 4 nt used when necessary;
- native base pairs are used to define bead boundaries;
- no bead contains nucleotides from different strands;
- each bead has at most one real partner bead;
- bead pairing preserves nucleotide ordering;
- native contacts are not split across multiple bead pairs.

The algorithm automatically handles irregular patterns such as bulges, loops, and partial pairing.

## Output

The bead-division script can generate:

- terminal output: bead list and partners (color-coded)
- `system.json`: ANNaMo-ready input file (`--json-out`)
- `native_bonds.txt`: bead-level native bonds (`--bonds-out`)

It also prints a final report including:

- total number of beads,
- counts of 2 nt / 3 nt / 4 nt beads,
- dimensional matches and mismatches between paired beads,
- warnings for possible pathological cases.

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

## Output files (ANNaMo)

| File | Description |
|---|---|
| `topology.top` | bead topology |
| `dHdS_matrix.dat` | nearest-neighbour interaction matrix |
| `input.dat` | oxDNA input file |
| `init_conf.dat` | random initial configuration |
| `trajectory.dat` | simulation trajectory (`run` only) |
| `energy.dat` | energy output (`run` only) |

## Workflow

A typical workflow is:

```bash
python bead_division.py pairs.dat -t topology.top \
  --json-template template.json \
  --json-out system.json

python annamo.py prepare system.json
python annamo.py run system.json
```

This allows starting from a nucleotide-level native-contact map and automatically generating an ANNaMo-ready simulation system while reusing predefined simulation parameters.
