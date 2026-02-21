# OAT command line documentation

## Align

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/align.py
    :func: cli_parser
    :prog: oat align
```

## ANM parameterize

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/anm_parameterize.py
    :func: cli_parser
    :prog: oat anm_parameterize
```

## Backbone flexibility

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/backbone_flexibility.py
    :func: cli_parser
    :prog: oat backbone_flexibility
```

## Bond analysis

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/bond_analysis.py
    :func: cli_parser
    :prog: oat bond_analysis
```

## Centroid

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/centroid.py
    :func: cli_parser
    :prog: oat centroid
```

## Clustering

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/clustering.py
    :func: cli_parser
    :prog: oat clustering
```

## Config

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/config.py
    :func: cli_parser
    :prog: oat config
```

## Contact map

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/contact_map.py
    :func: cli_parser
    :prog: oat contact_map
```
## Decimate

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/decimate.py
    :func: cli_parser
    :prog: oat decimate
```

## Deviations

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/deviations.py
    :func: cli_parser
    :prog: oat deviations
```

## Distance

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/distance.py
    :func: cli_parser
    :prog: oat distance
```

## Dot-bracket to forces

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/db2forces.py
    :func: cli_parser
    :prog: oat db_to_forces
```

## Duplex angle plotter

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/duplex_angle_plotter.py
    :func: cli_parser
    :prog: oat duplex_angle_plotter
```

## Duplex finder

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/duplex_finder.py
    :func: cli_parser
    :prog: oat duplex_finder
```

## File info
```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/file_info.py
    :func: cli_parser
    :prog: oat file_info
```

## Forces to dot-bracket
```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/forces2db.py
    :func: cli_parser
    :prog: oat forces2db
```

## Forces to pairs

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/forces2pairs.py
    :func: cli_parser
    :prog: oat forces2pairs
```

## Generate force

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/generate_forces.py
    :func: cli_parser
    :prog: oat generate_forces
```

## IDconvert

Translate nucleotide IDs between two oxDNA topology files.

```
oat IDconvert <reference.top> <compare.top> <ids.txt>
```

**Positional arguments:**

| Argument | Description |
|---|---|
| `reference.top` | Topology file whose nucleotide IDs are provided as input. |
| `compare.top` | Topology file to map the IDs into. |
| `ids.txt` | Plain-text file containing a comma-separated list of nucleotide IDs from the reference topology. |

Both old-style (`<N_nucs> <N_strands>`) and new-style (`<N_nucs> <N_strands> 5->3`) topology formats are supported in any combination.

**Output:**

- **stdout** – comma-separated list of the corresponding nucleotide IDs in the compare topology. Nucleotides present in the reference but absent from the compare file are listed as `DELETED`.
- **stderr** – human-readable summary: nucleotide and strand counts for both files, per-strand sequence-similarity scores, and the total number of deleted/added nucleotides.

**Example:**

```bash
oat IDconvert before.top after.top selected_ids.txt
```

## Mean

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/mean.py
    :func: cli_parser
    :prog: oat mean
```

## Minify

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/minify.py
    :func: cli_parser
    :prog: oat minify
```

## Multidimensional scaling mean

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/multidimensional_scaling_mean.py
    :func: cli_parser
    :prog: oat multidimensional_scaling_mean
```

## Output bonds

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/output_bonds.py
    :func: cli_parser
    :prog: oat output_bonds
```

## oxDNA -> PDB

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/oxDNA_PDB.py
    :func: cli_parser
    :prog: oat oxDNA_PDB
```

## Pairs to dot-bracket
```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/pairs2db.py
    :func: cli_parser
    :prog: oat pairs2db
```

## PDB -> oxDNA

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/PDB_oxDNA.py
    :func: cli_parser
    :prog: oat PDB_oxDNA
```

## Persistence Length

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/persistence_length.py
    :func: cli_parser
    :prog: oat persistence_length
```

## Plot energy

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/plot_energy.py
    :func: cli_parser
    :prog: oat plot_energy
```

## Principal component analysis

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/pca.py
    :func: cli_parser
    :prog: oat pca
```

## Subset trajectory

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/subset_trajectory.py
    :func: cli_parser
    :prog: oat subset_trajectory
```

## Superimpose

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/superimpose.py
    :func: cli_parser
    :prog: oat superimpose
```