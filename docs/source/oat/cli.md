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

## ANM fluctuation solver

Builds the elastic network of a protein for the ANM-based interactions (`ANM`, `ANMT`, `DNANM`, `DNACT`, `RNANM`, `RNACT`) from a PDB or mmCIF file. Every pair of C-alpha atoms closer than the cutoff is joined by a spring, and the spring constant is fitted so that the B-factors predicted by normal mode analysis match the experimental ones (by default the C-alpha B-factors of the structure). It writes the parameter file (`<output>.par`), a new-style topology (`<output>.top`, one `type=peptide` strand per chain) and a starting configuration (`<output>.dat`).

With `--model anmt`, every spring between consecutive particles also carries the equilibrium cosines used by the bending/torsion potentials of the DNACT/RNACT interactions. The bending and torsion terms (`--kb`/`--kt`, default 50) are included in the normal mode analysis and only the spring constant is fitted, so the same values must be used in the simulation: either set `bending_k` and `torsion_k` in the input file or use `--write-kbkt` to write them on every backbone line of the parameter file. The particle indices in the parameter file are global indices: if the protein is not the first molecule of a hybrid system, set `--index-offset` to the number of particles that precede it. `--plot` compares the predicted and experimental B-factors.

For more information, see [the ANM model publication](https://doi.org/10.1039/d0sm01639j).

```bash
oat anm_fluctuation_solver protein.pdb -o protein -m anmt --plot protein_bfactors.png
```

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/anm_fluctuation_solver.py
    :func: cli_parser
    :prog: oat anm_fluctuation_solver
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

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/IDconvert.py
    :func: cli_parser
    :prog: oat IDconvert
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

When converting protein-containing structures, reference PDB/mmcif files for the proteins must be provided. The ANM model does not contain any orientation information, so a crude alignment between the reference all-atom files and the ANM model is performed before translating the Cα atoms to match the coarse-grained bead positions. We do not recommend using proteins converted from the ANM model for rigirous analysis which relies on sidechains or all-atom simulations.

For nucleic acid structures converted to all-atom formats using `oxDNA_PDB`, it is often beneficial to optimize the all-atom geometry using a pairing-aware tool like [QRNAS](https://genesilico.pl/software/stand-alone/qrnas) before continuing with analysis or further simulation.

### All-atom formats
`oxDNA_PDB` contains writers for both the classic PDB the modern mmcif format. For structures with over 100000 atoms or 1000 strands, the mmcif format is automatically generated due to limitations in the PDB format's fixed-width columns. If you want to force a classic PDB format file with hybrid-36 encoding for column overflows, use the `--format pdb` flag.

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

## Skeleton

```{eval-rst}
.. argparse::
    :filename: ../analysis/src/oxDNA_analysis_tools/skeleton.py
    :func: cli_parser
    :prog: oat skeleton
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