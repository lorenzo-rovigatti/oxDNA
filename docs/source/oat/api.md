# OAT scripting interface documentation

Each script has been broken down into two types of functions: 

1. An importable module with the same name as the corresponding command line invocation which runs the entire script and returns the result as a Python object.
2. A series of sub-modules which implement individual steps of the overall process.

The first type of import can be used to compose multiple analysis types into a single script or Jupyter notebook.  For example, to first align a trajectory, then decimate the resulting trajectory to contain only every 10th trajectory while skipping the first 200 confgiruations, each using 5 processes, you would run:

```python
from oxDNA_analysis_tools.align import align
from oxDNA_analysis_tools.decimate import decimate

traj = 'traj.dat'
align_output = 'aligned.dat'
decimate_output = 'decimated.dat'

align(traj, align_output, ncpus=5)
decimate(align_output, decimate_output, ncpus=5, start=200, stride=10)
```

## Align

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    align.align
    align.svd_align
    
.. autofunction:: oxDNA_analysis_tools.align.align
.. autofunction:: oxDNA_analysis_tools.align.svd_align
```

## ANM parameterize

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    anm_parameterize.anm_parameterize
    
.. autofunction:: oxDNA_analysis_tools.anm_parameterize.anm_parameterize
```

## Backbone flexibility

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    backbone_flexibility.backbone_flexibility
    
.. autofunction:: oxDNA_analysis_tools.backbone_flexibility.backbone_flexibility
```

## Bond analysis

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    bond_analysis.bond_analysis
    
.. autofunction:: oxDNA_analysis_tools.bond_analysis.bond_analysis
```

## Centroid

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    centroid.centroid
    
.. autofunction:: oxDNA_analysis_tools.centroid.centroid
```

## Clustering

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    clustering.split_trajectory
    clustering.get_centroid
    clustering.perform_DBSCAN
    
.. autofunction:: oxDNA_analysis_tools.clustering.split_trajectory
.. autofunction:: oxDNA_analysis_tools.clustering.get_centroid
.. autofunction:: oxDNA_analysis_tools.clustering.perform_DBSCAN
```

## Config

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    config.check
    config.set_chunk_size
    config.get_chunk_size
    
.. autofunction:: oxDNA_analysis_tools.config.check
.. autofunction:: oxDNA_analysis_tools.config.set_chunk_size
.. autofunction:: oxDNA_analysis_tools.config.get_chunk_size
```

## Contact map

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

contact_map.contact_map
    
.. autofunction:: oxDNA_analysis_tools.contact_map.contact_map
```

## Decimate

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

decimate.decimate
    
.. autofunction:: oxDNA_analysis_tools.decimate.decimate
```

## Deviations

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

deviations.deviations
deviations.output
    
.. autofunction:: oxDNA_analysis_tools.deviations.deviations
.. autofunction:: oxDNA_analysis_tools.deviations.output
```

## Distance

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

distance.min_image
distance.vectorized_min_image
distance.distance
    
.. autofunction:: oxDNA_analysis_tools.distance.min_image
.. autofunction:: oxDNA_analysis_tools.distance.vectorized_min_image
.. autofunction:: oxDNA_analysis_tools.distance.distance
```

## Dot-bracket to force

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

db2forces.parse_dot_bracket
db2forces.db_to_forcelist
    
.. autofunction:: oxDNA_analysis_tools.db2forces.parse_dot_bracket
.. autofunction:: oxDNA_analysis_tools.db2forces.db_to_forcelist
```

## Duplex angle plotter

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

duplex_angle_plotter.get_angle_between
    
.. autofunction:: oxDNA_analysis_tools.duplex_angle_plotter.get_angle_between
```

## Duplex finder

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

duplex_finder.Duplex
duplex_finder.find_duplex
duplex_finder.duplex_finder
    
.. autoclass:: oxDNA_analysis_tools.duplex_finder.Duplex
.. autofunction:: oxDNA_analysis_tools.duplex_finder.find_duplex
.. autofunction:: oxDNA_analysis_tools.duplex_finder.duplex_finder
```

## File info

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

file_info.file_info
    
.. autofunction:: oxDNA_analysis_tools.file_info.file_info
```

## Forces to dot-bracket

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

forces2db.forces2db
    
.. autofunction:: oxDNA_analysis_tools.forces2db.forces2db
```

## Mean

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

mean.mean
    
.. autofunction:: oxDNA_analysis_tools.mean.mean
```

## Minify

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

minify.minify
    
.. autofunction:: oxDNA_analysis_tools.minify.minify
```

## Multidimensional scaling mean

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

multidimensional_scaling_mean.multidimensional_scaling_mean
multidimensional_scaling_mean.distance_deviations
    
.. autofunction:: oxDNA_analysis_tools.multidimensional_scaling_mean.multidimensional_scaling_mean
.. autofunction:: oxDNA_analysis_tools.multidimensional_scaling_mean.distance_deviations
```

## Output bonds

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

output_bonds.output_bonds
    
.. autofunction:: oxDNA_analysis_tools.output_bonds.output_bonds
```

## OxDNA to PDB

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

oxDNA_PDB.oxDNA_PDB
oxDNA_PDB.choose_reference_nucleotides
    
.. autofunction:: oxDNA_analysis_tools.oxDNA_PDB.oxDNA_PDB
.. autofunction:: oxDNA_analysis_tools.oxDNA_PDB.choose_reference_nucleotides
```

## Pairs to dot-bracket

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

pairs2db.pairs2db
    
.. autofunction:: oxDNA_analysis_tools.pairs2db.pairs2db
```

## PDB to oxDNA
```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

PDB_oxDNA.PDB_oxDNA
    
.. autofunction:: oxDNA_analysis_tools.PDB_oxDNA.PDB_oxDNA
```

## Persistence length

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

persistence_length.persistence_length
persistence_length.get_r
persistence_length.fit_PL
    
.. autofunction:: oxDNA_analysis_tools.persistence_length.persistence_length
.. autofunction:: oxDNA_analysis_tools.persistence_length.get_r
.. autofunction:: oxDNA_analysis_tools.persistence_length.fit_PL
```

## Principle component analysis

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

pca.align_positions
pca.map_confs_to_pcs
pca.make_heatmap
pca.pca
    
.. autofunction:: oxDNA_analysis_tools.pca.align_positions
.. autofunction:: oxDNA_analysis_tools.pca.map_confs_to_pcs
.. autofunction:: oxDNA_analysis_tools.pca.make_heatmap
.. autofunction:: oxDNA_analysis_tools.pca.pca
```

## Subset trajectory

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

subset_trajectory.subset
    
.. autofunction:: oxDNA_analysis_tools.subset_trajectory.subset
```

## Superimpose

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

superimpose.superimpose
    
.. autofunction:: oxDNA_analysis_tools.superimpose.superimpose
```