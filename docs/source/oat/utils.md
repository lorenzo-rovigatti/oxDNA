# OAT utilities documentation

## Data structures

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    UTILS.data_structures.Chunk
    UTILS.data_structures.ConfInfo
    UTILS.data_structures.TrajInfo
    UTILS.data_structures.Configuration
    UTILS.data_structures.TopInfo
    UTILS.data_structures.System
    UTILS.data_structures.Strand 
    UTILS.data_structures.Monomer
    
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.Chunk
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.ConfInfo
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.TrajInfo
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.Configuration
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.TopInfo
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.System
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.Strand
.. autoclass:: oxDNA_analysis_tools.UTILS.data_structures.Monomer

```

## Geometry utilities

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    UTILS.geom.fit_plane
    UTILS.geom.get_RNA_axis
    UTILS.geom.get_DNA_axis

.. autofunction:: oxDNA_analysis_tools.UTILS.geom.fit_plane
.. autofunction:: oxDNA_analysis_tools.UTILS.geom.get_RNA_axis
.. autofunction:: oxDNA_analysis_tools.UTILS.geom.get_DNA_axis
```

## iPython oxView plugin

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    UTILS.oxview.display_files
    UTILS.oxview.from_path
    UTILS.oxview.oxdna_conf
    UTILS.oxview.loro_patchy_conf
    UTILS.oxview.flro_patchy_conf

.. autofunction:: oxDNA_analysis_tools.UTILS.oxview.display_files
.. autofunction:: oxDNA_analysis_tools.UTILS.oxview.from_path
.. autofunction:: oxDNA_analysis_tools.UTILS.oxview.oxdna_conf
.. autofunction:: oxDNA_analysis_tools.UTILS.oxview.loro_patchy_conf
.. autofunction:: oxDNA_analysis_tools.UTILS.oxview.flro_patchy_conf
```

## Multiprocesser

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    UTILS.oat_multiprocesser.oat_multiprocesser

.. autofunction:: oxDNA_analysis_tools.UTILS.oat_multiprocesser.oat_multiprocesser
```

## Rye reader

```{eval-rst}
.. toctree::
   :maxdepth: 2

.. currentmodule:: oxDNA_analysis_tools

.. autosummary::
    :nosignatures:

    UTILS.RyeReader.Chunker
    UTILS.RyeReader.linear_read
    UTILS.RyeReader.get_confs
    UTILS.RyeReader.get_top_info
    UTILS.RyeReader.get_top_info_from_traj
    UTILS.RyeReader.get_traj_info
    UTILS.RyeReader.describe
    UTILS.RyeReader.strand_describe
    UTILS.RyeReader.get_input_parameter
    UTILS.RyeReader.inbox
    UTILS.RyeReader.write_conf
    UTILS.RyeReader.conf_to_str
    UTILS.RyeReader.get_top_string

.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.Chunker
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.linear_read
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.get_confs
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.get_top_info
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.get_top_info_from_traj
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.get_traj_info
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.describe
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.strand_describe
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.get_input_parameter
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.inbox
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.write_conf
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.conf_to_str
.. autofunction:: oxDNA_analysis_tools.UTILS.RyeReader.get_top_string

