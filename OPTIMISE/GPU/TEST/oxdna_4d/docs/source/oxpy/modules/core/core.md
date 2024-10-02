# The core module

This module contains the Python bindings of the C++ code. Every `object` described here is contained in the same module and can be accessed as `oxpy.core.object` or as `oxpy.object`.

```{eval-rst}
.. toctree::
   :maxdepth: 2
   
   forces.md
   observables.md

.. currentmodule:: oxpy.core

.. autosummary::
    :nosignatures:

    Context
    OxpyManager
    BaseParticle
    DNANucleotide
    RNANucleotide
    Molecule
    BaseInteraction
    ConfigInfo
    FlattenedConfigInfo
    FlattenedVectorArray
    BaseBox
    InputFile
    
.. autoclass:: Context
    
.. autoclass:: OxpyManager
    :inherited-members:
    
.. autoclass:: BaseParticle

.. autoclass:: DNANucleotide

.. autoclass:: RNANucleotide

.. autoclass:: Molecule

.. autoclass:: BaseInteraction

.. autoclass:: ConfigInfo

.. autoclass:: FlattenedConfigInfo

.. autoclass:: FlattenedVectorArray

.. autoclass:: BaseBox

.. autoclass:: InputFile

```
