# oxDNA

oxDNA is a simulation code that was initially conceived as an implementation of the coarse-grained DNA model introduced by [T. E. Ouldridge, J. P. K. Doye and A. A. Louis](http://dx.doi.org/10.1063/1.3552946). It has been since reworked and it is now an extensible simulation+analysis framework. The following nucleic-acid-related force fields are implemented:

* [oxDNA1](https://pubs.aip.org/aip/jcp/article-abstract/137/13/135101/191221/Sequence-dependent-thermodynamics-of-a-coarse?redirectedFrom=fulltext)
* [oxDNA2](https://pubs.aip.org/aip/jcp/article-abstract/142/23/234901/193554/Introducing-improved-structural-properties-and?redirectedFrom=fulltext)
* [oxRNA](https://pubs.aip.org/aip/jcp/article-abstract/140/23/235102/73414/A-nucleotide-level-coarse-grained-model-of-RNA?redirectedFrom=fulltext)
* [oxNA](https://pubs.aip.org/aip/jcp/article/160/11/115101/3275728)

oxDNA can perform both molecular dynamics (MD) and Monte Carlo (MC) simulations of the oxDNA and oxRNA models. MD simulations can be run on single CPUs or single CUDA-enabled GPUs, while MC simulations, which can only be run serially, can exploit the Virtual Move Monte Carlo algorithm to greatly speed-up equilibration and sampling, and [Umbrella Sampling biasing](umbrella_sampling.md) to efficiently obtain free-energy profiles. The package also features a [Forward-Flux Sampling interface](ffs.md) to study the kinetics of rare events, and makes it possible to alter the behaviour of the systems by adding [*external forces*](forces.md) that can be used, for instance, to pull on or apply torques to strands or confine nucleotides within semi-planes or spheres.

The package also includes [the `oxpy` Python module](oxpy/index.md), which makes it possible to control the behavior of the simulation using Python scripts. The repository contains examples that demonstrate how to leverage `oxpy` to write backends to run replica-exchange and well-tempered metadynamics simulations, which are popular techniques in modern molecular dynamics to improve sampling efficiency.

The package can be used to compute [nucleic-acid-related quantities](observables.md) such as the energy due to hydrogen bonding or stacking, the distance between groups of nucleotides, the list of hydrogen-bonded nucleotides or of over-stretched bonds, and much more. The analysis can be performed while the simulation is running or on generated trajectory files.

The simulation engine is complemented by an updated version of `oxDNA_analysis_tools` (`oat`), a Python library aimed at facilitating the analysis of oxDNA/oxRNA trajectories. `Oat` provides numerous common simulation trajectory analysis tools including alignment, mean structures, subsetting trajectories, distances between nucleotides, interduplex angles, and comparisons in hydrogen bonding patterns between the trajectory and an idealized structure.

## Citing oxDNA

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04693/status.svg)](https://doi.org/10.21105/joss.04693)

Please cite these publications for any work that uses the oxDNA simulation package:

- for the code:
  * [E. Poppleton et al., J. Open Source Softw. 8, 4693 (2023)](https://doi.org/10.21105/joss.04693)
- for the CUDA-powered code:
  * [L. Rovigatti et al., J. Comput. Chem. 36, 1 (2015)](https://doi.org/10.1002/jcc.23763)
- for oxDNA analysis tools:
  * [E. Poppleton et al., Nucleic Acids Research e72 (2020)](https://doi.org/10.1093/nar/gkab324)

```{eval-rst}
.. toctree::
   :caption: Table of Contents
   :maxdepth: 2
   
   install.md
   usage.md
   input.md
   configurations.md
   relaxation.md
   performance.md
   scaling.md
   forces.md
   observables.md
   umbrella_sampling.md
   ffs.md
   events.md
```

```{eval-rst}
.. ifconfig:: with_oxpy

Analysis
--------

.. toctree::
   :maxdepth: 3
   
   oat/index.md

Python Bindings
---------------

.. toctree::
   :maxdepth: 3
   
   oxpy/index.md
```
