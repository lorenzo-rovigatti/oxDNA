# oxDNA

oxDNA is a simulation code that was initially conceived as an implementation of the coarse-grained DNA model introduced by [T. E. Ouldridge, J. P. K. Doye and A. A. Louis](http://dx.doi.org/10.1063/1.3552946). It has been since reworked and it is now an extensible simulation+analysis framework. 

oxDNA can perform both molecular dynamics (MD) and Monte Carlo (MC) simulations of the oxDNA and oxRNA models. MD simulations can be run on single CPUs or single CUDA-enabled GPUs, while MC simulations, which can only be run serially, can exploit the Virtual Move Monte Carlo algorithm[@VMMC] to greatly speed-up equilibration and sampling, and Umbrella Sampling biasing to efficiently obtain free-energy profiles. The package also features a Forward-Flux Sampling interface to study the kinetics of rare events, and makes it possible to alter the behaviour of the systems by adding *external forces* that can be used, for instance, to pull on or apply torques to strands or confine nucleotides within semi-planes or spheres.

The package also includes `oxpy`, a Python library which makes it possible to control the behavior of the simulation using Python scripts. The repository contains examples that demonstrate how to leverage `oxpy` to write backends to run replica-exchange and well-tempered metadynamics simulations, which are popular techniques in modern molecular dynamics to improve sampling efficiency.

The package can also be used to compute nucleic-acid-related quantities such as the energy due to hydrogen bonding or stacking, the distance between groups of nucleotides, the list of hydrogen-bonded nucleotides or of over-stretched bonds, and much more. The analysis can be performed while the simulation is running or on generated trajectory files.

The simulation engine is complemented by an updated version of `oxDNA_analysis_tools` (`oat`), a Python library aimed at facilitating the analysis of oxDNA/oxRNA trajectories. `Oat` provides numerous common simulation trajectory analysis tools including alignment, mean structures, subsetting trajectories, distances between nucleotides, interduplex angles, and comparisons in hydrogen bonding patterns between the trajectory and an idealized structure.

```{eval-rst}
.. toctree::
   :caption: Table of Contents
   :maxdepth: 2
   
   install.md
   usage.md
   input.md
   configurations.md
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
