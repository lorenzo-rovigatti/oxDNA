---
title: 'oxDNA: coarse-grained simulations of nucleic acids made simple'
tags:
  - C++
  - Python
  - DNA
  - DNA nanotechnology
  - Molecular simulations
authors:
  - name: Erik Poppleton
	orcid: 
    affiliation: 1
  - name: Michael Matthies
  	orcid: 
    affiliation: 1
  - name: Flavio Romano
  	orcid: 
  	affiliation: 2
  - name: Petr Šulc
   	orcid:
   	affiliation: 1
  - name: Lorenzo Rovigatti^[Corresponding author]
    orcid: 0000-0001-5017-2829
    affiliation: "3, 4"
affiliations:
 - name: School of Molecular Sciences and Center for Molecular Design and Biomimetics, The Biodesign Institute, Arizona State University, USA
   index: 1
 - name: Dipartimento di Scienze Molecolari e Nanosistemi, Universitá Ca Foscari di Venezia, Italy
   index: 2
 - name: Department of Physics, Sapienza University of Rome, Italy
   index: 3
 - name: CNR-ISC UoS Sapienza, Rome, Italy
   index: 4
date: 08 May 2022
bibliography: paper.bib
---

# Summary

The fields of DNA and RNA nanotechnology have progressed from pioneering, proof-of-principle experiments to fully-fledged applications in material science, biology and medicine. These applications exploit the intrinsic programmability of nucleic acids to generate nano- and even micro-scale structures with tailored properties. However, the design of the DNA/RNA sequences that self-assemble into a desired structure is not straightforward and often relies on expensive trial-and-error experimental protocols. A complementary approach is provided by computer simulations, which can model biomacromolecules at different levels of detail, ranging from atomistic to continuous, and can be leveraged to investigate the whole range of time- and length-scales relevant for applications.

# Statement of need

The oxDNA and oxRNA models, which describe nucleic acids at the nucleotide level, has become popular choices to investigate the dynamics, thermodynamics and self-assembly kinetics of DNA and RNA systems[@oxDNA_review; @oxRNA; @oxDNA_primer]. Here we present `oxDNA`, an efficient, multi-technique simulation package written in C++ and specifically developed to carry out simulations of coarse-grained nucleic acid systems. The package, which has been used in more than a hundred publications to date, can perform both molecular dynamics (MD) and Monte Carlo (MC) simulations of oxDNA and oxRNA. MD simulations can be run on single CPUs or single CUDA-enabled GPUs, while MC simulations, which can only be run serially, can exploit the Virtual Move Monte Carlo algorithm[@VMMC] to greatly speed-up equilibration and sampling, as well as Umbrella Sampling biasing to efficiently obtain free-energy profiles. It is also possible to alter the behaviour of the systems by adding *external forces* that can be used, for instance, to pull on or apply torques to strands or confine nucleotides within semi-planes or spheres. Finally, the code is equipped with Python bindings that makes it possible to control the behaviour of the simulation at a high level: The repository contains examples that shows how to leverage `oxDNA`'s Python bindings to write backends to run replica-exchange or metadynamics simulations.

The simulation engine is complemented by `oxDNA_analysis_tools` (`oat`), a Python library aimed at facilitating the analysis of oxDNA/oxRNA trajectories. Specifically `oat` ... @Erik, Michael, Petr, can you add a few details about what `oat` can do?

We note that molecular dynamics simulations of oxDNA and oxRNA can also be performed with the LAMMPS software package, whose efficient parallel-computing algorithms make it possible to simulate large systems on HPC clusters[@oxDNA_LAMMPS]. However, the generic and hyper-customizable nature of the LAMMPS package does not lend itself well to being used by non specialists, and it lacks many of the tools required to ease the burden of working with coarse-grained nucleic acids.


# Acknowledgements

We acknowledge contributions from Debesh Mandal, Ferdinando Randisi and Benedict E. K. Snodin. 

# References

