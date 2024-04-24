# Persistence Length

Authors: Erik Popppleton & Petr Šulc  
Last updated: Jan 2023

This example shows how to calculate a persistence length of a DNA duplex. The configuration and topology files in the "classic" format are in this directory, while files with the "new" format can be found in the `NEW_TOPOLOGY` directory (see documentation for details).

### Running the simulation

This directory contains initial configuration files for a 201 base pair dsDNA segment.  Run the simulation with

```
oxDNA input_persistence
```

In order to compute the persistence length, one needs a large number of decorrelated states, so this simulation takes a while.  It can be sped up by running on a GPU rather than the default CPU.

### Analysis

There is a script for calculating the persistence length of a paired sequence of DNA in oxDNA_analysis_tools.  Make sure that you compiled oxDNA with `-DPython=1`, and then run
```
oat persistence_length trajectory.dat input_persistences 10 50 -p5
```

This script will calculate the persistence length, which should be around 120 nucleotides.

### Persistence length

Persistence length can be calculated from the average correlation between vectors tangent to a polymer.  In this case, we use vectors pointing from the midpoint of one base pair to the midpoint of the next base pair along the helix (note that this only works for DNA, not RNA). The persistence length can be found by performing an exponential fit to:
$$\langle \bf{n_x} ~ \dot ~ \bf{n_0} \rangle = B * exp(-x/L_{ps} )$$


Where $x$ is the offset and $\langle \bf{n_x} ~ \dot ~ \bf{n_0} \rangle$ is the correlation between reference vector $\bf{n_0}$ and a vector at offset $x$, $\bf{n_x}$
