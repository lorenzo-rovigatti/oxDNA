# Relaxing initial configurations

OxDNA structures generated with external tools such as [tacoxDNA](http://tacoxdna.sissa.it/), [oxview](https://sulcgroup.github.io/oxdna-viewer/) or [Adenita](https://doi.org/10.1093/nar/gkaa593) will generally contain regions that are locally stressed: over-bending or twisting, steric clashes, or excessively large separations of sites along the backbone. These local stresses could lead to numerical instabilities in an MD simulation and must be eliminated with preliminary relaxation runs. Note that some of these tools (notably tacoxDNA and oxview) can optionally provide a [force-file](forces.md) that can be used to keep nucleotide pairs bonded during the relaxation.

With oxDNA, a standard relaxation procedure starts with running short ($10^2$ to $10^4$ steps) Monte Carlo simulations on a CPU, followed by order of $10^6$ steps of an MD relaxation with a maximum-value of the cutoff for the backbone potential (ideally on CUDA, as simulations of large structures are prohibitively slow on CPU).

```{warning}
Many cadnano designs exported to oxDNA are near impossible to relax using the relaxation protocols described here. In these cases we advise to first use [oxview's rigid body simulator](https://github.com/sulcgroup/oxdna-viewer#rigid-body-simulations) and then proceed as described below.
```

## First stage: Monte Carlo relaxation

The aim of the first stage of the relaxation is to remove those issues that only require small displacements of the nucleotides (*e.g.* particle overlaps) and it is run on a single CPU core. As such, it is advisable to simulate a small number of time steps, especially for large configurations. The following snippet contains the options required to enable this type of relaxation:

```
sim_type = MC
ensemble = NVT

# these 2 options can be changed to improve performance, but it is usually not necessary
delta_translation = 0.1
delta_rotation = 0.1

steps = 1000 # you will very rarely need more than this number of steps

max_backbone_force = 5
max_backbone_force_far = 10

# set this to true to use a force-file (as given by some of the tacoxDNA/oxview tools)
external_forces = false
external_forces_file = YOURFILE.force
```

## Second stage: molecular dynamics relaxation

If only small local structural changes are required, this first stage may be sufficient to produce a configuration that is suitable as a starting point for an MD simulation. The aim of the second stage is to remove issues that require larger scale displacements of the nucleotides. In this stage, an MD simulation is run with the modified backbone potential; a typical duration might be $10^6$ MD steps,
but this will depend on the extent of the displacements required. To deal with the potentially large forces, it is also advantageous to increase the damping compared to a standard oxDNA MD simulation (by reducing the `diff_coeff` parameter). For large DNA structures, this more computationally demanding step is ideally performed on GPU-based architectures (CUDA).

The most important options that are linked to this type of relaxation are listed below:

```
sim_type = MD

dt = 0.002 # reduce this value to increase numerical stability
thermostat = bussi
bussi_tau = 1000
newtonian_steps = 53 # use 1 if you want the thermostat to act every step: 
				   # may be useful for very stressed configurations 

steps = 1e6 # the number of steps required may be 1 or 2 orders of magnitude 
            # larger than this figure, depending on the particular structure

max_backbone_force = 5
max_backbone_force_far = 10

# set this to true to use a force-file (as given by some of the tacoxDNA/oxview tools)
external_forces = false
external_forces_file = YOURFILE.force
```

