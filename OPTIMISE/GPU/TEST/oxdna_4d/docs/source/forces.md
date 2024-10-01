# External forces

The code implements several types of external forces that can be imposed on the system that can be used either to simulate tension exerted on DNA or simply to accelerate the formation of secondary or tertiary structure. External forces can be tricky to treat, especially in a dynamics simulation, since they are an external source of work. Care should be taken in adjusting the time step, thermostat parameters and such.

To enable external forces, one needs to specify `external_forces = true` in the input file and also supply an external force file to read from with the key `external_forces_file = <file>`.

The syntax of the external forces file is quite simple. Examples of such files can be found in the hairpin formation (`examples/HAIRPIN`) and pseudoknot formation (`examples/PSEUDOKNOT`) examples. Each force is specified within a block contained in curly brackets. Empty lines and lines beginning with an hash symbol (#) are ignored. Different forces require different keys to be present. If the file has the wrong syntax, oxDNA should print out a sensible error message while parsing the file.

````{note}
All forces act on the centre of mass of the nucleotides. Currently, it is not possible to exert forces on specific nucleotide sites (*e.g.* base or backbone) or to exert torques.
````

Forces of different kinds can be combined in the same simulation, but there is a maximum number of 15 external forces per particle for memory reasons. This can be manually overridden recompiling the code with a different value of the macro `MAX_EXT_FORCES` in `src/defs.h`.

The main external forces implemented in the code, together with their associated options, are presented below. The most important ones are also accompanied by a brief description and an example.

````{note}
If `external_forces_as_JSON = true`, the external forces file will be parsed as [JSON](https://www.json.org/). This can be useful when using tools that automatise the generation of the force file. Here is a working example:

```
{
   "force_1":{
      "type":"trap",
      "particle":"0",
      "stiff":"1.00",
      "pos0":"46.3113780977, 11.1604626391, 26.8730311801",
      "rate":"0.0",
      "dir":"0.,0.,-1."
   },
   "force_2":{
      "type":"trap",
      "particle":"99",
      "pos0":"83.1532046751, 15.950789638, 37.3071701142",
      "stiff":"1.0",
      "rate":"0.0",
      "dir":"0.,0.,1."
   }
}
```

Note that **all fields** should be escaped with double quotes.
````

## Common options

The following two options are supported by all forces and can be used by observables or by the [oxpy](oxpy/index.md) library to retrieve information from or set options of specific external forces:

* `[group_name = <string>]`: an identifier that can be used to group forces together. Used, for instance, by the `force_energy` observable.
* `[id = <string>]`: a unique identifier that can be used to refer to this specific external force.

## String

A string is implemented as a force that does not depend on the particle position. Its value can be constant or can change linearly with time. It is useful as it does not fluctuate with time. A force of this kind is specified with `type = string`. The relevant keys are: 

* `particle = <int>`: comma-separated list of indices of particles to apply the force to. A value of `-1` or `all` applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. For instance, `particle = 1,2,5-7` applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.
* `F0 = <float>`: initial force.
* `rate = <float>`: growth rate of the force, in units of [oxDNA energy units / (oxDNA distance units * (MD/MC) steps].
* `dir = <float>,<float>,<float>`: direction of the force (automatically normalised by the code).
* `[dir_as_centre = <bool>]` if true the "dir" parameter will be interpreted as the origin of the force, so that the true direction will be `dir - p->pos`.

````{admonition} Example

The following bit of code will create an external force on the first nucleotide in the system starting at 1 simulation units (48.6 pN for DNA) and growing linearly with time at the rate of 1 simulation unit every {math}`10^6` time steps. The force will pull the nucleotide along the z direction.

	{
	type = string
	particle = 0
	F0 = 1.
	rate = 1e-6
	dir = 0., 0., 1.
	}

````

## Mutual Trap

This force is useful to form initial configurations. It is a harmonic force that at every moment pulls a particle towards a reference particle. It is possible to specify the separation at which the force will be 0.

````{warning}
Please note that the reference particle (`ref_particle` below) will not feel any force, thus making the name mutual trap somewhat misleading. If you want to have an actual *mutual* trap you will need to add a force on the reference particle.
````

A force of this kind is specified with `type = mutual_trap`. The relevant keys are: 

* `particle = <int>`: the particle on which to exert the force.
* `ref_particle = <int>`: particle to pull towards.
* `stiff = <float>`: stiffness of the trap.
* `r0 = <float`: equilibrium distance of the trap.

````{admonition} Example

Here is an example, extracted from the pseudoknot formation example (`examples/PSEUDOKNOT`). This will pull particle 14 towards particle 39, favouring an equilibrium distance of 1.4 (which corresponds roughly to the minimum of the hydrogen bonding potential, not a coincidence). The same force with opposite sign is exerted on particle 39 through a separate force. It is not necessary to have both particles feel the force, but it usually works much better. 

	{
	type = mutual_trap
	particle = 14
	ref_particle = 39
	stiff = 1.
	r0 = 1.2
	}
	
	{
	type = mutual_trap
	particle = 39
	ref_particle = 14
	stiff = 1.
	r0 = 1.2
	}

````

## Harmonic trap

This type of force implements an harmonic trap, of arbitrary stiffness, that can move linearly with time. It can be useful to fix the position of the nucleotides to simulate attachment to something or to implement (quasi) constant extension simulations.

A force of this kind is specified with `type = trap`. The relevant keys are:

* `particle = <int>`: comma-separated list of indices of particles to apply the force to. A value of `-1` or `all` applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. For instance, `particle = 1,2,5-7` applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.
* `pos0 = <float, float, float>`: rest position of the trap.
* `stiff = <float>`: stiffness of the trap (the force is stiff * dx).
* `rate = <float>`: speed of the trap (length simulation units/time steps).
* `dir = <float>,<float>,<float>`: direction of movement of the trap. Used only if rate is non zero.

````{admonition} Example

Here is an example input for a harmonic trap acting on the third nucleotide constraining it to stay close to the origin. In this example the trap does not move (rate=0), but one could have it move at a constant speed along the direction specified by dir, in this case the x direction. 

	{
	type = trap
	particle = 2
	pos0 = 0., 0., 0.
	stiff = 1.0
	rate = 0.
	dir = 1.,0.,0.
	}

````

## Rotating harmonic trap

Same as the harmonic trap, with the exception that the trap position rotates in space with constant angular velocity. Several of these can be used e.g. to twist DNA.

A force of this kind is specified with `type = twist`. The relevant keys are:

* `particle = <int>`: comma-separated list of indices of particles to apply the force to. A value of `-1` or `all` applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. For instance, `particle = 1,2,5-7` applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.
* `pos0 = <float>,<float>,<float>`: position of the trap when the rotation angle equals 0.
* `stiff = <float>`: stiffness of the trap (the force is stiff * dx).
* `rate = <float>`: angular velocity of the trap (length simulation units/time steps).
* `base = <float>`: initial phase of the trap.
* `axis = <float>,<float>,<float>`: rotation axis of the trap.
* `mask = <float>,<float>,<float>`: masking vector of the trap: the force vector will be element-wise multiplied by the masking vector.

````{admonition} Example

The following is an example input for a rotating trap acting on the first nucleotide forcing it to stay close to a point that starts at pos0 and then rotates around an axis containing the center point and parallel to the z axis. In this case we want to neglect the force component along the z-axis, so we set mask accordingly. 

	{
	type = twist
	particle = 0
	stiff = 1.00
	rate = 1e-5
	base = 0.
	pos0 = 15, 0.674909093169, 18.6187733563
	center = 13., 0.674909093169, 18.6187733563
	axis = 0, 0, 1
	mask = 1, 1, 0
	}

````

## Repulsion plane

This kind of external force implements a repulsion plane that constrains particles to stay on one side of it. It is implemented as a harmonic repulsion, but the stiffness can be made arbitrarily high to mimic a hard repulsion.

A force of this kind is specified with `type = repulsion_plane`. The relevant keys are:

* `particle = <int>`: comma-separated list of indices of particles to apply the force to. A value of `-1` or `all` applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. For instance, `particle = 1,2,5-7` applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.
* `stiff = <float>`: stiffness of the repulsion.
* `dir = <float>,<float>,<float>`: the vector normal to the plane: it should point towards the half-plane where the repulsion is not acting.
* `position = <float>`: defines the position of the plane along the direction identified by the plane normal.

If direction is `dir`{math}` = (u,v,w)`, then the plane contains all the points {math}`(x,y,z)` that satisfy the equation: {math}`u x + v y + w z + {\rm position} = 0`.
Only nucleotides  with coordinates {math}`(x,y,z)` that satisfy {math}`u x + v y + w z + {\rm position} \lt 0` will feel the force.
The force exerted on a nucleotide is equal to `stiff` \* *D*, where {math}`D = | ux + vy + wz + \mbox{position}| / \sqrt{u^2 + v^2 + w^2 }` is the distance of the nucleotide from the plane.
For nucleotides for which {math}`u x + v y + w z + {\rm position} \geq 0`, no force will be exerted.

````{admonition} Example

The following snippet defines a plane that acts on the whole system and will not exert any force on nucleotides with a positive x coordinate. A force proportional to 1 simulation unit \* x (48.6 pN \* x for DNA) will be exerted on all particles . 

	{
	type = twist
	particle = 0
	stiff = 1.00
	rate = 1e-5
	base = 0.
	pos0 = 15, 0.674909093169, 18.6187733563
	center = 13., 0.674909093169, 18.6187733563
	axis = 0, 0, 1
	mask = 1, 1, 0
	}

````

## Repulsive sphere

This force encloses particle(s) in a repulsive sphere. The repulsion force is harmonic.

A force of this kind is specified with `type = sphere`. The relevant keys are: 

* `particle = <int>`: comma-separated list of indices of particles to apply the force to. A value of `-1` or `all` applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. For instance, `particle = 1,2,5-7` applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.
* `stiff = <float>`: stiffness of the repulsion.
* `r0 = <float>`: radius of the sphere, in simulation units.
* `rate = <float>`: rate of growth of the radius. Note that the growth is linear in timesteps/MC steps, not reduced time units.
* `[center = <float>,<float>,<float>]`: centre of the sphere, defaults to `0,0,0`.

````{admonition} Example

The following snippet defines a repulsive sphere that acts on the first 50 nucleotides, confining them within a sphere of radius 6 centred in {math}`(10, 10, 10)`.

	{
	type = sphere
	particle = 0-50
	center = 10.0,10.0,10.0
	stiff = 10.0
	rate = 0.
	r0 = 6
	}

````

## `type = com`

    stiff = <float>
        stiffness of the spring
    r0 = <float>
        equilibrium elongation of the spring
    com_list = <string>
        comma-separated list containing the ids of all the particles whose
        centre of mass is subject to the force
    ref_list = <string>
        comma-separated list containing the ids of all the particles whose
        centre of mass is the reference point for the force acting on the
        other group of particles

## `type = LJ_wall`

    dir = <float>,<float>,<float>
        the vector normal to the plane: it should point towards the half-plane
        where the repulsion is not acting.
    position = <float>
        defines the position of the plane along the direction identified by
        the plane normal.
    particle = <int>
        index of the particle on which the force shall be applied. If -1, the
        force will be exerted on all the particles.
    [stiff = <float>]
        stiffness of the repulsion. Defaults to 1.
    [sigma = <float>]
        "Diameter" of the wall. It effectively rescales the distance between
        particle and wall. Defaults to 1.
    [n = <int>]
        Exponent of the 2n-n generalized Lennard-Jones expression. Defaults to
        6.
    [only_repulsive = <bool>]
        If true, the interactio between particle and wall gets cut-off at the
        minimum, resulting in a purely-repulsive potential. Defaults to false.
    [generate_inside = <bool>]
        If true the wall-particle interaction may not diverge, even for
        negative distances. Useful when generating the starting configuration.
        Defaults to false

## `type = sawtooth`

    particle = <int>
        particle to apply the force to. -1 applies it to all particles.
    F0 = <float>
        Initial force
    wait_time = <float>
        time interval over which the force is constant. Units are (MD/MC)
        steps.
    increment = <float>
        amount by which to increment the force every wait_time steps.

## `type = repulsion_plane_moving`

    stiff = <float>
        stiffness of the repulsion.
    dir = <float>,<float>,<float>
        the vector normal to the plane: it should point towards the half-plane
        where the repulsion is not acting.
    particle = <int>
        index(es) of the particle(s) on which the force shall be applied. Can
        be a list of comma-separated indexes. If -1, the force will be exerted
        on all the particles.
    ref_particle = <int>
        index(es) of the particle(s) whose position along the normal will be
        used to define the repulsive plane(s). Can be a list of comma-
        separated indexes.

## `type = hard_wall`

    dir = <float>,<float>,<float>
        the vector normal to the plane: it should point towards the half-plane
        that is allowed to the particles.
    position = <float>
        defines the position of the plane along the direction identified by
        the plane normal.
    particle = <int>
        index of the particle on which the force shall be applied. If -1, the
        force will be exerted on all the particles.
    [sigma = <float>]
        "Diameter" of the wall. It effectively rescales the distance between
        particle and wall. Defaults to 1.
        