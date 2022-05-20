# Observables

The observable infrastructure was devised to help building customized output from *oxDNA* (and *DNAnalysis*) without having to dive in the simulation code itself. 

The relevant keys in the input file are `analysis_data_output_*` (used by *DNAnalysis*) and `data_output_*` (used by *oxDNA*). These take as an argument, between curly brackets, a series of lines that is interpreted as an input file and whose options dictate the way the resulting output file is printed. 

````{warning}
`data_output_i` (with `i` integer) will not be taken into consideration if all `data_output_j` for `j < i` are present. The same holds true for `analysis_data_output_*`.
````

The following options are supported:

* `name`: the name of the output file. If `stdout` or `stderr` than the standard output or the standard error will be used, respectively. **Mandatory**.
* `linear`: if set to "true" then the output will be printed every `print_every` time steps. If set to "false" a log-linear scale will be used instead, where the output will be printed in linearly-spaced cycles made of log-spaced times. More specifically, the {math}`i`-th output time will be given by {math}`t_i = t_{jp} + n_0 (i - jp)^k`, where {math}`p` is the number of points per cycle, {math}`j = \lfloor i / p\rfloor` and {math}`n_0` and {math}`k` are additional parameters. Defaults to true.
* `print_every`: the number of time steps every which the output is written. **Mandatory if `linear = true`**.
* `log_ppc`: The number of points per cycle, {math}`p`. **Mandatory if `linear = false`**.
* `log_n0`: The parameter {math}`n_0`. **Mandatory if `linear = false`**.
* `log_fact`: The parameter {math}`k`. **Mandatory if `linear = false`**.
* `col_*`: an option that expects text enclosed by curly brackets that specifies the observable to be printed in the given column `*`. Note that `col_i` requires that all `col_j` with `j < i` are also present.
* `start_from`: the number of steps beyond which the output will be written. Defaults to 0.
* `stop_at`: the number of steps beyond which the output will stop being written. Defaults to -1 ("never").
* `only_last`: overwrite the content of the output with the current one.
* `update_name_with_time`: change the name of the output file every time the output is printed by appending the current simulation time step to `name`.

An example is:

```text
data_output_1 = {
	name = prova.dat
	print_every = 100
	col_1 = {
		type = step
		units = MD
	}
	col_2 = {
		type = potential_energy
	}
}
```

this will print in *prova.dat* two colums, the first with the steps in MD units (`dt`-aware) and the second with the potential energy. 


````{note}
oxDNA provides a plugin infrastructure to manage additional Observables. See the doxygen documentation for `src/PluginManagement/PluginManager.h`. The `contrib` folder contains a few working plugins you can use as starting points to write your owns.
````

## External observable file

Output files and observables can also be specified in an external [JSON](https://www.json.org/) file by setting the non-mandatory option `observables_file` in the [main input file](input.md#core-options). The following snippet shows an example of such a file:

```text
{
  "output_1" : {
    "print_every" : "10000",
    "name" : "hb_energy.dat",
    "cols" : [
      {
        "type" : "step",
        "units" : "MD"
      },
      {
        "type" : "hb_energy"
      }
    ]
  },
  "output_2" : {
    "print_every" : "1000",
    "name" : "pot_energy.dat",
    "cols" : [
      {
        "type" : "step"
      },
      {
        "type" : "potential_energy"
      }
    ]
  }
}
```

## `type = hb_energy`

    [pairs_file = <string>]
        OrderParameter file containing the list of pairs whose HB energy is to
        be computed
    [base_file = <string>]
        file containing a list of nucleotides whose HB energy is to be
        computed, one nucleotide per line

## `type = particle_position`

    particle_id = <int>
        particle id
    [orientation = <bool>]
        defaults to false. If 1, it also prints out the orientation
    [absolute = <bool>]
        defaults to false. If 1, does not use periodic boundaries and it
        prints out the absolute position of the center of mass

## `type = stretched`

    print_list = <bool>
        Whether to print the indexes of the particles that have stretched
        bonds. If set to false, only the total number of streched bonds is
        printed. Defaults to true
    [threshold = <float>]
        Threshold above which to report a stretched bond, in energy units.
        Default is 1.

## `type = vector_angle`

    [first_particle_index = <int>]
        defaults to 0. index of the first particle on which to compute the
        angle with the next particle.
    [last_particle_index = <int>]
        defaults to the index of the first-but last bead in the same strand as
        the first particle. Therefore, if I have a strand of N beads, the last
        one will be the one with index N-2. This is because the last bead is
        atypical in the TEP model (e.g. it's aligned with the vector before it
        rather than the one in front of it.). index of the last particle of
        the chain on which to compute the angle.
    [angle_index = <int>]
        defaults to 1. Can be 1,2, or 3 depending on which orientation vector
        we want to compute the cosine, or 0. In that case it measures the
        amount of twist, defined as in the TEP model: (1 +
        v2.v2'+v3.v3')/(2*pi)*(1 + v1.v1') where v1, v2 and v3 are the
        orientation vectors of the first particle in a pair and v1', v2' and
        v3' are the orientation vectors of the following particle.
    [print_local_details = <bool>]
        defaults to true. If true, print the quantity relative to each pair of
        particles. Otherwise, print their average (for angle_index == 1,2,3)
        OR their sum (that is, the total twist, for angle_index = 0.

## `type = step`

    [units = steps|MD]
        units to print the time on. time in MD units = steps * dt, defaults to
        step

## `type = pressure`

    type = pressure
        an observable that computes the osmotic pressure of the system
    [stress_tensor = <bool>]
        if true, the output will contain 9 fields: the total pressure and the
        nine components of the stress tensor, xx, xy, xz, yx, yy, yz, zx, zy,
        zz

## `type = pitch`

    bp1a_id = <int>
        base pair 1 particle a id
    bp1b_id = <int>
        base pair 1 particle b id
    bp2a_id = <int>
        base pair 2 particle a id
    bp2b_id = <int>
        base pair 2 particle b id

## `type = distance`

    particle_1 = <int>
        index of the first particle or comma-separated list of particle
        indexes composing the first set
    particle_2 = <int>
        index of the second particle or comma-separated list of the particle
        indexes composing the second set. The distance is returned as r(2) -
        r(1)
    [PBC = <bool>]
        Whether to honour PBC. Defaults to True
    [dir = <float>, <float>, <float>]
        vector to project the distance along. Beware that it gets normalized
        after reading. Defaults to (1, 1, 1) / sqrt(3)

## `type = potential_energy`

    [split = <bool>]
        defaults to false, it tells the observable to print all the terms
        contributing to the potential energy

## `type = writhe`

    [first_particle_index = <int>]
        defaults to 0. index of the first particle on which to compute the
        angle with the next particle.
    [last_particle_index = <int>]
        defaults to the index of the first-but last bead in the same strand as
        the first particle. Therefore, if I have a strand of N beads, the last
        one will be the one with index N-2. This is because the last bead is
        atypical in the TEP model (e.g. it's aligned with the vector before it
        rather than the one in front of it.). index of the last particle of
        the chain on which to compute the angle.
    [subdomain_size = <int>]
        if locate_plectonemes is false, defaults to the entire subchain,
        therefore computing the total writhe,otherwise it defaults to 35. If
        smaller, the writhe will be computed only on the beads between i and
        i+subdomain_size, for every i between first_particle_index and
        last_particle_index (wrapping around if go_round is true, or not
        computing it if the end of the chain is reached otherwise.
    [go_around = <bool>]
        whether to assume periodic boundary conditions when building
        subdomains - see above. Defaults to true if the last particle is right
        before the first and if subdomain_size is not the entire subchain, and
        to false otherwise.
    [locate_plectonemes = <bool>]
        if this is true, the writhe will be used to locate plectonemes with
        the algorithm from Vologodskii et al. "Conformational and
        THermodynamic Properties of Supercoiled DNA (1992)" and the indices on
        beads at the center of a plectoneme loop will be printed. Defaults to
        false.
    [writhe_threshold = <double>]
        if the writhe exceeds this, then we mark a plectoneme. Only used if
        locate_plectonemes is true. Defaults to 0.28, since this is the value
        that works best with a subdomain_size of 35, which is the default one.
    [print_space_position = <bool>]
        defaults to false. Whether to print the position of the plectoneme tip
        segment in 3d space as well as its index. Only used if
        locate_plectonemes = true.
    [print_size = <bool>]
        defaults to false. Whether to print the plectoneme size compute with
        Ferdinando-Lorenzo's algorithm. Only used if locate_plectonemes =
        true.
    [contact_threshold = <number>]
        defaults to 5. Segments closer than this will be considered to be
        touching accourding to the plectoneme size algorithm.
    [size_outer_threshold = <int>]
        defaults to 30. Outer threshold parameter, which substantially is the
        maximum separation in indices between two points of contact of a
        plectoneme loop, for the plectoneme size algorithm.
    [minimum_plectoneme_size = <int>]
        defaults to 1. Plectonemes shorter than this wont' be reported.
    [bending_angle_number_segments = <int>]
        defaults to 0. When non-zero, the angle between that many segments
        surrounding the plectoneme tip will be averaged and printed on file.

## `type = coax_variables`

    particle1_id = <int>
        particle 1 id
    particle2_id = <int>
        particle 2 id

## `type = force_energy`

    [print_group = <string>]
        limits the energy computation to the forces belonging to a specific
        group of forces. This can be set by adding a group_name option to each
        force's input. By default ForceEnergy computes the energy due to all
        the forces.

## `type = pair_energy`

    particle1_id = <int>
        particle 1 id
    particle2_id = <int>
        particle 2 id

## `type = structure_factor`

    max_q = <float>
        maximum q to consider
    [type = <int>]
        particle species to consider. Defaults to -1, which means "all
        particles"

## `type = density_profile`

    max_value = <float>
        anything with a relevant coordinate grater than this will be ignored.
        Mind that the observable is PBC-aware.
    bin_size = <float>
        the bin size for the profile
    axis = <char>
        Possible values: x, y, z the axis along which to compute the profile

## `type = contacts`

    [first_particle_index = <int>]
        defaults to 0. index of the first particle to consider. All the
        particles coming before this one will be ignored.
    [last_particle_index = <int>]
        defaults to the index of the first-but last bead in the same strand as
        the first particle. Therefore, if I have a strand of N beads, the last
        one will be the one with index N-2. This is because the last bead is
        atypical in the TEP model (e.g. it's aligned with the vector before it
        rather than the one in front of it.). index of the last particle to
        consider. All the particles coming before this one will be ignored.
    [neighbours_to_ignore = <int>]
        defalts to 1. Number of neighbours to ignore before-after each
        particle. E.g., if equals to 1, contacts between given first-
        neighbours will never be reported, if equals to 2, contacts between
        second neighbours will never be reported, etc.
    [contact_distance = <number>]
        defaults to 1. A contact is defined if the centers of mass of the
        particles is lower than this value.
    [only_outermost_contacts = <bool>]
        defaults to false. if true, contacts nested within other contacts will
        not be reported. E.g. if the i-th monomer is linked to both the i-1-th
        and the i+1-th monomer, and the contacts are 10-40, 10-25, 13-32,
        12-48 and 45-60, only 10-40, 12-48 and 45-60 will be reported, since
        10-25 and 13-32 are both nested inside 10-40. This is only get one
        result per plectoneme. Whatch out though, since this will report
        clashes between a plectoneme and the chain/other plectonemes. Telling
        a plectoneme and a plectoneme contact just by using the contact map
        might be non-banal.

## `type = pair_force`

    [particle_id = <int>]
        Optional argument. particle id.

## `type = rdf`

    max_value = <float>
        maximum r to consider
    bin_size = <float>
        bin size for the g(r)
    [axes = <string>]
        Possible values: x, y, z, xy, yx, zy, yz, xz, zx. Those are the axes
        to consider in the computation. Mind that the normalization always
        assumes 3D sytems for the time being.

## `type = hb_list`

    only_count = <bool>
        if True, don't report the detailed binding profile but just count the
        bonds. Defaults to False.

## `type = configuration`

    [back_in_box = <bool>]
        if true the particle positions will be brought back in the box,
        defaults to false
    [show = <int>,<int>,...]
        list of comma-separated particle indexes whose positions will be put
        into the final configuration
    [hide = <int>,<int>,...]
        list of comma-separated particle indexes whose positions won't be put
        into the final configuration
    [reduced = <bool>]
        if true only the strand centres of mass will be printed, defaults to
        false
