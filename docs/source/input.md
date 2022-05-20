# Input file

The behaviour of oxDNA (and of the other executables) is controlled by a set of options that are usually listed in the input file, but can also be set from the [command line](usage.md).

## General syntax

* Options are set with the `key = value` syntax.
* Whitespace before `key`, around the equal sign and after `value` is never taken into account.
* Each line can contain one or more options separated by semicolons (`key1 = value1; key2 = value2`).
* Empty lines or lines starting with # are skipped.
* The value of a key can be referenced in another option with the `$(key)` syntax. For instance, if the input file contains `T = 300k`, then `energy_file = energy_T$(T).dat` is equivalent to `energy_file = energy_T300k.dat`. The order with which these options are set is not important.
* Mathematical expressions enclosed by `${...}` are evaluated and substituted into values. For instance, `steps = ${10 * exp(0)}` is equivalent to `steps = 10`.

````{note}
Mathematical expressions can contain references to other options. Here is an example:

```
steps = 1000
print_energy_every = ${$(steps) / 10}
```
````

The list of the main options supported by oxDNA is reported below. Square brackets mean optional, angular brackets specify the type of the value expected.

## Core options

These are the options that control the overall behaviour of the simulation and of the most common input/output operations.

* `T = <float>`: temperature of the simulation. It can be expressed in simulation units, kelvin (append a k or K after the value) or celsius (append a c or C after the value).
* `restart_step_counter = <bool>`: if `true` oxDNA will reset the step counter to 0, otherwise it will start from the step counter found in the initial configuration. Defaults to `false`.
* `steps = <int>`: length of the simulation, in time steps.
* `conf_file = <path>`: path to the starting configuration.
* `topology = <path>`: path to the file containing the system's topology.
* `trajectory_file = <path>`: path to the file which will contain the output trajectory of the simulation.
* `time_scale = linear/log_lin`: a linear time_scale will make oxDNA print linearly-spaced configurations. a log_lin will make it print linearly-spaced cycles of logarithmically-spaced configurations.
* `print_conf_interval = <int>`: if the time scale is linear, this is the number of time steps between the outputing of configurations, otherwise this is just the first point of the logarithmic part of the log_lin time scale.
* `print_energy_every = <int>`: number of time steps between the outputing of the energy (and of the other default observables such as acceptance ratios in Monte Carlo simulations).
* `[interaction_type = DNA|DNA2|RNA|RNA2|LJ|...]`: Particle-particle interaction of choice. Check the documentation relative to the specific interaction for more details and the supported interaction-specific options. Defaults to `DNA`.
* `[max_io = <float>]`: the maximum rate at which the output is printed, in MB/s. This is a useful option to avoid filling up the disk too quickly. Increase the default value (1 MB/s) at your own risk! 
* `[fix_diffusion = <bool>]`: if true, particles that leave the simulation box are brought back in via periodic boundary conditions. Defaults to `true`.
* `[fix_diffusion_every = <int>]`: number of time steps every which the diffusion is fixed. Used only if `fix_diffusion = true`, defaults to 100000 ({math}`10^5`).
* `[seed = <int>]`: seed for the random number generator. On Unix systems, defaults to either a number from /dev/urandom (if it exists and it's readable) or to time(NULL)
* `[confs_to_skip = <int>]`: how many configurations should be skipped before using the next one as the initial configuration, defaults to `0`.
* `[external_forces = <bool>]`: specifies whether there are external forces acting on the nucleotides or not. If it is set to `true`, then a file which specifies the external forces' configuration has to be provided (see below).
* `[external_forces_file = <path>]`:  specifies the file containing all the external forces' configurations. See [here](forces.md) for more details.
* `[back_in_box = <bool>]`: whether particles should be brought back into the box when a configuration is printed or not, defaults to false.
* `[lastconf_file = <path>]`: path to the file where the last configuration will be dumped. Defaults to `last_conf.dat`.
* `[binary_initial_conf = <bool>]`: whether the initial configuration is a binary configuration or not, defaults to `false`.
* `[lastconf_file_bin = <path>]`: path to the file where the last configuration will be printed in binary format. If not specified no binary configurations will be printed.
* `[print_reduced_conf_every = <int>]`: every how many time steps configurations containing only the centres of mass of the strands should be printed. If `0` (or not set), no reduced configurations will be printed.
* `[reduced_conf_output_dir = <path>]`: path to the folder where reduced configurations will be printed
* `[no_stdout_energy = <bool>]`: if `true `oxDNA will not print the default simulation output, including the energy, to stdout. Defaults to `false`.
* `[output_prefix = <string>]`: the name of all output files will be preceded by this prefix, defaults to an empty string.
* `[checkpoint_every = <int>]`: if > 0, it enables the production of checkpoints, which have a binary format. Beware that trajectories that do have this option enabled will differ from trajectories that do not. If this key is specified, at least one of `checkpoint_file` and `checkpoint_trajectory` needs to be specified.
* `[checkpoint_file = <string>]`: File name for the last checkpoint. If not specified, the last checkpoint will not be printed separately.
* `[checkpoint_trajectory = <string>]`: File name for the checkpoint trajectory. If not specified, only the last checkpoint will be printed.
* `[reload_from = <string>]`: checkpoint to reload from. This option is incompatible with the keys `conf_file` and `seed`, and requires `restart_step_counter = false` as well as `binary_initial_conf = true`. Note that this option is incompatible with `backend = CUDA`.
* `[print_input = <bool>]`: make oxDNA write the input key=value pairs used by the simulation in a file named input.pid, with pid being the oxDNA pid. Defaults to `false`.
* `[equilibration_steps = <int>]`: number of equilibration steps. During equilibration, oxDNA does not generate any output. Defaults to `0`.
* `[print_conf_ppc = <int>]`: this is the number of printed configurations in a single logarithmic cycle. Mandatory if `time_scale = log_line`.
* `[list_type = verlet|cells|no]`: type of neighbouring list to be used in CPU simulations. `no` implies a O(N^2) computational complexity. Defaults to `verlet`.
* `[verlet_skin = <float>]`: width of the skin that controls the maximum displacement after which Verlet lists need to be updated. mandatory if `list_type = verlet`.

## Molecular dynamics options

These options control the behaviour of MD simulations.

* `sim_type = MD|FFS_MD`: run either an MD or an FFS simulation.
* `backend = CPU|CUDA`: MD simulations can be run either on single CPU cores or on single CUDA-enabled GPUs.
* `backend_precision = <any>`: by default CPU simulations are run with `double` precision, CUDA with `mixed` precision (see [here](https://doi.org/10.1002/jcc.23763) for details). The CUDA backend also supports single precision (`backend_precision = float`), but we do not recommend to use it. Optionally, [by using CMake switches](install.md#CMake-options) it is possible to run CPU simulations in single precision or CUDA simulations in double precision.
* `[reset_initial_com_momentum = <bool>]`: if `true` the momentum of the centre of mass of the initial configuration will be set to 0. Defaults to `false` to enforce the reproducibility of the trajectory.
* `[reset_com_momentum = <bool>]`: if `true` the momentum of the centre of mass will be set to 0 each time fix_diffusion is performed. Defaults to `false` to enforce the reproducibility of the trajectory

### Constant-temperature simulations

* `[thermostat = no|refresh|brownian|langevin|DPD]`:  Select the simulation thermostat for MD simulations. `no` means constant-energy simulations.  `refresh` is the Anderson thermostat. `brownian` is an Anderson-like thermostat that refreshes momenta of randomly chosen particles. `langevin` implements a regular Langevin thermostat. `bussi` is the Bussi-Donadio-Parrinello thermostat. `DPD` is a Dissipative Particle Dynamics thermostat. The `no`, `brownian`, `langevin` and `bussi` thermostats are also available on CUDA. Defaults to `no`.
* `newtonian_steps = <int>`: number of integration timesteps after which the thermostat acts. Can be 1. Mandatory if there is a thermostat.
* `pt = <float>`: probability of refreshing the momenta of each particle. Used only if `thermostat = brownian`.
* `diff_coeff = <float>`:  base diffusion coefficient. Either this or `pt` should be specified if `thermostat = brownian`.
* `gamma_trans = <float>`: translational damping coefficient for the Langevin thermostat. Either this or `diff_coeff` should be specified in the input file if `thermostat = langevin`.
* `bussi_tau = <int>`: correlation time, in time steps, for the stochastic evolution of the kinetic energy for BDP thermostat. Mandatory if `thermostat = bussi`.
* `DPD_zeta = <float>`: translational damping coefficient for the DPD thermostat. Mandatory if `thermostat = DPD`.
* `DPD_rcut = <float>`: radial cut-off used by the DPD thermostat. Mandatory if `thermostat = DPD`.

### Constant-pressure simulations

* `[use_barostat = <bool>]`: apply an MC-like barostat to the simulation to keep the pressure constant. Defaults to `false`.
* `[P = <float>]`: the taget pressure of the simulation. Mandatory if `use_barostat = true`.
* `[delta_L = <float>]`: the extent of the box side change performed by the MC-like barostat. Mandatory if `use_barostat = true`.

## CUDA options

The following options require `backend = CUDA`.

* `[use_edge = <bool>]`: parallelise computations over interacting pairs rather than particles. It often results in a performance increase. Defaults to `false`.
* `[CUDA_list = no|verlet]`: neighbour lists for CUDA simulations. Defaults to `verlet`.
* `[cells_auto_optimisation = <bool>`: increase the size of the cells used to build Verlet lists if the total number of cells exceeds two times the number of nucleotides. Sometimes disabling this option increases performance. Used only if `CUDA_list = verlet`, defaults to `true`.
* `[max_density_multiplier = <float>]`: scale the size of data structures that store neighbours and cell lists. it is sometime necessary to increase this value (which also increases the memory footprint of the simulation) if the local density of nucleotides is high and the simulation crashes. defaults to `3`.
* `[print_problematic_ids = <bool>]`: if `true`, the code will print the indexes of particles that have very large coordinates (which may be caused by incorrectly-defined external forces and/or large time steps) before exiting. Useful for debugging purposes. Defaults to `false`.
* `[CUDA_device = <int>]`: CUDA-enabled device to run the simulation on. If it is not specified or it is given a negative number, a suitable device will be automatically chosen.
* `[CUDA_sort_every = <int>]`: sort particles according to a 3D Hilbert curve after the lists have been updated `CUDA_sort_every` times. This will greatly enhnance performances for some types of interaction. Defaults to `0`, which disables sorting.
* `[threads_per_block = <int>]`: Number of threads per block on the CUDA grid. defaults to 2 * the size of a warp.
* `[CUDA_avoid_cpu_calculations = <bool>]`: Do not run any computations on the CPU. If set to `true`, the energy will not be printed. It may speed up the simulation of very large systems. Defaults to `false`.
* `[CUDA_barostat_always_refresh = <bool>]`: Refresh the momenta of all particles after a successful volume move. Used only if `use_barostat = true`, defaults to `false`.
* `[CUDA_print_energy = <bool>]`: print the potential energy as computed on the GPU to the standard output. Useful for debugging purposes, since the "regular" potential energy is computed on the CPU.

## Monte Carlo options

The following options control the behaviour of MC simulations.

* `sim_type = MC|VMMC|MC2`: run regular (`MC`), Virtual Move (`VMMC`) or custom (`MC2`) Monte Carlo simulations.
* `ensemble = nvt|npt`:  ensemble of the simulation. It can be set to perform either constant temperature, volume and number of particles simulations (`nvt`) or constant-temperature, pressure and number of particles simulations (`npt`).
* `delta_translation = <float>`: set the maximum translational displacement, which is a randomly chosen number between -0.5\*delta and 0.5\*delta for each direction.
* `delta_rotation = <float>`: set the maximum angular rotational displacement, given by a randomly chosen angle between -0.5\*delta and 0.5\*delta radians.
* `delta_volume = <float>`: set the maximum change of the box edge in volume moves, given by a randomly chosen length between -0.5\*delta and 0.5\*delta radians. Mandatory if `ensemble = npt`.
* `P = <float>`: the target pressure of the simulation. Used only if `ensemble = npt`.
* `[check_energy_every = <int>]`: oxDNA will compute the energy from scratch, compare it with the current energy and throw an error if the difference is larger then `check_energy_threshold`. Defaults to `10`.
* `[check_energy_threshold = <float>]`: threshold for the energy check. Defaults to 0.01 for single precision and {math}`10^{-6}` for double precision.
* `[adjust_moves = <bool>]`: if `true`, oxDNA will run for `equilibration_steps` time steps while changing the delta of the moves in order to have an optimal acceptance ratio. It does not make sense if `equilibration_steps = 0` or not given. Defaults to `false`.
* `[maxclust = <int>]`: maximum number of particles to be moved together if `sim_type = VMMC`. Defaults to the size of the whole system.
* `[small_system = <bool>]`: whether to use an interaction computation suited for small systems. Defaults to `false`.
* `[preserve_topology = <bool>]`: set a maximum size for the move attempt to 0.5, which guarantees that the topology of the system is conserved. Also prevents very large moves and might speed up simulations of larger systems, while suppressing diffusion. Defaults to `false`.
* `[umbrella_sampling = <bool>]`: whether to use umbrella sampling. Defaults to `false`.
* `[op_file = <string>]`: path to file with the description of the order parameter. Mandatory if `umbrella_sampling = true`.
* `[weights_file = <string>]`: path to file with the weights to use in umbrella sampling. Mandatory if `umbrella_sampling = true`.
* `[last_hist_file = <string>]`: path to file where the histograms associated with umbrella sampling will be stored. This is printed with the same frequency as the energy file. Used if `umbrella_sampling = true`, defaults to `last_hist.dat`.
* `[traj_hist_file = <string>]`: path to file where the series histograms associated with umbrella sampling will be stored, allowing to monitor the time evolution of the histogram and possibly to remove parts of the simulation. This is printed with the same frequency as the energy file. Used if `umbrella_sampling = true`, defaults to `traj_hist.dat`.
* `[init_hist_file = <string>]`: path to a file to load a previous histogram from,
    useful if one wants to continue a simulation to obtain more
    statistics. Used if `umbrella_sampling = true`.
* `[extrapolate_hist = <float>,<float>,..., <float>]`:  series of temperatures to which to extrapolate the histograms. They can be given as float in reduced units, or the units can be specified as in the `T` option. Used if `umbrella_sampling = true`.
* `[safe_weights = <bool>]`: whether to check consistency in between order parameter
    file and weight file. Used if `umbrella_sampling = true`, defaults to `true`.
* `[default_weight = <float>]`: default weight for states that have no specified weight assigned from the weights file. Mandatory if `safe_weights = true`.
* `[skip_hist_zeros = <bool>]`: whether to skip zero entries in the `traj_hist` file. Mandatory if `umbrella_sampling = true`, defaults to `false`.

## Common options for `DNA`, `DNA2`, `RNA` and `RNA2` simulations

* `[use_average_seq = <bool>]`: use the average-sequence parameters. Defaults to `true`.
* `[seq_dep_file = <path>]`: path to the location of the file containing the sequence-dependent parameters. Mandatory if `use_average_seq = false`.
* `[max_backbone_force = <float>]`: specify the maximum force (in reduced units) that the FENE bonds will exert. After the separation corresponding to the specified value, the potential has a form {math}`A x + B \log(x)`, where *A* and *B* are computed automatically to obtain a continuous, differentiable and monotonically increasing potential. The computation involves the value of `max_backbone_force` as well as the value of `max_backbone_force_far`, below. It should be larger than 0.
* `[max_backbone_force_far = <float>]`: limit the value of the force exerted by the bonded interactions for large separations, in reduced units. Used in the computation of the *A* and *B* constants (see above). The default value is set to be very weak (0.04 = ~2pN), so that two neighbours will eventually get close without ever breaking base pairs. Only used if `max_backbone_force` is set, should be larger than 0.

## Common options for `DNA2` and `RNA2` simulations

* `[dh_lambda = <float>]`: the value that lambda, which is a function of temperature (T) and salt concentration (I), should take when T = 300K and I = 1M. Defaults to the value from Debye-Huckel theory, 0.3616455.
* `[dh_strength = <float>]`: the value that scales the overall strength of the Debye-Huckel interaction. Defaults to 0.0543.
* `[dh_half_charged_ends = <bool>]`: if `false`, nucleotides at the end of a strand carry a full charge, if `true` their charge is halved. Defaults to `true`.

## Common options for `DNA` and `DNA2` simulations

* `[hb_multiplier = <float>]`: hydrogen-bond interaction multiplier applied to all the nucleotides having a custom numbered base whose magnitude is > 300, defaults to 1.0.

## Options for `DNA2` simulations

* `salt_concentration = <float>`: the salt concentration in molar (M).

## Options for `RNA` and `RNA2` simulations

* `[external_model = <path>]`: override default constants for the model, set in `src/rna_model.h`, with the values specified in the file specified here.
        
## Options for `RNA2` simulations

* `[salt = <float>]`: the salt concentration in molar (M). Defaults to 1.
* `[mismatch_repulsion = <bool>]`: if `true`, add a repulsion between mismatches. Defaults to `false`.
* `[mismatch_repulsion_strength = <float>]`: set the strength of the repulsion between mismatches. Used only if `mismatch_repulsion = true`, defaults to `1`. 

## Options for `LJ` (Lennard-Jones) simulations

* `[LJ_rcut = <float>]`: interaction cut-off. Defaults to 2.5.
* `[LJ_kob_andersen = <bool>]`: set to `true` to simulate a Kob-Andersen mixture. Defaults to `false`.
* `[LJ_n = <int>]`: Generalised LJ exponent. Defaults to 6, which is the classic LJ value.

## Forward Flux Sampling (FFS) options

* `backend = CPU/CUDA`: FFS simulations can be run either on CPU or GPU. Note that, unlike the CPU implementation, the CUDA implementation does not print extra columns with the current order parameter values whenever the energy is printed.
* `backend_precision = <any>/mixed`: CPU FFS may use any precision allowed for a normal CPU MD simulation, while CUDA FFS is currently only implemented for mixed precision.
* `sim_type = FFS_MD`: This must be set for an FFS simulation.
* `order_parameters_file = <path>`: path to the order parameters file.
* `ffs_file = <string>`: path to the file with the simulation stopping conditions. Optionally, one may use `master conditions` (CUDA FFS only), which allow one to more easily handle very high dimensional order parameters. See the `EXAMPLES/CUDA_FFS/README` file for more information.
* `[ffs_generate_flux = <bool>]`: **CUDA FFS only**. If `false`, the simulation will run until a stopping condition is reached; if `true`, a flux generation simulation will be run, in which case reaching a condition will cause a configuration to be saved but will not terminate the simulation. In the stopping condition file, the conditions must be labelled forward1, forward2, ... (for the forward conditions); and backward1, backward2, ... (for the backward conditions), ... instead of condition1, condition2, ... . To get standard flux generation, set the forward and backward conditions to correspond to crossing the same interface. As with the single shooting run mode, the name of the condition crossed will be printed to stderr each time. Defaults to `false`.
* `[gen_flux_save_every = <integer>]`: **CUDA FFS only**. Save a configuration each time this number of forward crossings is achieved. Mandatory if `ffs_generate_flux = `true`.
* `[gen_flux_total_crossings = <integer>]`: **CUDA FFS only**. Stop the simulation after this number of crossings is achieved. Mandatory if `ffs_generate_flux = `true`.
* `[gen_flux_conf_prefix = <string>]`: **CUDA FFS only**. the prefix used for the file names of configurations corresponding to the saved forward crossings. Counting starts at zero so the 3rd crossing configuration will be saved as `MY_PREFIX_N2.dat`. Mandatory if `ffs_generate_flux = `true`.
* `[gen_flux_debug = <bool>]`: **CUDA FFS only**. In a flux generation simulation, set to `true` to save backward-crossing configurations for debugging. Defaults to `false`.
* `[check_initial_state = <bool>]`: **CUDA FFS only**. In a flux generation simulation, set to `true` to turn on initial state checking. In this mode an initial configuration that crosses the forward conditions after only 1 step will cause the code to complain and exit. Useful for checking that a flux generation simulation does not start out of the A-state. Defaults to `false`
* `[die_on_unexpected_master = <bool>]`: **CUDA FFS only**. In a flux generation simulation that uses master conditions, set to true to cause the simulation to die if any master conditions except master_forward1 or master_backward1 are reached. Useful for checking that a flux generation simulation does not enter any unwanted free energy basins (i.e. other than the initial state and the desired final state). Defaults to `false`.
* `[unexpected_master_prefix = <string>]`: **CUDA FFS only**. The prefix used for the file names of configurations corresponding to reaching any unexpected master conditions (see
    `die_on_unexpected_master`). Mandatory if `die_on_unexpected_master = true`.

## *DNAnalysis* options

* `[analysis_confs_to_skip = <int>]`: number of configurations that should be excluded from the analysis. Note that these configurations are parsed, initialised and then discarded. For big systems this is a slow process, in which case it is better to use the option below. Defaults to 0.
* `[analysis_bytes_to_skip = <int>]`: jump to this position in the trajectory file before starting the analysis. Useful to quickly analyse only portions of large trajectories. Defaults to 0.
* `[confs_to_analyse = <int>]`: the maximum number of configurations that should be analysed. if not set, the whole trajectory will be analysed.
 

## Observables

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

### `type = hb_energy`

    [pairs_file = <string>]
        OrderParameter file containing the list of pairs whose HB energy is to
        be computed
    [base_file = <string>]
        file containing a list of nucleotides whose HB energy is to be
        computed, one nucleotide per line

### `type = particle_position`

    particle_id = <int>
        particle id
    [orientation = <bool>]
        defaults to false. If 1, it also prints out the orientation
    [absolute = <bool>]
        defaults to false. If 1, does not use periodic boundaries and it
        prints out the absolute position of the center of mass

### `type = stretched`

    print_list = <bool>
        Whether to print the indexes of the particles that have stretched
        bonds. If set to false, only the total number of streched bonds is
        printed. Defaults to true
    [threshold = <float>]
        Threshold above which to report a stretched bond, in energy units.
        Default is 1.

### `type = vector_angle`

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

### `type = step`

    [units = steps|MD]
        units to print the time on. time in MD units = steps * dt, defaults to
        step

### `type = pressure`

    type = pressure
        an observable that computes the osmotic pressure of the system
    [stress_tensor = <bool>]
        if true, the output will contain 9 fields: the total pressure and the
        nine components of the stress tensor, xx, xy, xz, yx, yy, yz, zx, zy,
        zz

### `type = pitch`

    bp1a_id = <int>
        base pair 1 particle a id
    bp1b_id = <int>
        base pair 1 particle b id
    bp2a_id = <int>
        base pair 2 particle a id
    bp2b_id = <int>
        base pair 2 particle b id

### `type = distance`

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

### `type = potential_energy`

    [split = <bool>]
        defaults to false, it tells the observable to print all the terms
        contributing to the potential energy

### `type = writhe`

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

### `type = coax_variables`

    particle1_id = <int>
        particle 1 id
    particle2_id = <int>
        particle 2 id

### `type = force_energy`

    [print_group = <string>]
        limits the energy computation to the forces belonging to a specific
        group of forces. This can be set by adding a group_name option to each
        force's input. By default ForceEnergy computes the energy due to all
        the forces.

### `type = pair_energy`

    particle1_id = <int>
        particle 1 id
    particle2_id = <int>
        particle 2 id

### `type = structure_factor`

    max_q = <float>
        maximum q to consider
    [type = <int>]
        particle species to consider. Defaults to -1, which means "all
        particles"

### `type = density_profile`

    max_value = <float>
        anything with a relevant coordinate grater than this will be ignored.
        Mind that the observable is PBC-aware.
    bin_size = <float>
        the bin size for the profile
    axis = <char>
        Possible values: x, y, z the axis along which to compute the profile

### `type = contacts`

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

### `type = pair_force`

    [particle_id = <int>]
        Optional argument. particle id.

### `type = rdf`

    max_value = <float>
        maximum r to consider
    bin_size = <float>
        bin size for the g(r)
    [axes = <string>]
        Possible values: x, y, z, xy, yx, zy, yz, xz, zx. Those are the axes
        to consider in the computation. Mind that the normalization always
        assumes 3D sytems for the time being.

### `type = hb_list`

    only_count = <bool>
        if True, don't report the detailed binding profile but just count the
        bonds. Defaults to False.

### `type = configuration`

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
