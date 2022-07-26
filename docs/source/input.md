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

The list of main options supported by oxDNA is reported below. Square brackets mean non-mandatory options, angular brackets specify the type of the value expected.

## Core options

These are the options that control the overall behaviour of the simulation and of the most common input/output operations.

* `T = <float>`: temperature of the simulation. It can be expressed in simulation units, kelvin (append a k or K after the value) or celsius (append a c or C after the value).
* `restart_step_counter = <bool>`: if `true` oxDNA will reset the step counter to 0, otherwise it will start from the step counter found in the initial configuration. Defaults to `false`.
* `steps = <int>`: length of the simulation, in time steps.
* `conf_file = <path>`: path to the starting configuration.
* `topology = <path>`: path to the file containing the system's topology.
* `trajectory_file = <path>`: path to the file which will contain the output trajectory of the simulation.
* `[trajectory_print_momenta = <bool>]`: print the linear and angular momenta of the particles to the trajectory. Set it to `false` to decrease the size of the trajectory by {math}`\approx 40\%`. Defaults to `true`.
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
* `[external_forces_file = <path>]`: specifies the file containing all the external forces' configurations. See [here](forces.md) for more details.
* `[external_forces_as_JSON = <bool>]`: if `true` the file specified by `external_forces_file` will be parsed as JSON. See [here](forces.md) for more details.
* `[observables_file = <path>]`: path to a file specifying additional outputs in JSON format. See [here](observables.md#external-observable-file) for more details.
* `[analysis_observables_file = <path>]`: as as above, but for `DNAnalysis`.
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
* `[print_conf_ppc = <int>]`: this is the number of printed configurations in a single logarithmic cycle. Mandatory if `time_scale = log_lin`.
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

## *confGenerator* options

* `[generate_consider_bonded_interactions = <bool>]`: if `true`, the generator will attempt to generate the position of a particle so that it is closer than `generate_bonded_cutoff` (see below) to its bonded neighbours. Defaults to `true`.
* `[generate_bonded_cutoff = <float>]`: the maximum distance at which the generator will put bonded neighbours. Defaults to `2.0`.
* `[energy_threshold = <float>]`: every time a particle is inserted in the system its total energy is computed, and if the resulting value is higher than this threshold than the insertion is cancelled and another trial position is generated. Increasing this value will make the generation of the initial configuration quicker, but its initial potential energy will be higher (on average). As a result, a more aggressive relaxation will be required.

## External forces

OxDNA supports several types of forces acting on and between specific nucleotides (see `external_forces_*` options above). See [here](forces.md) for additional details.

## Observables

In oxDNA the output can be customised with a set of *observables* that are described [here](observables.md).

## Plugins options

OxDNA provides a plugin infrastructure to make it easier to develop new observables, interactions and Monte Carlo moves. The following options are related to this functionality:

* `[plugin_search_path]`: colon-separated paths that will be searched when looking for plugins.
* `[plugin_do_cleanup]`: set to `false` to make it easier to use tools such as [GDB](https://www.sourceware.org/gdb/) or [valgrind](https://valgrind.org/) to debug the plugins.
* `[plugin_observable_entry_points]`: colon-separated names that will be used as entry points for observable plugins. Defaults to `make:make_observable`.
* `[plugin_interaction_entry_points]`: colon-separated names that will be used as entry points for interaction plugins. Defaults to `make:make_interaction`.
* `[plugin_move_entry_points]`: colon-separated names that will be used as entry points for Monte Carlo move plugins. Defaults to `make:make_move`.
