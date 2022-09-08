## Input file options

In the list below the expected type of the value is specified between < and >. A <bool> key can be specified with values 1|0, yes|no or true|false. The | symbol (pipe) is used to indicate the different values that can be used to specify a value for the key. Keys between [ and ] are optional.

Core options:

    T = <float>
        temperature of the simulation. It can be expressed in simulation units
        or kelvin (append a k or K after the value) or celsius (append a c or
        C after the value).
    [fix_diffusion = <bool>]
        if true, particles that leave the simulation box are brought back in
        via periodic boundary conditions. Defaults to true.
    [seed = <int>]
        seed for the random number generator. On Unix systems, defaults to
        either a number from /dev/urandom or to time(NULL)
    [confs_to_skip = <int>]
        how many configurations should be skipped before using the next one as
        the initial configuration, defaults to 0
    restart_step_counter = <boolean>/<bool>
        false means that the step counter will start from the value read in
        the configuration file, true means that the step counter will start
        from 0/if True oxDNA will reset the step counter to 0, otherwise it
        will start from the step counter found in the initial configuration.
        Defaults to False.
    [external_forces = <bool>]
        specifies whether there are external forces acting on the nucleotides
        or not. If it is set to 1, then a file which specifies the external
        forces' configuration has to be provided (see external_forces_file)
    [external_forces_file = <path>]
        specifies the file containing all the external forces' configurations.
        Currently there are six supported force types: string, twist, trap,
        repulsion_plane, repulsion_plane_moving and mutual_trap (see
        EXAMPLES/TRAPS for some examples)
    [back_in_box = <bool>]
        whether particles should be brought back into the box when a
        configuration is printed or not, defaults to false
    [lastconf_file = <path>]
        path to the file where the last configuration will be dumped
    trajectory_file = <path>
        path to the file which will contain the output trajectory of the
        simulation
    [binary_initial_conf = <bool>]
        whether the initial configuration is a binary configuration or not,
        defaults to false
    [lastconf_file_bin = <path>]
        path to the file where the last configuration will be printed in
        binary format, if not specified no binary configurations will be
        printed
    [print_reduced_conf_every = <int>]
        every how many time steps configurations containing only the centres
        of mass of the strands should be printed. If 0, no reduced
        configurations will be printed
    [reduced_conf_output_dir = <path>]
        path to the folder where reduced configurations will be printed
    [no_stdout_energy = <bool>]
        if true oxDNA will not print the default simulation output, including
        the energy, to stdout. Defaults to false
    [print_timings = <bool>]
        whether oxDNA should print out to a file performance timings at the
        end of the simulation or not, defaults to false
    [timings_filename = <path>]
        path to the file where timings will be printed
    [output_prefix = <string>]
        the name of all output files will be preceded by this prefix, defaults
        to an empty string
    [checkpoint_every = <int>]
        If > 0, it enables the production of checkpoints, which have a binary
        format. Beware that trajectories that do have this option enabled will
        differ from trajectories that do not. If this key is specified, at
        least one of checkpoint_file and checkpoint_trajectory needs to be
        specified
    [checkpoint_file = <string>]
        File name for the last checkpoint. If not specified, the last
        checkpoint will not be printed separately
    [checkpoint_trajectory = <string>]
        File name for the checkpoint trajectory. If not specified, only the
        last checkpoint will be printed
    [reload_from = <string>]
        checkpoint to reload from. This option is incompatible with the keys
        conf_file and seed, and requires restart_step_counter=0 as well as
        binary_initial_conf!=1
    [print_input = <bool>]
        make oxDNA write the input key=value pairs used by the simulation in a
        file named input.pid, with pid being the oxDNA pid. Defaults to False.
    conf_file = <string>
        path to the starting configuration
    steps = <int>
        length of the simulation, in time steps
    [equilibration_steps = <int>]
        number of equilibration steps. During equilibration, oxDNA does not
        generate any output. Defaults to 0
    time_scale = linear/log_lin
        a linear time_scale will make oxDNA print linearly-spaced
        configurations. a log_lin will make it print linearly-spaced cycles of
        logarithmically-spaced configurations.
    print_conf_interval = <int>
        if the time scale is linear, this is the number of time steps between
        the outputing of configurations, otherwise this is just the first
        point of the logarithmic part of the log_lin time scale
    print_conf_ppc = <int>
        mandatory only if time_scale == log_line. This is the number of
        printed configurations in a single logarithmic cycle.
    [print_energy_every = <int>]
        number of time steps between the outputing of the energy (and of the
        other default observables such as acceptance ratios in Monte Carlo
        simulations). Defaults to 0.
    verlet_skin = <float>
        width of the skin that controls the maximum displacement after which
        Verlet lists need to be updated.
    [list_type = verlet|cells|no]
        Type of neighbouring list to be used in CPU simulations. 'no' implies
        a O(N^2) computational complexity. Defaults to verlet.

-------------------------------------------------------------------------------

MD options:

    [reset_initial_com_momentum = <bool>]
        if true the momentum of the centre of mass of the initial
        configuration will be set to 0. Defaults to false to enforce the
        reproducibility of the trajectory
    [reset_com_momentum = <bool>]
        if true the momentum of the centre of mass will be set to 0 each time
        fix_diffusion is performed. Defaults to false to enforce the
        reproducibility of the trajectory
    [use_barostat = <bool>]
        apply an MC-like barostat to the simulation
    [P = <float>]
        the pressure of the simulation
    [delta_L = <float>]
        controls the box side change by the MC-like barostat
    backend = CPU
        For CPU FFS
    backend_precision = <any>
        CPU FFS may use any precision allowed for a normal CPU MD simulation
    sim_type = FFS_MD
        This must be set for an FFS simulation
    newtonian_steps = <int>
        number of integration timesteps after which the thermostat acts. Can
        be 1.
    pt = <float>
        probability of refreshing the momenta of each particle
    diff_coeff = <float>
        base diffusion coefficient. Either pt or diff_coeff should be
        specified in the input file
    gamma_trans = <float>
        translational damping coefficient for the Langevin thermostat. Either
        this or diff_coeff should be specified in the input file.
    bussi_tau = <int>
        correlation time, in time steps, for the stochastic evolution of the
        kinetic energy
    DPD_zeta = <float>
        translational damping coefficient for the DPD thermostat.
    [thermostat = no|refresh|brownian|langevin|srd]
        Select the simulation thermostat for MD simulations. 'no' means
        constant-energy simulations. 'refresh' is the Anderson thermostat.
        'brownian' is an Anderson-like thermostat that refreshes momenta of
        randomly chosen particles. 'langevin' implements a regular Langevin
        thermostat. 'srd' is an (experimental) implementation of a stochastic
        rotational dynamics algorithm. 'no' and 'brownian' are also available
        on CUDA. Defaults to 'no'.

-------------------------------------------------------------------------------

MC options:

    ensemble = nvt|npt
        ensemble of the simulation
    [check_energy_every = <float>]
        oxDNA will compute the energy from scratch, compare it with the
        current energy and throw an error if the difference is larger then
        check_energy_threshold. Defaults to 10.
    [check_energy_threshold = <float>]
        threshold for the energy check. Defaults to 1e-2f for single precision
        and 1e-6 for double precision.
    delta_translation = <float>
        controls the trial translational displacement, which is a randomly
        chosen number between -0.5*delta and 0.5*delta for each direction.
    delta_rotation = <float>
        controls the angular rotational displacement, given by a randomly
        chosen angle between -0.5*delta and 0.5*delta radians.
    delta_volume = <float>
        controls the volume change in npt simulations.
    P = <float>
        the pressure of the simulation. Used only if ensemble == npt.
    [adjust_moves = <bool>]
        if true, oxDNA will run for equilibration_steps time steps while
        changing the delta of the moves in order to have an optimal acceptance
        ratio. It does not make sense if equilibration_steps is 0 or not
        given. Defaults to false
    [maxclust = <int>]
        Default: N; maximum number of particles to be moved together. Defaults
        to the whole system
    [small_system = <bool>]
        Default: false; whether to use an interaction computation suited for
        small systems.
    [preserve_topology = <bool>]
        Default: false; sets a maximum size for the move attempt to 0.5, which
        guarantees that the topology of the system is conserved. Also prevents
        very large moves and might speed up simulations of larger systems,
        while suppressing diffusion
    [umbrella_sampling = <bool>]
        Default: false; whether to use umbrella sampling
    [op_file = <string>]
        Mandatory if umbrella_sampling is set to true; path to file with the
        description of the order parameter
    [weights_file = <string>]
        Mandatory if umbrella_sampling is set to true; path to file with the
        weights to use in umbrella sampling
    [last_hist_file = <string>]
        Optional if umbrella_sampling is set to true, otherwise ignored;
        Default: last_hist.dat; path to file where the histograms associated
        with umbrella sampling will be stored. This is printed with the same
        frequency as the energy file. Should become an observable sooner or
        later
    [traj_hist_file = <string>]
        Optional if umbrella_sampling is set to true, otherwise ignored;
        Default: traj_hist.dat; path to file where the series histograms
        associated with umbrella sampling will be stored, allowing to monitor
        the time evolution of the histogram and possibly to remove parts of
        the simulation. This is printed with the same frequency as the energy
        file. Should become an observable sooner or later
    [init_hist_file = <string>]
        Optional if umbrella_sampling is set to true, otherwise ignored;
        Default: none; path to a file to load a previous histogram from,
        useful if one wants to continue a simulation to obtain more
        statistics.
    [extrapolate_hist = <float>,<float>,..., <float>]
        Optional if umbrella_sampling is set to true, otherwise ignored;
        Default: none; series of temperatures to which to extrapolate the
        histograms. They can be given as float in reduced units, or the units
        can be specified as in the T option
    [safe_weights = <bool>]
        Default: true; whether to check consistency in between order parameter
        file and weight file. Only used if umbrella_sampling = true
    [default_weight = <float>]
        Default: none; mandatory if safe_weights = true; default weight for
        states that have no specified weight assigned from the weights file
    [skip_hist_zeros = <bool>]
        Default: false; Wether to skip zero entries in the traj_hist file
    [equilibration_steps = <int>]
        Default: 0; number of steps to ignore to allow for equilibration
    [type = rotation|traslation|possibly other as they get added]
        move to perform. No Defaults

-------------------------------------------------------------------------------

Interactions/DNA2Interaction.h options:

    salt_concentration = <float>
        sets the salt concentration in M
    [dh_lambda = <float>]
        the value that lambda, which is a function of temperature (T) and salt
        concentration (I), should take when T=300K and I=1M, defaults to the
        value from Debye-Huckel theory, 0.3616455
    [dh_strength = <float>]
        the value that scales the overall strength of the Debye-Huckel
        interaction, defaults to 0.0543
    [dh_half_charged_ends = <bool>]
        set to false for 2N charges for an N-base-pair duplex, defaults to 1

Interactions/DHSInteraction.h options:

    DHS_eps = <float>
        background dielectrci constant for reaction field treatment
    DHS_rcut = <float>
        cutoff for the reaction field treatment

Interactions/DNAInteraction.h options:

    [use_average_seq = <boolean>]
        defaults to yes
    [hb_multiplier = <float>]
        HB interaction multiplier applied to all the nucleotides having a
        custom numbered base whose magnitude is > 300, defaults to 1.0
    [max_backbone_force = <float>]
        defaults to nothing, has to be > 0) (if set to a float value, it
        specifies the maximum force (in reduced units) that the FENE bonds
        will exert. After the separation corresponding to the specified value,
        the potential has a form A x + B log(x), where A and B are computed
        automatically to obtain a continuous, differentiable and monotonically
        increasing potential. The computation involves the value of
        max_backbone_force as well as the value of max_backbone_force_far,
        below
    [max_backbone_force_far = <float>]
        defaults to 0.04, only read if max_backbone_force is set, has to be >
        0) (limit value of the force exerted by the bonded interactions for
        large separations, in reduced units. Used in the computation of A and
        B above. The default value is set to be very weak (0.04 = ~2pN), so
        that two neighbours will eventually get close without ever breaking
        base pairs

Interactions/KFInteraction.h options:

    KF_N = <int>
        number of patches
    KF_continuous = <bool>
        selects either the original KF model, valid only for MC simulations,
        or its continuous variant (see ACS Nano 10, 5459 (2016))
    KF_delta = <float>
        radial width of the patches
    KF_cosmax = <float>
        angular half-width of the patches
    [KF_N_B = <int>]
        number of patches on species B
    [KF_epsilon_AA = <float>]
        depth of the well of the patch-patch interaction between particles of
        species A
    [KF_epsilon_BB = <float>]
        depth of the well of the patch-patch interaction between particles of
        species B
    [KF_epsilon_AB = <float>]
        depth of the well of the patch-patch interaction between particles of
        unlike species
    [KF_sigma_AA = <float>]
        diameter controlling the repulsive interaction between particles of
        species A
    [KF_sigma_BB = <float>]
        diameter controlling the repulsive interaction between particles of
        species B
    [KF_sigma_AB = <float>]
        diameter controlling the repulsive interaction between particles of
        unlike species

Interactions/PatchyInteraction.h options:

    PATCHY_N = <int>
        number of patches
    [PATCHY_N_B = <int>]
        number of patches on species B
    [PATCHY_alpha = <float>]
        width of patches, defaults to 0.12
    [PATCHY_epsilon_AA = <float>]
        depth of the well of the patch-patch interaction between particles of
        species A
    [PATCHY_epsilon_BB = <float>]
        depth of the well of the patch-patch interaction between particles of
        species B
    [PATCHY_epsilon_AB = <float>]
        depth of the well of the patch-patch interaction between particles of
        different species
    [PATCHY_sigma_AA = <float>]
        diameter controlling the repulsive interaction between particles of
        species A
    [PATCHY_sigma_BB = <float>]
        diameter controlling the repulsive interaction between particles of
        species B
    [PATCHY_sigma_AB = <float>]
        diameter controlling the repulsive interaction between particles of
        different species

Interactions/BoxInteraction.h options:

    box_sides = <float>, <float>, <float>
        sides of the box

Interactions/TSPInteraction.h options:

    TSP_rfene = <float>
        FENE length constant for bonded interactions
    TSP_sigma[type] = <float>
        particle diameter associated to each interaction
    TSP_epsilon[type] = <float>
        energy scale associated to each interaction
    TSP_attractive[type] = <float>
        whether the interaction contains an attractive tail or not
    TSP_n[type] = <int>
        exponent for the generalised LJ potential for each interaction
    [TSP_attractive_anchor = <bool>]
        set to true if you want the anchor monomer to be of type B instead of
        type A. Defaults to false
    [TSP_only_chains = <bool>]
        if true the system will be composed of chains only. The topology will
        be interpreted accordingly by ignoring the first value of each line
        (which, in the case of TSPs, is the number of arms). Defaults to false
    [TSP_only_intra = <bool>]
        if true monomers belonging to different stars will not interact.
        Defaults to false

Interactions/HardSpheroCylinderInteraction.h options:

    length = <float>
        length of the spherocylinder

Interactions/InteractionFactory.h options:

    [interaction_type = DNA|RNA|HS|LJ|patchy|patchyDan|TSP|DNA_relax|DNA_nomesh|Box|HardCylinder|HardSpheroCylinder|DHS|Dirk]
        Particle-particle interaction of choice. Check the documentation
        relative to the specific interaction for more details. Defaults to
        dna.

Interactions/JordanInteraction.h options:

    JORDAN_N_patches = <int>
        number of patches
    [JORDAN_s = <float>]
        sigma of the gaussian modulation, defaults to 0.3
    [JORDAN_m = <float>]
        exponent to the 2m-m lennard-jones part
    [JORDAN_phi = <float>]
        angle below the equator for the rest position of the patches, defaults
        to PI/6
    [JORDAN_int_k = <float>]
        stiffness of the internal spring, defaults to 0., i.e., free patches

Interactions/RNAInteraction.h options:

    [use_average_seq = <boolean>]
        defaults to yes
    [seq_dep_file = <string>]
        sets the location of the files with sequence-dependent parameters
    [external_model = <string>]
        overrides default constants for the model, set in rna_model.h), by
        values specified by this option

Interactions/LJInteraction.h options:

    LJ_rcut = <float>
        interaction cutoff
    [LJ_kob_andersen = <bool>]
        Simulate a Kob-Andersen mixture. Defaults to false.
    [LJ_n = <int>]
        Generalised LJ exponent. Defaults to 6, which is the classic LJ value.

Interactions/RNAInteraction2.h options:

    [use_average_seq = <boolean>]
        defaults to yes
    [seq_dep_file = <string>]
        sets the location of the files with sequence-dependent parameters
    [external_model = <string>]
        overrides default constants for the model, set in rna_model.h), by
        values specified by this option
    [salt = <float>]
        sets the salt concentration in M, defaults to 1
    [mismatch_repulsion = <boolean>]
        defaults to no
    [mismatch_repulsion_strength = <float>]
        defaults to 1, sets the strength of repulsion if mismatch_repulsion is
        true

Interactions/HardCylinderInteraction.h options:

    height = <float>
        cylinder length

Interactions/TEPInteraction.h options:

    [_prefer_harmonic_over_fene = <bool>]
        if True, neighbouring beads are bound by an harmonic potential instead
        of a FENE one. Defaults to false.
    [_allow_broken_fene = <bool>]
        if True, the code won't die when the beads are beyond the acceptable
        FENE range. Defaults to True.

Interactions/DNAInteraction_relax.h options:

    relax_type = <string>
        Possible values: constant_force, harmonic_force; Relaxation algorithm
        used
    relax_strength = <float>
        Force constant for the replacement of the FENE potential

Interactions/DNAInteraction_relax2.h options:

    relax_type = <string>
        Possible values: constant_force, harmonic_force; Relaxation algorithm
        used
    relax_strength = <float>
        Force constant for the replacement of the FENE potential

Interactions/RNAInteraction_relax.h options:

    relax_type = <string>
        Possible values: constant_force, harmonic_force; Relaxation algorithm
        used
    relax_strength = <float>
        Force constant for the replacement of the FENE potential

-------------------------------------------------------------------------------

CUDA options:

    [CUDA_list = no|verlet]
        Neighbour lists for CUDA simulations. Defaults to 'no'.
    backend = CUDA
        For CUDA FFS -- NB unlike the CPU implementation, the CUDA
        implementation does not print extra columns with the current order
        parameter values whenever the energy is printed
    backend_precision = mixed
        CUDA FFS is currently only implemented for mixed precision
    sim_type = FFS_MD
        This must be set for an FFS simulation
    order_parameters_file = <string>
        path to the order parameters file
    ffs_file = <string>
        path to the file with the simulation stopping conditions. Optionally,
        one may use 'master conditions' (CUDA FFS only), which allow one to
        more easily handle very high dimensional order parameters. See the
        EXAMPLES/CUDA_FFS/README file for more information
    [ffs_generate_flux = <bool>]
        CUDA FFS only. Default: False; if False, the simulation will run until
        a stopping condition is reached; if True, a flux generation simulation
        will be run, in which case reaching a condition will cause a
        configuration to be saved but will not terminate the simulation. In
        the stopping condition file, the conditions must be labelled forward1,
        forward2, ... (for the forward conditions); and backward1, backward2,
        ... (for the backward conditions), ... instead of condition1,
        condition2, ... . To get standard flux generation, set the forward and
        backward conditions to correspond to crossing the same interface (and
        use conditions corresponding to different interfaces for Tom's flux
        generation). As with the single shooting run mode, the name of the
        condition crossed will be printed to stderr each time.
    [gen_flux_save_every = <integer>]
        CUDA FFS only. Mandatory if ffs_generate_flux is True; save a
        configuration for 1 in every N forward crossings
    [gen_flux_total_crossings = <integer>]
        CUDA FFS only. Mandatory if ffs_generate_flux is True; stop the
        simulation after N crossings achieved
    [gen_flux_conf_prefix = <string>]
        CUDA FFS only. Mandatory if ffs_generate_flux is True; the prefix used
        for the file names of configurations corresponding to the saved
        forward crossings. Counting starts at zero so the 3rd crossing
        configuration will be saved as MY_PREFIX_N2.dat
    [gen_flux_debug = <bool>]
        CUDA FFS only. Default: False; In a flux generation simulation, set to
        true to save backward-crossing configurations for debugging
    [check_initial_state = <bool>]
        CUDA FFS only. Default: False; in a flux generation simulation, set to
        true to turn on initial state checking. In this mode an initial
        configuration that crosses the forward conditions after only 1 step
        will cause the code to complain and exit. Useful for checking that a
        flux generation simulation does not start out of the A-state
    [die_on_unexpected_master = <bool>]
        CUDA FFS only. Default: False; in a flux generation simulation that
        uses master conditions, set to true to cause the simulation to die if
        any master conditions except master_forward1 or master_backward1 are
        reached. Useful for checking that a flux generation simulation does
        not enter any unwanted free energy basins (i.e. other than the initial
        state and the desired final state)
    [unexpected_master_prefix = <string>]
        CUDA FFS only. Mandatory if die_on_unexpected_master is True; the
        prefix used for the file names of configurations corresponding to
        reaching any unexpected master conditions (see
        die_on_unexpected_master).
    [CUDA_device = <int>]
        CUDA-enabled device to run the simulation on. If it is not specified
        or it is given a negative number, a suitable device will be
        automatically chosen.
    [CUDA_sort_every = <int>]
        sort particles according to a 3D Hilbert curve every CUDA_sort_every
        time steps. This will greatly enhnance performances for some types of
        interaction. Defaults to 0, which disables sorting.
    [threads_per_block = <int>]
        Number of threads per block on the CUDA grid. defaults to 2 * the size
        of a warp.

-------------------------------------------------------------------------------

Analysis options:

    [analysis_confs_to_skip = <int>]
        number of configurations that should be excluded from the analysis.
    analysis_data_output_<n> = {
    ObservableOutput
    }
        specify an analysis output stream. <n> is an integer number and should
        start from 1. The setup and usage of output streams are documented in
        the ObservableOutput class.

-------------------------------------------------------------------------------

Observables/HBEnergy.h options:

    [pairs_file = <string>]
        OrderParameter file containing the list of pairs whose HB energy is to
        be computed
    [base_file = <string>]
        file containing a list of nucleotides whose HB energy is to be
        computed, one nucleotide per line

Observables/ParticlePosition.h options:

    particle_id = <int>
        particle id
    [orientation = <bool>]
        defaults to false. If 1, it also prints out the orientation
    [absolute = <bool>]
        defaults to false. If 1, does not use periodic boundaries and it
        prints out the absolute position of the center of mass

Observables/StretchedBonds.h options:

    print_list = <bool>
        Whether to print the indexes of the particles that have stretched
        bonds. If set to false, only the total number of streched bonds is
        printed. Defaults to true
    [threshold = <float>]
        Threshold above which to report a stretched bond, in energy units.
        Default is 1.

Observables/VectorAngle.h options:

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

Observables/Step.h options:

    [units = steps|MD]
        units to print the time on. time in MD units = steps * dt, defaults to
        step

Observables/Pressure.h options:

    type = pressure
        an observable that computes the osmotic pressure of the system
    [stress_tensor = <bool>]
        if true, the output will contain 9 fields: the total pressure and the
        nine components of the stress tensor, xx, xy, xz, yx, yy, yz, zx, zy,
        zz

Observables/Pitch.h options:

    bp1a_id = <int>
        base pair 1 particle a id
    bp1b_id = <int>
        base pair 1 particle b id
    bp2a_id = <int>
        base pair 2 particle a id
    bp2b_id = <int>
        base pair 2 particle b id

Observables/SaltExtrapolation.h options:

    salts = <float>, <float>, ...
        list of salt concentration to extrapolate to
    temps = <T>, <T>, ...
        list of temperatures to extrapolate to, separated with commas.
        Temperatures can be specified in reduced units, Kelvin, Celsius as
        0.10105, 30C, 30c, 30 c, 303.15 k, 303.15K, 303.15k
    [op_file = <string>]
        order parameter file. If not found, it will use the one from the input
        file
    [weights_file = <string>]
        weights file. If not found, the one from the input file will be used.

Observables/Distance.h options:

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

Observables/TEPPlectonemePosition.h options:

    [bead_minimum_distance = <int>]
        the minimum integer separation between beads whose relative distance
        will be checked. Defaults to 7
    [distance_threshold = <float>]
        a plectoneme is identified when any two beads with indices farther
        away than bead_minimum_distance are closer than this distance.
        Defaults to 2

Observables/MeanVectorCosine.h options:

    chain_id = <int>
        chain id
    first_particle_position = <int>
        defaults to 0. position along the chain of  the first particle on
        which to compute the vector's cosine with the next particle
    last_particle_position = <int>
        defaults to N-2, where N is the number of elements of the chain.
        Position along the chain of the last particle over which to compute
        the vector's cosine with the next particle
    vector_to_average = <int>
        defaults to 1. Can be 1,2, or 3 depending on the vectors we wish to
        consider, or 0. In that case it measures the quantity (v2*v2')(v3*v3')
        - |v2 ^ v2||v3 ^ v3|

Observables/PotentialEnergy.h options:

    [split = <bool>]
        defaults to false, it tells the observable to print all the terms
        contributing to the potential energy

Observables/Writhe.h options:

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

Observables/CoaxVariables.h options:

    particle1_id = <int>
        particle 1 id
    particle2_id = <int>
        particle 2 id

Observables/ForceEnergy.h options:

    [print_group = <string>]
        limits the energy computation to the forces belonging to a specific
        group of forces. This can be set by adding a group_name option to each
        force's input. By default ForceEnergy computes the energy due to all
        the forces.

Observables/PairEnergy.h options:

    particle1_id = <int>
        particle 1 id
    particle2_id = <int>
        particle 2 id

Observables/StructureFactor.h options:

    max_q = <float>
        maximum q to consider
    [type = <int>]
        particle species to consider. Defaults to -1, which means "all
        particles"

Observables/DensityProfile.h options:

    max_value = <float>
        anything with a relevant coordinate grater than this will be ignored.
        Mind that the observable is PBC-aware.
    bin_size = <float>
        the bin size for the profile
    axis = <char>
        Possible values: x, y, z the axis along which to compute the profile

Observables/Contacts.h options:

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

Observables/PairForce.h options:

    [particle_id = <int>]
        Optional argument. particle id.

Observables/ObservableOutput.h options:

    name = <string>
        name of the output stream. stdout or stderr are accepted values
    print_every = <integer>
        frequency of output, in number steps for oxDNA, in number of
        configurations for DNAnalysis
    [start_from = <integer>]
        start outputing from the given step, defaults to 0
    [stop_at = <integer>]
        stop outputing at this step, defaults to -1 (which means never)
    [only_last = <bool>]
        if true, the output will not be appended to the stream, but it will
        overwrite the previous output each time, defaults to false
    [binary = <bool>]
        if true, the output will be printed in binary, defaults to false
    [linear = <bool>]
        if true the OutputObservable will save in linear scale, otherwise will
        use the logline scale by FS. Defaults to true
    [update_name_with_time = <bool>]
        if true the output filename will be changed by using the 'name' key as
        a prefix and the current step as a suffix. Defaults to false
    col_<n> = {
    type = name of the first observable
    [other observable options as lines of 'key = value']
    }
        this syntax specifies the column of the output file. Note that <n> is
        the column index and should start from 1

Observables/Rdf.h options:

    max_value = <float>
        maximum r to consider
    bin_size = <float>
        bin size for the g(r)
    [axes = <string>]
        Possible values: x, y, z, xy, yx, zy, yz, xz, zx. Those are the axes
        to consider in the computation. Mind that the normalization always
        assumes 3D sytems for the time being.

Observables/HBList.h options:

    type = hb_list
        name of  the observable
    only_count = <bool>
        if True, don't report the detailed binding profile but just count the
        bonds. Defaults to False.

Observables/Configurations/TEPtclOutput.h options:

    [back_in_box = <bool>]
        Default: true; if true the particle positions will be brought back in
        the box
    [show = <int>,<int>,...]
        Default: all particles; list of comma-separated indexes of the
        particles that will be shown. Other particles will not appear
    [hide = <int>,<int>,...]
        Default: no particles; list of comma-separated indexes of particles
        that will not be shown
    [print_labels = <bool>]
        Default: false; if true labels with the strand id are printed next to
        one end of the strand.
    [resolution = <int>]
        Default: 20; resolution set in the tcl file.
    [ref_particle = <int>]
        Default: -1, no action; The particle with the id specified (starting
        from 0) is set at the centre of the box. Overriden if ref_strands is
        specified. Ignored if negative or too large for the system.
    [ref_strand = <int>]
        Default: -1, no action; The strand with the id specified, starting
        from 1, is set at the centre of the box. Ignored if negative or too
        large for the system.

Observables/Configurations/TEPxyzOutput.h options:

    [back_in_box = <bool>]
        Default: true; if true the particle positions will be brought back in
        the box
    [show = <int>,<int>,...]
        Default: all particles; list of comma-separated indexes of the
        particles that will be shown. Other particles will not appear
    [hide = <int>,<int>,...]
        Default: no particles; list of comma-separated indexes of particles
        that will not be shown
    [ref_particle = <int>]
        Default: -1, no action; The particle with the id specified (starting
        from 0) is set at the centre of the box. Overriden if ref_strands is
        specified. Ignored if negative or too large for the system.
    [ref_strand = <int>]
        Default: -1, no action; The strand with the id specified, starting
        from 1, is set at the centre of the box. Ignored if negative or too
        large for the system.

Observables/Configurations/Configuration.h options:

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

Observables/Configurations/PdbOutput.h options:

    [back_in_box = <bool>]
        Default: true; if true the particle positions will be brought back in
        the box
    [show = <int>,<int>,...]
        Default: all particles; list of comma-separated indexes of the
        particles that will be shown. Other particles will not appear
    [hide = <int>,<int>,...]
        Default: no particles; list of comma-separated indexes of particles
        that will not be shown
    [ref_particle = <int>]
        Default: -1, no action; The nucleotide with the id specified (starting
        from 0) is set at the centre of the box. Overriden if ref_strands is
        specified. Ignored if negative or too large for the system.
    [ref_strand = <int>]
        Default: -1, no action; The strand with the id specified (starts from
        1) is set at the centre of the box. Ignored if negative or too large
        for the system.

Observables/Configurations/ChimeraOutput.h options:

    [colour_by_sequece = <bool>]
        Default: false; whether to coulour the bases according to the base
        type (A, C, G, T

Observables/Configurations/TclOutput.h options:

    [back_in_box = <bool>]
        Default: true; if true the particle positions will be brought back in
        the box
    [show = <int>,<int>,...]
        Default: all particles; list of comma-separated indexes of the
        particles that will be shown. Other particles will not appear
    [hide = <int>,<int>,...]
        Default: no particles; list of comma-separated indexes of particles
        that will not be shown
    [print_labels = <bool>]
        Default: false; if true labels with the strand id are printed next to
        one end of the strand.
    [resolution = <int>]
        Default: 20; resolution set in the tcl file.
    [ref_particle = <int>]
        Default: -1, no action; The nucleotide with the id specified (starting
        from 0) is set at the centre of the box. Overriden if ref_strands is
        specified. Ignored if negative or too large for the system.
    [ref_strand = <int>]
        Default: -1, no action; The strand with the id specified, starting
        from 1, is set at the centre of the box. Ignored if negative or too
        large for the system.

-------------------------------------------------------------------------------

Forces/COMForce.h options:

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

Forces/RepulsionPlane.h options:

    stiff = <float>
        stiffness of the repulsion.
    dir = <float>,<float>,<float>
        the vector normal to the plane: it should point towards the half-plane
        where the repulsion is not acting.
    position = <float>
        defines the position of the plane along the direction identified by
        the plane normal.
    particle = <int>
        index of the particle on which the force shall be applied. If -1, the
        force will be exerted on all the particles.

Forces/LJWall.h options:

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

Forces/GenericCentralForce.h options:

    particle = <int>
        comma-separated list of indices of particles to apply the force to. -1
        applies it to all particles. Entries separated by a dash "-" get
        expanded in a list of all the particles on a same strand comprised
        between the two indices. E.g., particle= 1,2,5-7 applies the force to
        1,2,5,6,7 if 5 and 7 are on the same strand.
    center = <float>,<float>,<float>
        the centre from which the force originates from

Forces/ConstantRateForce.h options:

    particle = <int>
        comma-separated list of indices of particles to apply the force to. -1
        applies it to all particles. Entries separated by a dash "-" get
        expanded in a list of all the particles on a same strand comprised
        between the two indices. E.g., particle= 1,2,5-7 applies the force to
        1,2,5,6,7 if 5 and 7 are on the same strand.
    F0 = <float>
        Initial force.
    rate = <float>
        growth rate of the force. It is [oxDNA energy units / (oxDNA distance
        units * (MD/MC) steps].
    [dir_as_centre = <bool>]
        if true the "dir" parameter will be interpreted as the origin of the
        force, so that the true direction will be dir - p->pos

Forces/RepulsiveSphere.h options:

    stiff = <float>
        stiffness of the repulsion.
    r0 = <float>
        radius of the sphere, in simulation units.
    rate = <float>
        rate of growth of the radius. Note that the growth is linear in
        timesteps/MC steps, not reduced time units.
    particle = <int>
        index of the particle on which the force shall be applied. If -1, the
        force will be exerted on all the particles
    [center = <float>,<float>,<float>]
        centre of the sphere, defaults to 0,0,0

Forces/SawtoothForce.h options:

    particle = <int>
        particle to apply the force to. -1 applies it to all particles.
    F0 = <float>
        Initial force
    wait_time = <float>
        time interval over which the force is constant. Units are (MD/MC)
        steps.
    increment = <float>
        amount by which to increment the force every wait_time steps.

Forces/RepulsionPlaneMoving.h options:

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

Forces/HardWall.h options:

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

-------------------------------------------------------------------------------

Forward Flux Sampling (FFS) options:

    backend = CPU/CUDA
        For CPU FFS/For CUDA FFS -- NB unlike the CPU implementation, the CUDA
        implementation does not print extra columns with the current order
        parameter values whenever the energy is printed
    backend_precision = <any>/mixed
        CPU FFS may use any precision allowed for a normal CPU MD
        simulation/CUDA FFS is currently only implemented for mixed precision
    sim_type = FFS_MD
        This must be set for an FFS simulation
    order_parameters_file = <string>
        path to the order parameters file
    ffs_file = <string>
        path to the file with the simulation stopping conditions. Optionally,
        one may use 'master conditions' (CUDA FFS only), which allow one to
        more easily handle very high dimensional order parameters. See the
        EXAMPLES/CUDA_FFS/README file for more information
    [ffs_generate_flux = <bool>]
        CUDA FFS only. Default: False; if False, the simulation will run until
        a stopping condition is reached; if True, a flux generation simulation
        will be run, in which case reaching a condition will cause a
        configuration to be saved but will not terminate the simulation. In
        the stopping condition file, the conditions must be labelled forward1,
        forward2, ... (for the forward conditions); and backward1, backward2,
        ... (for the backward conditions), ... instead of condition1,
        condition2, ... . To get standard flux generation, set the forward and
        backward conditions to correspond to crossing the same interface (and
        use conditions corresponding to different interfaces for Tom's flux
        generation). As with the single shooting run mode, the name of the
        condition crossed will be printed to stderr each time.
    [gen_flux_save_every = <integer>]
        CUDA FFS only. Mandatory if ffs_generate_flux is True; save a
        configuration for 1 in every N forward crossings
    [gen_flux_total_crossings = <integer>]
        CUDA FFS only. Mandatory if ffs_generate_flux is True; stop the
        simulation after N crossings achieved
    [gen_flux_conf_prefix = <string>]
        CUDA FFS only. Mandatory if ffs_generate_flux is True; the prefix used
        for the file names of configurations corresponding to the saved
        forward crossings. Counting starts at zero so the 3rd crossing
        configuration will be saved as MY_PREFIX_N2.dat
    [gen_flux_debug = <bool>]
        CUDA FFS only. Default: False; In a flux generation simulation, set to
        true to save backward-crossing configurations for debugging
    [check_initial_state = <bool>]
        CUDA FFS only. Default: False; in a flux generation simulation, set to
        true to turn on initial state checking. In this mode an initial
        configuration that crosses the forward conditions after only 1 step
        will cause the code to complain and exit. Useful for checking that a
        flux generation simulation does not start out of the A-state
    [die_on_unexpected_master = <bool>]
        CUDA FFS only. Default: False; in a flux generation simulation that
        uses master conditions, set to true to cause the simulation to die if
        any master conditions except master_forward1 or master_backward1 are
        reached. Useful for checking that a flux generation simulation does
        not enter any unwanted free energy basins (i.e. other than the initial
        state and the desired final state)
    [unexpected_master_prefix = <string>]
        CUDA FFS only. Mandatory if die_on_unexpected_master is True; the
        prefix used for the file names of configurations corresponding to
        reaching any unexpected master conditions (see
        die_on_unexpected_master).

-------------------------------------------------------------------------------
