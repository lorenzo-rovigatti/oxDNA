# External forces

## `type = string`

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

## `type = repulsion_plane`

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

## `type = repulsive_string`

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
        