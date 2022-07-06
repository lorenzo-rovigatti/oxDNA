#see: https://dna.physics.ox.ac.uk/index.php/Documentation#External_Forces


def mutual_trap(particle, ref_particle, stiff, r0, PBC):
    """
    A spring force that pulls a particle towards the position of another particle

    Parameters:
        particle (int): the particle that the force acts upon
        ref_particle (int): the particle that the particle will be pulled towards
        stiff (float): the force constant of the spring (in simulation units)
        r0 (float): the equlibrium distance of the spring
        PBC (bool): does the force calculation take PBC into account (almost always 1)
    """
    PBC = int(PBC)
    return({
        "type" : "mutual_trap",
        "particle" : particle,
        "ref_particle" : ref_particle,
        "stiff" : stiff, 
        "r0" : r0,
        "PBC" : PBC
    })


def string(particle, f0, rate, direction):
    """
    A linear force along a vector

    Parameters:
        particle (int): the particle that the force acts upon
        f0 (float): the initial strength of the force at t=0 (in simulation units)
        rate (float or SN string): growing rate of the force (simulation units/timestep)
        dir ([float, float, float]): the direction of the force
    """
    return({
        "type" : "string",
        "particle" : particle, 
        "f0" : f0, 
        "rate" : rate, 
        "dir" : direction 
    })


def harmonic_trap(particle, pos0, stiff, rate, direction):
    """
    A linear potential well that traps a particle

    Parameters:
        particle (int): the particle that the force acts upon
        pos0 ([float, float, float]): the position of the trap at t=0
        stiff (float): the stiffness of the trap (force = stiff * dx)
        rate (float): the velocity of the trap (simulation units/time step)
        direction ([float, float, float]): the direction of movement of the trap
    """
    return({
        "type" : "trap",
        "particle" : particle, 
        "pos0" : pos0,
        "rate" : rate,
        "dir" : direction
    })


def rotating_harmonic_trap(particle, stiff, rate, base, pos0, center, axis, mask):
    """
    A harmonic trap that rotates in space with constant angular velocity

    Parameters:
        particle (int): the particle that the force acts upon
        pos0 ([float, float, float]): the position of the trap at t=0
        stiff (float): the stiffness of the trap (force = stiff * dx)
        rate (float): the angular velocity of the trap (simulation units/time step)
        base (float): initial phase of the trap
        axis ([float, float, float]): the rotation axis of the trap
        mask([float, float, float]): the masking vector of the trap (force vector is element-wise multiplied by mask)
    """
    return({
        "type" : "twist", 
        "particle" : particle,
        "stiff" : stiff,
        "rate" : rate,
        "base" : base,
        "pos0" : pos0,
        "center" : center,
        "axis" : axis,
        "mask" : mask
    })


def repulsion_plane(particle, stiff, direction, position):
    """
    A plane that forces the affected particle to stay on one side.

    Parameters:
        particle (int): the particle that the force acts upon.  -1 will act on whole system.
        stiff (float): the stiffness of the trap (force = stiff * distance below plane)
        dir ([float, float, float]): the normal vecor to the plane
        position(float): position of the plane (plane is d0*x + d1*y + d2*z + position = 0)
    """
    return({
        "type" : "repulsion_plane",
        "particle" : particle,
        "stiff" : stiff,
        "dir" : direction,
        "position" : position
    })


def repulsion_sphere(particle, center, stiff, r0, rate=1):
    """
    A sphere that encloses the particle
    
    Parameters:
        particle (int): the particle that the force acts upon
        center ([float, float, float]): the center of the sphere
        stiff (float): stiffness of trap
        r0 (float): radius of sphere at t=0
        rate (float): the sphere's radius changes to r = r0 + rate*t
    """
    return({
        "type" : "sphere",
        "particle" : particle,
        "center" : center,
        "stiff" : stiff,
        "r0" : r0,
        "rate" : rate
    })
