#see: https://dna.physics.ox.ac.uk/index.php/Documentation#External_Forces
from typing import Dict, List, Literal

def mutual_trap(particle:int, ref_particle:int, stiff:float, r0:float, PBC:bool, rate:float=0, stiff_rate:float=0) -> Dict:
    """
    A spring force that pulls a particle towards the position of another particle

    Parameters:
        particle (int): the particle that the force acts upon
        ref_particle (int): the particle that the particle will be pulled towards
        stiff (float): the force constant of the spring (in simulation units)
        r0 (float): the equlibrium distance of the spring
        PBC (0 or 1): does the force calculation take PBC into account (almost always 1)
        rate (float): changes r0 by this much every time step
        stiff_rate (float): changes stiff by this much every time step
    """
    return({
        "type" : "mutual_trap",
        "particle" : particle,
        "ref_particle" : ref_particle,
        "stiff" : stiff,
        "stiff_rate" : stiff_rate, 
        "r0" : r0,
        "rate" : rate,
        "PBC" : int(PBC)
    })


def string(particle:int, f0:float, rate:float, direction:List[float]) -> Dict:
    """
    A linear force along a vector

    Parameters:
        particle (int): the particle that the force acts upon
        f0 (float): the initial strength of the force at t=0 (in simulation units)
        rate (float): growing rate of the force (simulation units/timestep)
        direction ([float, float, float]): the direction of the force
    """
    return({
        "type" : "string",
        "particle" : particle, 
        "f0" : f0, 
        "rate" : rate, 
        "dir" : direction 
    })


def harmonic_trap(particle:int, pos0:List[float], stiff:float, rate:float, direction:List[float]) -> Dict:
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
        "stiff" : stiff,
        "rate" : rate,
        "dir" : direction
    })


def rotating_harmonic_trap(particle:int, pos0:List[float], stiff:float, rate:float, base:float, center:List[float], axis:List[float], mask:List[float]) -> Dict:
    """
    A harmonic trap that rotates in space with constant angular velocity

    Parameters:
        particle (int): the particle that the force acts upon
        pos0 ([float, float, float]): the position of the trap at t=0
        stiff (float): the stiffness of the trap (force = stiff * dx)
        rate (float): the angular velocity of the trap (simulation units/time step)
        base (float): initial phase of the trap
        center ([float, float, float]): the center of the circle
        axis ([float, float, float]): the rotation axis of the trap
        mask ([float, float, float]): the masking vector of the trap (force vector is element-wise multiplied by mask)
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


def repulsion_plane(particle:int, stiff:float, direction:List[float], position:List[float]) -> Dict:
    """
    A plane that forces the affected particle to stay on one side.

    Parameters:
        particle (int): the particle that the force acts upon.  -1 will act on whole system.
        stiff (float): the stiffness of the trap (force = stiff * distance below plane)
        direction ([float, float, float]): the normal vecor to the plane
        position ([float, float, float]): position of the plane (plane is d0*x + d1*y + d2*z + position = 0)
    """
    return({
        "type" : "repulsion_plane",
        "particle" : particle,
        "stiff" : stiff,
        "dir" : direction,
        "position" : position
    })


def repulsion_sphere(particle:int, center:List[float], stiff:float, r0:float, rate:float) -> Dict:
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
