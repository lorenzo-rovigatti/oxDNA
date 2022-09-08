try:
    import numpy as np
except:
    import mynumpy as np

from base import FLT_EPSILON

def get_angle(a, b):
    """
    Get angle between a,b

    >>> a = [0, 1, 0]
    >>> b = [0, 0, 1]
    >>> round(get_angle(a,b),3)
    1.571

    """
    ab = np.dot(a, b)
    if ab > (1.-FLT_EPSILON): return 0
    elif ab < (-1.+FLT_EPSILON): return np.pi
    else: return np.arccos(ab)

def get_orthonormalized_base(v1, v2, v3):
    v1_norm2 = np.dot(v1, v1)
    v2_v1 = np.dot(v2, v1)

    v2 -= (v2_v1/v1_norm2) * v1

    v3_v1 = np.dot(v3, v1)
    v3_v2 = np.dot(v3, v2)
    v2_norm2 = np.dot(v2, v2)

    v3 -= (v3_v1/v1_norm2) * v1 + (v3_v2/v2_norm2) * v2

    v1 /= np.sqrt(v1_norm2)
    v2 /= np.sqrt(v2_norm2)
    v3 /= np.sqrt(np.dot(v3, v3))

    return v1, v2, v3

def get_random_vector_in_sphere(r=1):
    r2 = r*r
    v = np.random.uniform(-r, r, 3)

    while np.dot(v, v) > r2:
        v = np.random.uniform(-r, r, 3)

    return v

def get_random_vector():
    ransq = 1.

    while ransq >= 1.:
        ran1 = 1. - 2. * np.random.random()
        ran2 = 1. - 2. * np.random.random()
        ransq = ran1*ran1 + ran2*ran2

    ranh = 2. * np.sqrt(1. - ransq)
    return np.array([ran1*ranh, ran2*ranh, 1. - 2. * ransq])

def get_random_rotation_matrix():
    v1, v2, v3 = get_orthonormalized_base(get_random_vector(), get_random_vector(), get_random_vector())

    R = np.array([v1, v2, v3])
    # rotations have det == 1
    if np.linalg.det(R) < 0: R = np.array([v2, v1, v3])

    return R

def get_rotation_matrix(axis, anglest):
    """
    The argument anglest can be either an angle in radiants
    (accepted types are float, int or np.float64 or np.float64)
    or a tuple [angle, units] where angle a number and
    units is a string. It tells the routine whether to use degrees,
    radiants (the default) or base pairs turns

    axis --- Which axis to rotate about
        Ex: [0,0,1]
    anglest -- rotation in radians OR [angle, units]
        Accepted Units:
            "bp"
            "degrees"
            "radiants"
        Ex: [np.pi/2] == [np.pi/2, "radians"]
        Ex: [1, "bp"]

    """
    if not isinstance (anglest, (np.float64, np.float32, float, int)):
        if len(anglest) > 1:
            if anglest[1] in ["degrees", "deg", "o"]:
                angle = (np.pi / 180.) * anglest[0]
                #angle = np.deg2rad (anglest[0])
            elif anglest[1] in ["bp"]:
								# Notice that the choice of 35.9 DOES NOT correspond to the minimum free energy configuration.
								# This is usually not a problem, since the structure will istantly relax during simulation, but it can be
								# if you need a configuration with an equilibrium value for the linking number.
								# The minimum free energy angle depends on a number of factors (salt, temperature, length, and possibly more),
								# so if you need a configuration with 0 twist make sure to carefully choose a value for this angle
								# and force it in some way (e.g. by changing the angle value below to something else in your local copy).
                # Allow partial bp turns
                angle = float(anglest[0]) * (np.pi / 180.) * 35.9
                # Older versions of numpy don't implement deg2rad()
                #angle = int(anglest[0]) * np.deg2rad(35.9)
            else:
                angle = float(anglest[0])
        else:
            angle = float(anglest[0])
    else:
        angle = float(anglest) # in degrees, I think

    axis = np.array(axis)
    axis /= np.sqrt(np.dot(axis, axis))

    ct = np.cos(angle)
    st = np.sin(angle)
    olc = 1. - ct
    x, y, z = axis

    return np.array([[olc*x*x+ct, olc*x*y-st*z, olc*x*z+st*y],
                    [olc*x*y+st*z, olc*y*y+ct, olc*y*z-st*x],
                    [olc*x*z-st*y, olc*y*z+st*x, olc*z*z+ct]])
