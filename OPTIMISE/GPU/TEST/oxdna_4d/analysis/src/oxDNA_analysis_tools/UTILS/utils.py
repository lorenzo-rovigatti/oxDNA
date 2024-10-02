import numpy as np

try:
    FLT_EPSILON = np.finfo(np.float).eps
except:
    FLT_EPSILON = 2.2204460492503131e-16


def get_angle(a, b):
    """
    Get angle between a,b

    >>> a = [0, 1, 0]
    >>> b = [0, 0, 1]
    >>> round(get_angle(a,b),3)
    1.571

    """
    ab = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    if ab > (1. - FLT_EPSILON): 
        return 0
    elif ab < (-1. + FLT_EPSILON): 
        return np.pi
    else: 
        return np.arccos(ab)


def get_orthonormalized_base(v1, v2, v3):
    v1_norm2 = np.dot(v1, v1)
    v2_v1 = np.dot(v2, v1)

    v2 -= (v2_v1 / v1_norm2) * v1

    v3_v1 = np.dot(v3, v1)
    v3_v2 = np.dot(v3, v2)
    v2_norm2 = np.dot(v2, v2)

    v3 -= (v3_v1 / v1_norm2) * v1 + (v3_v2 / v2_norm2) * v2

    v1 /= np.sqrt(v1_norm2)
    v2 /= np.sqrt(v2_norm2)
    v3 /= np.sqrt(np.dot(v3, v3))

    return v1, v2, v3


def get_random_vector_in_sphere(r=1):
    r2 = r * r
    v = np.random.uniform(-r, r, 3)

    while np.dot(v, v) > r2:
        v = np.random.uniform(-r, r, 3)

    return v


def get_random_vector():
    ransq = 1.

    while ransq >= 1.:
        ran1 = 1. - 2. * np.random.random()
        ran2 = 1. - 2. * np.random.random()
        ransq = ran1 * ran1 + ran2 * ran2

    ranh = 2. * np.sqrt(1. - ransq)
    return np.array([ran1 * ranh, ran2 * ranh, 1. - 2. * ransq])


def get_random_rotation_matrix():
    v1, v2, v3 = get_orthonormalized_base(get_random_vector(), get_random_vector(), get_random_vector())

    R = np.array([v1, v2, v3])
    # rotations have det == 1
    if np.linalg.det(R) < 0: R = np.array([v2, v1, v3])

    return R


def get_rotation_matrix(axis, angle):
    ct = np.cos(angle)
    st = np.sin(angle)
    olc = 1. - ct
    x, y, z = axis / np.linalg.norm(axis)

    return np.array([[olc * x * x + ct, olc * x * y - st * z, olc * x * z + st * y],
                    [olc * x * y + st * z, olc * y * y + ct, olc * y * z - st * x],
                    [olc * x * z - st * y, olc * y * z + st * x, olc * z * z + ct]])
