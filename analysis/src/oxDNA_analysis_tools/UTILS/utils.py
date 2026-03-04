import numpy as np

try:
    FLT_EPSILON = np.finfo(np.float).eps
except:
    FLT_EPSILON = 2.2204460492503131e-16

POS_BACK =  -0.4
POS_MM_BACK1 =  -0.3400
POS_MM_BACK2 =  0.3408
POS_STACK =  0.34
POS_BASE =  0.4
RNA_POS_BACK_a1 = -0.4
RNA_POS_BACK_a3 = 0.2
RNA_POS_BACK_a2 = 0.0

def get_pos_base(cm_pos, a1, a3 = None, type='DNA'):
        """
        Returns the position of the base-pairing site
        """
        return cm_pos + a1 * POS_BASE

def get_pos_stack(cm_pos, a1, a3=None, type='DNA'):
    return cm_pos + a1 * POS_STACK

def get_pos_back(cm_pos, a1, a3, type='DNA'):
    """
    Returns the position of the backbone centroid
    Note that cm_pos is the centrod of the backbone and base.
    """
    # In DNANucleotide.cpp, the definition of a1, a3, a2 is:
    # [1, 0, 0]
    # [0, 0, 1]
    # [0, 1, 0]
    # But because the frame is defined by a1, a3, to get a postive a2, that's a left-handed coordinate system.
    a2 = np.cross(a3, a1)
    if type=='DNA':
        return cm_pos + a1 * POS_MM_BACK1 + a2 * POS_MM_BACK2
    elif type=='RNA':
        return cm_pos + a1 * RNA_POS_BACK_a1 + a3 * RNA_POS_BACK_a3
    else:
        raise ValueError("Unknown type: {}".format(type))

def kabsch_align(arr:np.ndarray, ref:np.ndarray, center=True, inplace:bool=False, return_rot:bool=False):
    """
    Perform singular-value decomposition to align two 2D matricies

    Same procedure as align.svd_align, but not specialized to Configuration objects.

    Parameters:
        arr (np.ndarray) : The array to align
        ref (np.ndarray) : The array to align to
        center (bool) : Mean-center the arrays? If false, assume relative position is already correct.
        inplace (bool) : Modify the arrays in-place? If True, both arrays will be mean-centered, otherwise, copy both arrays and mean center arr on ref. (default False)
        return_rot (bool) : Return the rotation array? (default False)

    Returns:
        (np.ndarray|None|tuple[np.ndarray,np.ndarray]) : Returns None if `inplace==True` and `return_rot==False`, else returns the aligned array or (aligned, rot).
    """
    # Either operate on the arrays in-place or make copies
    if inplace:
        a = arr
        r = ref
    else:
        a = arr.copy()
        r = ref.copy()

    # mean-center both arrays
    if center:
        a_com = np.mean(a, axis=0)
        r_com = np.mean(r, axis=0)
        a -= a_com
        r -= r_com
        trans = r_com
    # Assume arrays are already "centered" on the desired point
    else:
        trans = np.zeros_like(a[0])

    # Compute the SVD of the covariance matrix and use that to get a rotation matrix
    cov = np.dot(a.T, r)
    u, _, vt = np.linalg.svd(cov)
    rot = np.dot(vt.T, u.T).T
    if np.linalg.det(rot) < 0:
        vt[2] = -vt[2]
        rot = np.dot(vt.T, u.T).T

    if inplace:
        np.dot(a, rot, out=a)
        if return_rot:
            return rot
        else:
            return None
    else:
        if return_rot:
            return (np.dot(a, rot) + trans, rot)
        else:
            return np.dot(a, rot) + trans

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
