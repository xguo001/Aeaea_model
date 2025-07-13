import numpy as np

def detect_photon_v2(photon_start, photon_end, cone_axis, alpha, R):
    """
    Parameters
    ----------
    photon_start, photon_end : (3,) array-like
        Cartesian coordinates of the segment's start and end.
    cone_axis : (3,) array-like
        Unit-length axis vector that defines the cone (pointing outward from origin).
    alpha : float
        Half-aperture of the cone in **radians**.
    R : float, optional
        Sphere radius.
    tolerance : float, optional
        Numerical tolerance.

    Returns
    -------
    bool
        True if the segment intersects the sphere at a point that lies inside
        the cone (i.e. within `alpha` of `cone_axis`); otherwise False.
    """
    photon_start, photon_end, cone_axis = map(np.asarray, (photon_start, photon_end, cone_axis))
    cone_axis = cone_axis / np.linalg.norm(cone_axis)  # ensure unit length

    d = photon_end - photon_start
    a = np.dot(d, d)
    b = 2.0 * np.dot(d, photon_start)
    c = np.dot(photon_start, photon_start) - R**2

    #need to catch 9 near zero to avoid error in the sqrt below
    if a < 1e-9:

        return False, None

    disc = b*b - 4*a*c
    if disc < 0:
        return False, None

    sqrt_disc = np.sqrt(disc)
    ts = [(-b - sqrt_disc)/(2*a), (-b + sqrt_disc)/(2*a)]

    cos_alpha = np.cos(alpha)

    for t in ts:
        if 0 <= t <= 1 :
            P = photon_start + t * d
            cos_theta = np.dot(P, cone_axis) / np.linalg.norm(P)  # |P| == R
            if cos_theta >= cos_alpha - 1e-9:
                return True, t

    return False, None

def detect_boundary(photon_end, r):
    #!!!!!!! <- boundary currently set to have the same radius as where detector is located

    photon_end = np.array(photon_end)
    distance = np.linalg.norm(photon_end)
    return distance > r

