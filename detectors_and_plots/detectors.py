import numpy as np
from photon_journey.computations import RFresnel, compute_ca1, compute_transmitted_direction, compute_reflected_direction,fresnel_mueller_matrices
from initialize.set_parameters import get_material
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

def photon_roulette(energy, chance):
    #implementing page 12 of the energy conservation paper
    #if W <= threshold
    #    if RND < chance, W = W / chance
    #    else W=0 (kill photon)
    #else do_nothing (photon continues)

    if np.random.rand() < chance:
        energy = energy / chance
    else:
        energy = 0

    return energy

def reflection_transmission(photon):
    if photon.died_detected==False:
        position= photon.position_hit_boundary
        direction = photon.direction
        radius = get_material('r')
        ca1,normal=compute_ca1(position,direction,radius)
        n=get_material('n')
        n1=get_material('n1')
        M_R, M_T, r, ca2 = fresnel_mueller_matrices(n,n1,ca1)
        r,ca2=RFresnel(n,n1,ca1)
        if r<=np.random.rand():
            direction_new= compute_reflected_direction(direction,normal,n,n1)
            photon.direction = direction_new
            photon.stokes = M_R @ photon.stokes
            photon.died_detected=True
            return photon
        else:
            direction_new = compute_transmitted_direction(direction, normal ,n, n1, ca2,ca1)
            photon.direction = direction_new
            photon.stokes = M_T @ photon.stokes
            return photon

    else:
        return photon


