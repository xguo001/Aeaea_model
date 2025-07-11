from computations import energy_decay, compute_phi, rotation_matrix_phi, rotation_matrix_glucose, mie_scattering_matrix_rayleigh
from set_parameters import initialize_photon, setup_detector, define_material
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


def change_direction(g):
    """ Sample scattering direction using Henyey-Greenstein phase function """
    cos_theta = (1 / (2 * g)) * (1 + g ** 2 - ((1 - g ** 2) / (1 - g + 2 * g * np.random.rand())) ** 2)
    sin_theta = np.sqrt(1 - cos_theta ** 2)
    phi = 2 * np.pi * np.random.rand()
    dx = sin_theta * np.cos(phi)
    dy = sin_theta * np.sin(phi)
    dz = cos_theta
    return np.array([dx, dy, dz])

# -----------------------------
# PHOTON PROPAGATION AND MATRIX OPERATIONS (Eq. 5)
# -----------------------------
def simulate_one_photon(GC):

    # Initialize variables
    material = define_material(GC)
    detector = setup_detector()
    pos, dir, stokes, energy = initialize_photon()
    mu_t = material["mu_s"] + material["mu_a"]
    total_path_length = 0
    alive = True
    detected = False
    step_counter = 0

    while alive:
        pos_start = pos.copy()

        # Travel step
        s = -np.log(np.random.rand()) / mu_t
        total_path_length += s
        pos += dir * s
        step_counter += 1

        # Energy decay
        energy = energy_decay(energy,mu_t, s)

        # Look to see if unalived
        if energy <= 1e-4:
            alive = False
            break

        # Look to see if detected
        detected_photon, t_value = detect_photon_v2(pos_start,pos,detector["cone_axis"],detector["alpha"],detector["r"])

        if detected_photon:
            #Recalculate path to only pre-detector path
            total_path_length -= (1-t_value)*s

        # Look to see if out of bound
        # !!!!!!! <- boundary currently set to have the same radius as where detector is located
        #boundary detection must come after photon detection because we have same radius for both.
        if not detected_photon and detect_boundary(pos,detector["r"]):
            alive = False
            break

        # Rotate Angle and Polarize due to Glucose
        theta_glucose = material["alpha"] * material["glucose_conc"] * s
        D = rotation_matrix_glucose(theta_glucose)
        stokes = D @ stokes

        #If detected, we record. If not detected, we scatter and travel again
        if detected_photon:

            return alive, step_counter, total_path_length, stokes

        # Scatter and Polarize due to scattering particle
        dir1= change_direction(g=0.9)
        phi= compute_phi(dir, dir1)
        dir=dir1
        M= mie_scattering_matrix_rayleigh(phi)
        D= rotation_matrix_phi(phi) @ M @ rotation_matrix_phi(phi)
        stokes = D @ stokes

    #If we reach here, alive = false, so the rest of return doesn't matter
    return alive, step_counter, total_path_length, stokes

