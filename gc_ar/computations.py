import gc_ar.set_parameters as set_parameters
import numpy as np


# -----------------------------
# ENERGY DECAY
# -----------------------------
def energy_decay(energy,mu_t,s):
    energy= energy*np.exp(-mu_t*s)
    return energy

def compute_phi(d_in, d_out):
    """
    Compute azimuthal rotation angle phi between incoming and outgoing photon directions. unit:radians
    """
    n_ref = np.array([0.0, 1.0, 0.0])  # y-axis as reference normal
    n_scat = np.cross(d_in, d_out)
    norm_ref = np.linalg.norm(n_ref)
    norm_scat = np.linalg.norm(n_scat)
    if norm_scat == 0:
        return 0.0
    dot = np.dot(n_ref, n_scat)
    cross = np.linalg.norm(np.cross(n_ref, n_scat))
    phi = np.arctan2(cross, dot)
    return phi

# -----------------------------
# ROTATION MATRIX AROUND Z (Eq. 3)
# -----------------------------
def rotation_matrix_phi(phi):
    return np.array([
        [1, 0, 0, 0],
        [0, np.cos(2*phi), np.sin(2*phi), 0],
        [0, -np.sin(2*phi), np.cos(2*phi), 0],
        [0, 0, 0, 1]
    ])

# -----------------------------
# ROTATION MATRIX FOR GLUCOSE EFFECT (Eq. 4)
# -----------------------------
def rotation_matrix_glucose(theta):
    thetar= np.radians(theta)
    return np.array([
        [1, 0, 0, 0],
        [0, np.cos(2*thetar), np.sin(2*thetar), 0],
        [0, -np.sin(2*thetar), np.cos(2*thetar), 0],
        [0, 0, 0, 1]
    ])

# -----------------------------
# MIE SCATTERING MATRIX (Eq. 1 with elements from Eq. 2)
# -----------------------------
def mie_scattering_matrix_rayleigh(theta_s):
    cos_theta = np.cos(theta_s)
    cos2_theta = cos_theta ** 2

    a = (3 / (16 * np.pi)) * (1 + cos2_theta)
    b = (3 / (16 * np.pi)) * (-1 + cos2_theta)
    d = (3 / (8 * np.pi)) * cos_theta
    e = 0.0

    return np.array([
        [a, b, 0, 0],
        [b, a, 0, 0],
        [0, 0, d, -e],
        [0, 0, e, d]
    ])

# -----------------------------
# TAKE THE TOTAL PATH LENGTH FROM A PHOTON'S JOURNEY AND RETURN THE ROTATION ANGLE (EQ.1 OF CHICKEN FINGER PAPER)
# -----------------------------
def rotation_angle_calculation(GC,total_path_length):
    angle=GC*(total_path_length/10)*set_parameters.get_material("alpha") #degree
    return angle

def mid_point(position_0,position_1):
    return position_0 + (position_0 - position_1)/2

def cut_path_at_boundary(photon_start, photon_end, R):
    photon_start, photon_end = map(np.asarray, (photon_start, photon_end))

    d = photon_end - photon_start
    a = np.dot(d, d)
    b = 2.0 * np.dot(d, photon_start)
    c = np.dot(photon_start, photon_start) - R ** 2

    if a < 1e-9:
        print ("warning: bug in cut_path_boundary()")
        return 0.0

    disc = b * b - 4 * a * c
    if disc < 0:
        print ("warning: bug in cut_path_boundary()")
        return 0.0  # no intersection

    sqrt_disc = np.sqrt(disc)
    ts = sorted([(-b - sqrt_disc) / (2 * a), (-b + sqrt_disc) / (2 * a)])

    # Clip to segment [0, 1]
    t1 = np.clip(ts[0], 0, 1)
    t2 = np.clip(ts[1], 0, 1)

    P1 = photon_start + t1 * d
    P2 = photon_start + t2 * d

    length_inside = np.linalg.norm(P2 - P1)
    total_length = np.linalg.norm(d)

    return length_inside / total_length


def change_direction():
    """ Sample scattering direction using Henyey-Greenstein phase function """
    g = set_parameters.get_material("g")
    cos_theta = (1 / (2 * g)) * (1 + g ** 2 - ((1 - g ** 2) / (1 - g + 2 * g * np.random.rand())) ** 2)
    sin_theta = np.sqrt(1 - cos_theta ** 2)
    phi = 2 * np.pi * np.random.rand()
    dx = sin_theta * np.cos(phi)
    dy = sin_theta * np.sin(phi)
    dz = cos_theta
    return np.array([dx, dy, dz])

def compute_ca1(position, direction, sphere_radius):
    """
    Compute ca1 for a photon near a spherical boundary centered at origin.

    Args:
        position: np.array([x, y, z]) – photon position
        direction: np.array([dx, dy, dz]) – photon direction (need not be normalized)
        sphere_radius: float – radius of spherical boundary

    Returns:
        ca1: cosine of angle of incidence (positive scalar)
    """
    normal = np.array(position) / sphere_radius


    theta = np.arctan2(normal[1], normal[0])
    #normal = position / sphere_radius  # since position is on the surface
    direction = direction / np.linalg.norm(direction)
    ca1 = abs(np.dot(direction, normal))
    return ca1, normal

def RFresnel(n1, n2, ca1):
    """
    Calculates Fresnel reflectance (r) and cosine of transmission angle (ca2)
    based on incident cosine ca1, and refractive indices n1 (current) and n2 (next).

    Returns:
        r   - Fresnel reflectance (probability of reflection)
        ca2 - Cosine of transmission angle (if transmitted)
    """

    sin_theta1 = np.sqrt(max(0.0, 1.0 - ca1 ** 2))
    print("n1", n1)
    print("n2", n2)
    sin_theta2 = (n1 / n2) * sin_theta1

    if sin_theta2 >= 1.0:
        # Total internal reflection
        return 1.0, None

    ca2 = np.sqrt(1.0 - sin_theta2 ** 2)

    # Fresnel equations (unpolarized)
    rs = ((n1 * ca1 - n2 * ca2) / (n1 * ca1 + n2 * ca2)) ** 2
    rp = ((n1 * ca2 - n2 * ca1) / (n1 * ca2 + n2 * ca1)) ** 2
    r = 0.5 * (rs + rp)

    return r, ca2


def compute_transmitted_direction(d_in, normal, n1, n2, ca2,ca1):
    """
    Computes transmitted direction vector given ca2 (cos theta2).

    Args:
        d_in: np.array([dx, dy, dz]) — incident direction (unit vector)
        normal: np.array([nx, ny, nz]) — surface normal (unit vector)
        n1: refractive index of current medium
        n2: refractive index of next medium
        ca2: cosine of transmission angle (theta2)

    Returns:
        d_out: transmitted direction vector (unit vector)
    """
    d_in = d_in / np.linalg.norm(d_in)
    normal = normal / np.linalg.norm(normal)

    eta = n1 / n2

    d_out = eta * d_in + (eta * ca1 - ca2) * normal
    d_out = d_out / np.linalg.norm(d_out)
    return d_out

def compute_reflected_direction(d_in, normal,n1,n2):
    """
    Computes the reflected direction vector when a photon hits a boundary.

    Args:
        d_in: np.array([dx, dy, dz]) — incident direction (unit vector)
        normal: np.array([nx, ny, nz]) — surface normal (unit vector)

    Returns:
        d_reflected: reflected direction vector (unit vector)
    """
    d_in = d_in / np.linalg.norm(d_in)
    normal = normal / np.linalg.norm(normal)
    ca1 = -np.dot(d_in, normal)  # cos(theta1)

    # Total internal reflection check
    sin_theta2_sq = (n1 / n2) ** 2 * (1.0 - ca1 ** 2)
    if sin_theta2_sq > 1.0:
        # TIR — must reflect
        d_reflected = d_in + 2 * ca1 * normal
        return d_reflected / np.linalg.norm(d_reflected)

    # Partial reflection — always reflect (ignore transmission)
    d_reflected = d_in + 2 * ca1 * normal
    return d_reflected / np.linalg.norm(d_reflected)

def fresnel_mueller_matrices(n1, n2, ca1):
    """
    Computes the Fresnel reflection and transmission Mueller matrices.

    Args:
        n1: refractive index of current medium
        n2: refractive index of next medium
        ca1: cosine of angle of incidence (theta1), scalar

    Returns:
        M_R: Mueller matrix for reflection (4x4)
        M_T: Mueller matrix for transmission (4x4)
        r: total unpolarized reflectance (scalar)
        ca2: cosine of transmitted angle (None if total internal reflection)
    """
    ca1 = float(np.clip(ca1, -1.0, 1.0))
    sin_theta1 = np.sqrt(max(0.0, 1.0 - ca1**2))
    sin_theta2 = (n1 / n2) * sin_theta1

    if sin_theta2 >= 1.0:
        # Total internal reflection
        return np.eye(4), np.zeros((4, 4)), 1.0, None

    ca2 = np.sqrt(1.0 - sin_theta2**2)

    # Fresnel coefficients for power reflectance
    rs = ((n1 * ca1 - n2 * ca2) / (n1 * ca1 + n2 * ca2))**2
    rp = ((n1 * ca2 - n2 * ca1) / (n1 * ca2 + n2 * ca1))**2
    r = 0.5 * (rs + rp)

    ts = 1 - rs
    tp = 1 - rp

    # Reflection Mueller matrix
    sqrt_rsrp = np.sqrt(rs * rp)
    M_R = 0.5 * np.array([
        [rs + rp, rs - rp,     0,          0],
        [rs - rp, rs + rp,     0,          0],
        [0,       0,       2*sqrt_rsrp,    0],
        [0,       0,           0,      2*sqrt_rsrp]
    ])

    # Transmission Mueller matrix
    sqrt_tstp = np.sqrt(ts * tp)
    M_T = 0.5 * np.array([
        [ts + tp, ts - tp,     0,          0],
        [ts - tp, ts + tp,     0,          0],
        [0,       0,       2*sqrt_tstp,    0],
        [0,       0,           0,      2*sqrt_tstp]
    ])

    return M_R, M_T, r, ca2

def mu_a_circular_dichroism(stokes):
    """
    Computes the effective absorption coefficient μ_a based on circular dichroism.

    Args:
        stokes: np.array([I, Q, U, V]) — Stokes vector
        mu_0: float — baseline absorption coefficient (isotropic component)
        delta_mu: float — differential absorption for circular polarization (CD strength)

    Returns:
        mu_a: float — effective absorption coefficient adjusted for polarization
    """
    mu_0=set_parameters.get_material("mu_a")
    delta_mu=set_parameters.get_material("dmu_a")
    I, Q, U, V = stokes

    if I <= 0:
        return mu_0  # avoid divide-by-zero or negative intensity

    circularity = V / I
    mu_at = mu_0 + delta_mu * circularity
    return mu_at