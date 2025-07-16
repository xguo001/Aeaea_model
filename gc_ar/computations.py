import gc_ar.set_parameters as set_parameters
import numpy as np

# -----------------------------
# ENERGY DECAY
# -----------------------------
def energy_decay(energy,mu_a,r):
    energy= energy*np.exp(-mu_a*r)
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
    return position_1 + (position_1 - position_0)/2

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
