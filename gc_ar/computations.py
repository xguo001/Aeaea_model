import gc_ar.set_parameters as set_parameters
import numpy as np

# -----------------------------
# ENERGY DECAY
# -----------------------------
def energy_decay(energy,mu_t,r):
    energy= energy*np.exp(-mu_t*r)
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
