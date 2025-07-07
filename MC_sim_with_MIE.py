import numpy as np

# -----------------------------
# INITIALIZE PHOTON
# -----------------------------
def initialize_photon():
    position = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    stokes = np.array([1.0, 1.0, 0.0, 0.0])  # linearly polarized along x
    return position, direction, stokes

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
    return np.array([
        [1, 0, 0, 0],
        [0, np.cos(2*theta), np.sin(2*theta), 0],
        [0, -np.sin(2*theta), np.cos(2*theta), 0],
        [0, 0, 0, 1]
    ])

# -----------------------------
# MIE SCATTERING MATRIX (Eq. 1 with elements from Eq. 2)
# -----------------------------
def mie_scattering_matrix_rayleigh(theta):
    cos_theta = np.cos(theta)
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
# DETECTOR SETUP
# -----------------------------
def setup_detector():
    return {
        "position": np.array([0.0, 0.0, 0.12]),
        "acceptance_angle": np.pi / 8
    }

def detect_photon(photon_position, detector):
    dz = photon_position[2] - detector["position"][2]
    radial_distance = np.linalg.norm(photon_position[:2] - detector["position"][:2])
    theta = np.arctan2(radial_distance, abs(dz))
    return theta < detector["acceptance_angle"]

# -----------------------------
# MATERIAL PARAMETERS
# -----------------------------
def define_material():
    return {
        "mu_s": 100.0,
        "mu_a": 10.0,
        "g": 0.9,
        "n": 1.37,
        "alpha": 1.0e-5,
        "glucose_conc": 80.0,
        "thickness": 0.04,
    }

# -----------------------------
# PHOTON PROPAGATION AND MATRIX OPERATIONS (Eq. 5)
# -----------------------------
def simulate_one_photon():
    material = define_material()
    detector = setup_detector()

    pos, dir, stokes = initialize_photon()
    mu_t = material["mu_s"] + material["mu_a"]
    s = -np.log(np.random.rand()) / mu_t

    pos += dir * s

    # Step 1: Rotation due to glucose
    theta_glucose = material["alpha"] * material["glucose_conc"] * s
    D = rotation_matrix_glucose(theta_glucose)

    # Step 2: Set phi (scattering azimuth)
    phi = 0.0  # or randomly sampled
    R_phi = rotation_matrix_phi(phi)
    R_neg_phi = rotation_matrix_phi(-phi)

    # Step 3: Apply Eq. 5: R(-phi) * D * R(phi)
    M = R_neg_phi @ D @ R_phi
    stokes = M @ stokes

    # Step 4: Apply Rayleigh scattering matrix using Eq. 1 and 2
    scattering_angle = np.pi / 4  # example scattering angle (45 deg)
    M_scat = mie_scattering_matrix_rayleigh(scattering_angle)
    stokes = M_scat @ stokes

    if detect_photon(pos, detector):
        return stokes
    else:
        return None
