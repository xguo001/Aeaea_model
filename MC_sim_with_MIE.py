import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
import os,time
from itertools import repeat

# -----------------------------
# INITIALIZE PHOTON
# -----------------------------
def initialize_photon():
    position = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    stokes = np.array([1.0, 1.0, 0.0, 0.0])  # linearly polarized along x
    energy = 1
    return position, direction, stokes, energy



# -----------------------------
# ENERGY DECAY
# -----------------------------
def energy_decay(energy,mu_t,r):
    energy= energy*np.exp(-mu_t*r)
    return energy


def compute_phi(d_in, d_out):
    """
    Compute azimuthal rotation angle phi between incoming and outgoing photon directions.
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
    return np.array([
        [1, 0, 0, 0],
        [0, np.cos(2*theta), np.sin(2*theta), 0],
        [0, -np.sin(2*theta), np.cos(2*theta), 0],
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
# MATERIAL PARAMETERS
# -----------------------------
def define_material(GC):
    return {
        "mu_s": 100.0,
        "mu_a": 10.0,
        "g": 0.9,
        "n": 1.37,
        "alpha": 1.0e-5,
        "glucose_conc": 180.0,
        "thickness": 0.04,
    }

# -----------------------------
# PHOTON PROPAGATION AND MATRIX OPERATIONS (Eq. 5)
# -----------------------------
def simulate_one_photon(GC):
    material = define_material(GC)
    detector = setup_detector()
    pos, dir, stokes, energy = initialize_photon()
    mu_t = material["mu_s"] + material["mu_a"]
    total_path_length = 0
    alive=True
    while alive:
        # Travel step
        s = -np.log(np.random.rand()) / mu_t
        total_path_length += s
        pos += dir * s
        # Energy decay
        energy = energy_decay(energy,mu_t, s)
        # Glucose rotation
        theta_glucose = material["alpha"] * material["glucose_conc"] * s
        D = rotation_matrix_glucose(theta_glucose)
        stokes = D @ stokes
        if energy <= 1e-4:
            alive = False
        if detect_photon(pos, detector):
            return total_path_length, stokes
        dir1= change_direction(g=0.9)
        phi= compute_phi(dir, dir1)
        M= mie_scattering_matrix_rayleigh(phi)
        D= rotation_matrix_phi(phi) @ M @ rotation_matrix_phi(phi)
        stokes = D @ stokes







    # Step 2: Set phi
    phi =0.0  #randomly sampled
    R_phi = rotation_matrix_phi(phi)


    # Step 3: Apply Eq. 5: R(-phi) * D * R(phi)
    ##M = R_phi @ D @ R_phi
    ##stokes = M @ stokes

    # Step 4: Apply Rayleigh scattering matrix using Eq. 1 and 2
    scattering_angle = np.pi / 4  # example scattering angle (45 deg)
    M_scat = mie_scattering_matrix_rayleigh(scattering_angle)
    stokes = M_scat @ stokes

    if detect_photon(pos, detector):
        return stokes
    else:
        return None
    
# -----------------------------
# TAKE THE TOTAL PATH LENGTH FROM A PHOTON'S JOURNEY AND RETURN THE ROTATION ANGLE (EQ.1 OF CHICKEN FINGER PAPER)
# -----------------------------
def rotation_angle_calculation(GC,total_path_length):

    return GC*total_path_length*define_material(GC)["alpha"]

# -----------------------------
# TAKE A GLUCOSE LEVEL AND SIMULATE N NUMBER OF PHOTONS AND RETURN THE MEAN
# -----------------------------
def simulate_multiple_photon(GC, n_photons):
    results = []

    for _ in range(n_photons):
        stokes, total_path_length = simulate_one_photon(GC)
        results.append(rotation_angle_calculation(GC,total_path_length))

    return np.mean(results)

# -----------------------------
# WRAPPER AROUND SINGLE GLUCOSE LEVEL SIMUATION WITH MULTITHREADING AND ANGLE ROTATION PRINTOUT
# -----------------------------
def simulate_one_gc(GC, n_photons):
    start = time.perf_counter()
    pid = os.getpid()
    print(f"[PID {pid}] processing GC = {GC} begins", flush=True)  # live announcement

    angle = simulate_multiple_photon(GC, n_photons)

    elapsed = (time.perf_counter() - start)/60
    print(f"[PID {pid}] processing GC = {GC} ended in {elapsed:.2f} minutes", flush=True)  # live announcement
    print(f"GC = {GC}, Angle = {angle:.2f}")
    return angle

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":
    n_cores = 1
    n_photons = 1
    GC_a=[2,6,10,18,26]
    A=[]

    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        A = list(pool.map(simulate_one_gc, GC_a, repeat(n_photons)))

    plt.plot(GC_a,A)
    plt.xlabel("Glucose concentration")
    plt.ylabel("Rotation Angle")
    plt.title("GC vs. Angle")
    plt.show()
#    plt.savefig(r"/home/ubuntu/results/gc_angles.png",dpi=300)
