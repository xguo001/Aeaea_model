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
    # A detector is defined by an angle and the r
    return {
        "cone_axis": np.array([0.0, 0.0, 1]),
        "alpha": np.pi / 8,
        "r": 1
    }

def detect_photon(photon_start, photon_end, cone_axis, alpha, R):
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
        return False

    disc = b*b - 4*a*c
    if disc < 0:
        return False

    sqrt_disc = np.sqrt(disc)
    ts = [(-b - sqrt_disc)/(2*a), (-b + sqrt_disc)/(2*a)]

    cos_alpha = np.cos(alpha)

    for t in ts:
        if 0 <= t <= 1 :
            P = photon_start + t * d
            cos_theta = np.dot(P, cone_axis) / np.linalg.norm(P)  # |P| == R
            if cos_theta >= cos_alpha - 1e-9:
                return True

    return False

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
        "mu_s": 1.0,
        "mu_a": .10,
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

    # Initialize variables
    material = define_material(GC)
    detector = setup_detector()
    pos, dir, stokes, energy = initialize_photon()
    mu_t = material["mu_s"] + material["mu_a"]
    total_path_length = 0
    alive = True
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

        # Glucose rotation
        theta_glucose = material["alpha"] * material["glucose_conc"] * s
        D = rotation_matrix_glucose(theta_glucose)
        stokes = D @ stokes

        # Look to see if detected
        if detect_photon(pos_start,pos,detector["cone_axis"],detector["alpha"],detector["r"]):

            return alive, step_counter, total_path_length, stokes

        # Look to see if unalived
        if energy <= 1e-4:
            alive = False

        # Scatter if undetected
        dir1= change_direction(g=0.9)
        phi= compute_phi(dir, dir1)
        dir=dir1
        M= mie_scattering_matrix_rayleigh(phi)
        D= rotation_matrix_phi(phi) @ M @ rotation_matrix_phi(phi)
        stokes = D @ stokes

    return alive, step_counter, total_path_length, stokes

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
    step_counters = []
    death_counters = 0

    for _ in range(n_photons):
        alive, step_counter, total_path_length, stokes = simulate_one_photon(GC)

        if alive:
            results.append(rotation_angle_calculation(GC,total_path_length))
            step_counters.append(step_counter)
        else:
            death_counters += 1

    print ("This many died", death_counters)

    return np.mean(results), np.mean(step_counters)

# -----------------------------
# WRAPPER AROUND SINGLE GLUCOSE LEVEL SIMUATION WITH MULTITHREADING AND ANGLE ROTATION PRINTOUT
# -----------------------------
def simulate_one_gc(GC, n_photons):
    start = time.perf_counter()
    pid = os.getpid()
    print(f"[PID {pid}] processing GC = {GC} begins", flush=True)  # live announcement

    angle, steps = simulate_multiple_photon(GC, n_photons)

    elapsed = (time.perf_counter() - start)/60
    print(f"[PID {pid}] processing GC = {GC} ended in {elapsed:.2f} minutes", flush=True)  # live announcement
    print(f"GC = {GC}, Angle = {angle}, Average steps = {steps}")
    return angle

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":
   n_cores = 1
   n_photons = 1000
   GC_a=[2]
   A=[]

   with ProcessPoolExecutor(max_workers=n_cores) as pool:
       A = list(pool.map(simulate_one_gc, GC_a, repeat(n_photons)))

   plt.plot(GC_a,A)
   plt.xlabel("Glucose concentration")
   plt.ylabel("Rotation Angle")
   plt.title("GC vs. Angle")
   plt.show()
#    plt.savefig(r"/home/ubuntu/results/gc_angles.png",dpi=300)
#
#     photon_start = np.array([0, 0, 0])
#     photon_end = np.array([0, 0, 0.45])
#     cone_axis = np.array([0, 0, 1])
#     alpha = np.pi
#     R = 0.1
#
#     print("Detected:", detect_photon(photon_start, photon_end, cone_axis, alpha, R))