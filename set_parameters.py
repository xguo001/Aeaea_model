import numpy as np

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
# DETECTOR SETUP
# -----------------------------
def setup_detector():
    # A detector is defined by an angle and the r
    return {
        "cone_axis": np.array([0.0, 0.0, 1]),
        "alpha": np.pi / 8,
        "r": 2
    }

# -----------------------------
# MATERIAL PARAMETERS
# -----------------------------
def define_material(GC):
    return {
        "mu_s": 1,
        "mu_a": .10,
        "g": 0.9,
        "n": 1.37,
        "alpha": 1.0e-5,
        "glucose_conc": GC,
        "thickness": 0.04,
    }

# -----------------------------
# MONTE CARLO SIMULATION PARAMETERS
# -----------------------------

def set_simulation_parameters():

    n_cores = 1
    n_photons = 10
    GC_a = [2]

    return n_cores, n_photons, GC_a