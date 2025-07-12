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
        "r": 1
    }

# -----------------------------
# MATERIAL PARAMETERS
# -----------------------------

parameters = {
        "mu_s": 1, #cm^-1
        "mu_a": .1, #cm^-1
        "g": 0.9, #unitless
        "n": 1.37, #unitless
        "alpha": 1e-5, #degree per decimeter per g/ml
        "GC": 3 # g/ml
        }

def set_material(key, value):
    if key not in parameters:
        raise KeyError(f"{key} is not a valid config key.")
    parameters[key] = value



def get_material(key):
    return parameters.get(key)

# -----------------------------
# MONTE CARLO SIMULATION PARAMETERS
# -----------------------------

def set_simulation_parameters():

    n_cores = 1
    n_photons = 1000
    GC_a = [2]

    return n_cores, n_photons, GC_a