import numpy as np

# -----------------------------
# MATERIAL PARAMETERS
# -----------------------------

parameters = {
        #medium and scattering particle parameters
        "mu_s": 1, #cm^-1
        "mu_a": .1, #cm^-1
        "g": 0.9, #unitless
        "n": 1.37, #unitless
        "alpha": 52.7, #degree per decimeter per g/ml
        "GC": 2, # g/ml

        #detector parameters
        # A detector is defined by an angle and the r

        "cone_axis": np.array([0.0, 0.0, 1]),
        "cone_center": np.pi / 8,
        "r": 10,

        #photon roulette parameter
        "chance": 0.1, #implementing p12. of energy conservation paper
        "energy_threshold": 1e-4 #threshold energy level for photon to leave
}

# -----------------------------
# INITIALIZE PHOTON
# -----------------------------
def initialize_photon():
    position = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    stokes = np.array([1.0, 1.0, 0.0, 0.0])  # linearly polarized along x
    energy = 1
    x,y,z= position
    energy_m = [energy,x,y,z]
    return position, direction, stokes, energy, energy_m

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
    n_photons = 10

    return n_cores, n_photons