import numpy as np

# -----------------------------
# STORING MODEL OUTPUT
# -----------------------------

results = {
    "absorption_matrix": np.empty((0,4)),
    "detected_energy": np.empty((0,1)),
    "out_of_bound_energy": np.empty((0, 1)),
}

def conc_to_absorption_matrix(energy_position):
    #takes energy level, x, y, z and concatenante into the results

    results["absorption_matrix"] = np.vstack([results["absorption_matrix"], energy_position])
    print("adding to absorption matrix", energy_position[0])

def conc_to_detected_energy(energy):
    #takes a scalar and concatenante

    results["detected_energy"] = np.vstack([results["detected_energy"], energy])
    print("adding to energy",energy)

def conc_to_out_of_bound_energy(energy):
    #takes a scalar and concatenante

    results["out_of_bound_energy"] = np.vstack([results["out_of_bound_energy"], energy])
    print ("adding to out-of-bound energy",energy)

def return_absorption_matrix():
    return results["absorption_matrix"]

def return_detected_energy():
    return results["detected_energy"]

def return_out_of_bound_energy():
    return results["out_of_bound_energy"]