import numpy as np

# -----------------------------
# STORING MODEL OUTPUT
# -----------------------------

results = {

    #Energy matrices
    "absorption_matrix": np.empty((0,4)),
    "detected_energy": np.empty((0,1)),
    "out_of_bound_energy": np.empty((0, 1)),

    #Sweeped variable vs. output matrix
    #variable, np.mean(results), np.mean(step_counters), np.mean(path_lengths_collector), death_counters
    "variable_vs_output": np.empty((0, 5)),
}

def conc_to_absorption_matrix(energy_position):
    #takes energy level, x, y, z and concatenante into the results

    results["absorption_matrix"] = np.vstack([results["absorption_matrix"], energy_position])

def conc_to_detected_energy(energy):
    #takes a scalar and concatenante

    results["detected_energy"] = np.vstack([results["detected_energy"], energy])

def conc_to_out_of_bound_energy(energy):
    #takes a scalar and concatenante

    results["out_of_bound_energy"] = np.vstack([results["out_of_bound_energy"], energy])

def conc_to_variable_vs_output(one_row):
    #takes [0] variable being sweeped, [1] average rotation angles, [2] average path length [3] total number of death

    results["variable_vs_output"] = np.vstack([results["variable_vs_output"], one_row])

def return_absorption_matrix():
    return results["absorption_matrix"]

def return_detected_energy():
    return results["detected_energy"]

def return_out_of_bound_energy():
    return results["out_of_bound_energy"]

def return_variable_vs_output():
    return results["variable_vs_output"]