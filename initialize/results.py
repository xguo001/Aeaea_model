import numpy as np

# -----------------------------
# STORING MODEL OUTPUT
# -----------------------------

results = {

    #Absorption matrices
    "absorption_matrix": np.empty((0,5)),
    "detected_energy": np.empty((0,2)),
    "out_of_bound_energy": np.empty((0, 2)),

    #Energy level matrix
    #The energy matrix is not currently fully implemented. It was used for debugging.
    'energy_matrix': np.empty((0, 4)),

    #Sweeped variable vs. output matrix
    #variable, np.mean(results), np.mean(step_counters), np.mean(path_lengths_collector), death_counters
    "variable_vs_output": np.empty((0, 5)),

    #Global time
    "global_time": 0,

}

def update_global_time(added_time):
    results['global_time'] += added_time

def conc_to_energy_matrix(energy_m):
    #takes energy level, x, y, z and concatenante into the results
    results['energy_matrix'] = np.vstack([results['energy_matrix'], energy_m])

def conc_to_absorption_matrix(energy_position):
    #takes energy level, x, y, z and concatenante into the results
    new_row = np.array([energy_position[0], energy_position[1], energy_position[2], energy_position[3],results['global_time']])
    results["absorption_matrix"] = np.vstack([results["absorption_matrix"], new_row])

def conc_to_detected_energy(energy):
    new_row = np.array([energy, results['global_time']])
    results["detected_energy"] = np.vstack([results["detected_energy"], new_row])

def conc_to_out_of_bound_energy(energy):
    new_row = np.array([energy, results['global_time']])
    results["out_of_bound_energy"] = np.vstack([results["out_of_bound_energy"], new_row])

def conc_to_variable_vs_output(one_row):
    #takes [0] variable being sweeped, [1] average rotation angles, [2] average path length [3] total number of death

    results["variable_vs_output"] = np.vstack([results["variable_vs_output"], one_row])

def return_global_time():
    return results['global_time']

def return_absorption_matrix():
    return results["absorption_matrix"]

def return_energy_matrix():
    return results["energy_matrix"]

def return_detected_energy():
    return results["detected_energy"]

def return_out_of_bound_energy():
    return results["out_of_bound_energy"]

def return_variable_vs_output():
    return results["variable_vs_output"]