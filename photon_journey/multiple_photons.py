import initialize.set_parameters as set_parameters
from photon_journey.single_photon import simulate_one_photon
from photon_journey.computations import rotation_angle_calculation
import numpy as np


# -----------------------------
# TAKE A GLUCOSE LEVEL AND SIMULATE N NUMBER OF PHOTONS AND RETURN THE MEAN
# -----------------------------
def simulate_multiple_photon(n_photons):
    results = []
    step_counters = []
    path_lengths_collector =[]
    death_counters = 0
    gc = set_parameters.get_material("GC")

    for _ in range(n_photons):
        this_photon = simulate_one_photon()

        if this_photon.died_detected:
            results.append(rotation_angle_calculation(gc,this_photon.total_path_length))
            step_counters.append(this_photon.total_path_length)
            path_lengths_collector.append(this_photon.total_path_length)
        else:
            death_counters += 1

    print ("This many died", death_counters)
    print ("Average path length", np.mean(path_lengths_collector))

    return np.mean(results), np.mean(step_counters), np.mean(path_lengths_collector), death_counters
