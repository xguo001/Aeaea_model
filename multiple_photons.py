from single_photon import simulate_one_photon
from computations import rotation_angle_calculation
import numpy as np
import matplotlib.pyplot as plt
import os,time
#



# -----------------------------
# TAKE A GLUCOSE LEVEL AND SIMULATE N NUMBER OF PHOTONS AND RETURN THE MEAN
# -----------------------------
def simulate_multiple_photon(GC, n_photons):
    results = []
    step_counters = []
    path_lengths_collector =[]
    death_counters = 0

    for _ in range(n_photons):
        alive, step_counter, total_path_length, stokes = simulate_one_photon(GC)

        if alive:
            results.append(rotation_angle_calculation(GC,total_path_length))
            step_counters.append(step_counter)
            path_lengths_collector.append(total_path_length)
        else:
            death_counters += 1

    print ("This many died", death_counters)
    print ("Average path length", np.mean(path_lengths_collector))

    return np.mean(results), np.mean(step_counters), np.mean(path_lengths_collector)

# -----------------------------
# WRAPPER AROUND SINGLE GLUCOSE LEVEL SIMUATION WITH MULTITHREADING AND ANGLE ROTATION PRINTOUT
# -----------------------------
def simulate_one_gc(GC, n_photons):
    start = time.perf_counter()
    pid = os.getpid()
    print(f"[PID {pid}] processing GC = {GC} begins", flush=True)  # live announcement

    angle, steps, pathlength = simulate_multiple_photon(GC, n_photons)

    elapsed = (time.perf_counter() - start)/60
    print(f"[PID {pid}] processing GC = {GC} ended in {elapsed:.2f} minutes", flush=True)  # live announcement
    print(f"GC = {GC}, Angle = {angle}, Average steps = {steps}")

    return angle, steps, pathlength
