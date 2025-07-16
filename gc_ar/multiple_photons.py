import gc_ar.set_parameters as set_parameters
import gc_ar.photon as photon
from gc_ar.single_photon import simulate_one_photon
from gc_ar.computations import rotation_angle_calculation
import numpy as np
import os,time


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

    return np.mean(results), np.mean(step_counters), np.mean(path_lengths_collector)

# -----------------------------
# WRAPPER AROUND SINGLE GLUCOSE LEVEL SIMUATION WITH MULTITHREADING AND ANGLE ROTATION PRINTOUT
# -----------------------------
def simulate_one_gc(n_photons):
    gc = set_parameters.get_material("GC")
    start = time.perf_counter()
    pid = os.getpid()
    print(f"[PID {pid}] processing GC = {gc} begins", flush=True)  # live announcement

    angle, steps, pathlength = simulate_multiple_photon(n_photons)

    elapsed = (time.perf_counter() - start)/60
    print(f"[PID {pid}] processing GC = {gc} ended in {elapsed:.2f} minutes", flush=True)  # live announcement
    print(f"GC = {gc}, Angle = {angle}, Average steps = {steps}")

    return angle, steps, pathlength
