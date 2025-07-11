from single_photon import simulate_one_photon
from set_parameters import set_simulation_parameters
from computations import rotation_angle_calculation
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
import os,time
from itertools import repeat
#
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm


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

#    data = np.array(results)
#    sns.histplot(data, kde=True, stat="density", label="Data", bins=30)
#    xmin, xmax = plt.xlim()
#    x = np.linspace(xmin, xmax, 100)
#    p = norm.pdf(x, np.mean(data), np.std(data))
#    plt.plot(x, p, 'k', linewidth=2, label="Gaussian fit")
#    plt.legend()
#    plt.title("Histogram with Gaussian Fit")
#    plt.show()

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

    n_cores, n_photons, GC_a = set_simulation_parameters()

    angles = []

    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        angles = list(pool.map(simulate_one_gc, GC_a, repeat(n_photons)))

    plt.plot(GC_a, angles)
    plt.xlabel("Glucose concentration")
    plt.ylabel("Rotation Angle")
    plt.title("GC vs. Angle")
    plt.show()
#    plt.savefig(r"/home/ubuntu/results/gc_angles.png",dpi=300)