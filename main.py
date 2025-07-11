from set_parameters import set_simulation_parameters
from concurrent.futures import ProcessPoolExecutor
from multiple_photons import simulate_one_gc
from itertools import repeat
from plot_them import plot_GC_vs_angles_plot

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":

    n_cores, n_photons, GC_a = set_simulation_parameters()

    angles = []

    with ProcessPoolExecutor(max_workers=n_cores) as pool:
        angles = list(pool.map(simulate_one_gc, GC_a, repeat(n_photons)))

    plot_GC_vs_angles_plot(GC_a, angles)