import matplotlib.pyplot as plt
import gc_ar.set_parameters as set_parameters
from gc_ar.plots import plot_variable_vs_angle
from gc_ar.multiple_photons import simulate_multiple_photon
import numpy as np
import gc_ar.results as results
from gc_ar.plots import plot_energy_distribution, plot_photon_paths
from gc_ar.monitors import absorption_monitor

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":

    #setting simulation parameters
    start= 0.5
    end= 3.5
    loop= 1
    name="GC"

    def init_process(x):
        import gc_ar.set_parameters as set_parameters
        set_parameters.set_material(name, x)


    for i in range (loop):
        #increase constant we are sweeping by one
        x=start+i*(end-start)/loop
        set_parameters.set_material(name, x)

        #call one multi-photon run and write results into the results file
        output = simulate_multiple_photon(set_parameters.get_material("n_photons"))
        results.conc_to_variable_vs_output(np.hstack(([x],output)))

    print(results.return_absorption_matrix()[:, 0].sum())
    print(np.sum(results.return_detected_energy()))
    print(np.sum(results.return_out_of_bound_energy()))

    #Plot everything
    #plot_variable_vs_angle(name)
    detector = {
        "cone_axis": np.array([0.0, 0, 1]),
        "cone_center": np.pi / 8,
        "r": 0.6
    }
    path=np.array(results.return_energy_matrix())
    plot_variable_vs_angle(name)
    plot_photon_paths([path],detector,sphere_radius=0.6)
    print(len(path))
