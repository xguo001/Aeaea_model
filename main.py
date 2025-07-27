import initialize.set_parameters as set_parameters
from detectors_and_plots.plots import plot_variable_vs_angle
from photon_journey.multiple_photons import simulate_multiple_photon
import numpy as np
import initialize.results as results
from detectors_and_plots.plots import plot_photon_paths

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
        import initialize.set_parameters as set_parameters
        set_parameters.set_material(name, x)


    for i in range (loop):
        #increase constant we are sweeping by one
        x=start+i*(end-start)/loop
        set_parameters.set_material(name, x)

        #call one multi-photon run and write results into the results file
        output = simulate_multiple_photon(set_parameters.get_material("n_photons"))
        results.conc_to_variable_vs_output(np.hstack(([x], output)))


    #Plot everything
    #plot_variable_vs_angle(name)
    detector = {
        "cone_axis": np.array([0.0, 0, 1]),
        "cone_center": np.pi / 8,
        "r": 0.6
    }
    path=np.array(results.return_energy_matrix())
    plot_variable_vs_angle(name)
#    plot_photon_paths([path],detector,sphere_radius=0.6)
