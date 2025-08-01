import initialize.set_parameters as set_parameters
from photon_journey.photon_journey import simulate_multiple_photon
import numpy as np
import initialize.results as results
from detectors_and_plots.plots import plot_variable_vs_angle, plot_absorbed_energy_vs_time
from matplotlib.backends.backend_pdf import PdfPages
import os

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":

    #setting simulation parameters
    start= 0.5
    end= 3.5
    loop= 7
    name="GC"
    desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    pdf_filename = os.path.join(desktop_path, 'absorbed_energy.pdf')

    def init_process(x):
        import initialize.set_parameters as set_parameters
        set_parameters.set_material(name, x)


    with PdfPages(pdf_filename) as pdf:
        # Plot everything to save to pdf

        for i in range (loop):
            #increase constant we are sweeping by one
            x=start+i*(end-start)/loop
            set_parameters.set_material(name, x)

            #call one multi-photon run and write results into the results file
            output = simulate_multiple_photon(set_parameters.get_material("n_photons"))
            results.conc_to_variable_vs_output(np.hstack(([x], output)))
            plot_absorbed_energy_vs_time(pdf_pages=pdf)

        plot_variable_vs_angle(pdf_pages = pdf,name = name)

    #path=np.array(results.return_energy_matrix())
    #plot_photon_paths([path],detector,sphere_radius=0.6)
