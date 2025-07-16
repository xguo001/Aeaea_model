import matplotlib.pyplot as plt
import gc_ar.set_parameters as set_parameters
from gc_ar.multiple_photons import simulate_one_gc
import numpy as np
import gc_ar.results as results
from gc_ar.plots import plot_energy_distribution
from gc_ar.monitors import absorption_monitor

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":

    print( set_parameters.get_material("GC") )
    start= 0.5
    end= 3.5
    loop= 1
    name="GC"
    cons=[]
    output_angles=[]
    output_steps=[]
    output_length=[]
    output_paths = []



    def init_process(x):
        import gc_ar.set_parameters as set_parameters
        set_parameters.set_material(name, x)


    for i in range (loop):
        x=start+i*(end-start)/loop
        cons.append(x)
        print("set x"+str(x))
        set_parameters.set_material(name,x)
        print("get x"+str(set_parameters.get_material(name)))
        n_cores, n_photons = set_parameters.set_simulation_parameters()
        output=simulate_one_gc(n_photons)
        output_angles.append(output[0])
        output_steps.append(output[1])
        output_paths.append(output[2])

    print(results.return_absorption_matrix()[:, 0].sum())
    print(np.sum(results.return_detected_energy()))
    print(np.sum(results.return_out_of_bound_energy()))

    #print(results_angles)
    #plt.plot(cons,results_steps,"-o")
    #plt.ylabel("Steps")
    #plt.xlabel(name+" vs. Steps")
    #plt.savefig("/Users/xwguo/Results/steps_vs"+str(name)+'.png',dpi=300)
    #plt.close()

    plt.plot(cons, output_angles, "-o")
    plt.ylabel("angles")
    plt.xlabel(name + " vs. angles")
    plt.show()
    #plt.savefig("/Users/xwguo/Results/angles_vs" + str(name) + '.png', dpi=300)
    #plt.close()
    #plot_energy_distribution(results.return_absorption_matrix())

    #plt.plot(cons, results_paths, "-o")
    #plt.ylabel("pathlength")
    #plt.xlabel(name + " vs. pathlength")
    #plt.savefig("/Users/xwguo/Results/pathlength_vs" + str(name) + '.png', dpi=300)
    #plt.close()