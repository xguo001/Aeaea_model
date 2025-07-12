import matplotlib.pyplot as plt

import set_parameters
from concurrent.futures import ProcessPoolExecutor
from multiple_photons import simulate_one_gc
from itertools import repeat
from plot_them import plot_GC_vs_angles_plot

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":

    print( set_parameters.get_material("alpha") )
    start=1
    end=60
    np=1
    name="alpha"
    cons=[]
    results_angles=[]
    results_steps=[]
    results_length=[]
    results_paths = []


    def init_process(x):
        import set_parameters
        set_parameters.set_material(name, x)


    for i in range (np):
        x=start+i*(end-start)/np
        cons.append(x)
        print("set x"+str(x))
        set_parameters.set_material(name,x)
        print("get x"+str(set_parameters.get_material(name)))
        n_cores, n_photons, GC_a = set_parameters.set_simulation_parameters()
        results=simulate_one_gc(GC_a[0],n_photons)
        results_angles.append(results[0])
        results_steps.append(results[1])
        results_paths.append(results[2])
    print(results_angles)
    plt.plot(cons,results_steps,"-o")
    plt.ylabel("Steps")
    plt.xlabel(name+" vs. Steps")
    plt.savefig("/Users/xwguo/Results/steps_vs"+str(name)+'.png',dpi=300)
    plt.close()

    plt.plot(cons, results_angles, "-o")
    plt.ylabel("angles")
    plt.xlabel(name + " vs. angles")
    plt.savefig("/Users/xwguo/Results/angles_vs" + str(name) + '.png', dpi=300)
    plt.close()

    plt.plot(cons, results_paths, "-o")
    plt.ylabel("pathlength")
    plt.xlabel(name + " vs. pathlength")
    plt.savefig("/Users/xwguo/Results/pathlength_vs" + str(name) + '.png', dpi=300)
    plt.close()