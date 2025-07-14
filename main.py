import matplotlib.pyplot as plt

import gc_ar.set_parameters as set_parameters
from gc_ar.multiple_photons import simulate_one_gc

# -----------------------------
# RUN SIMULATION BASED ON NUMBER OF PHOTONS AND GLUCOSE LEVEL INPUTS, WRITE OUT TO A GRAPH
# -----------------------------
if __name__ == "__main__":

    print( set_parameters.get_material("GC") )
    start= 0.5
    end= 3.5
    loop= 7
    name="GC"
    cons=[]
    results_angles=[]
    results_steps=[]
    results_length=[]
    results_paths = []


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
        results=simulate_one_gc(n_photons)
        results_angles.append(results[0])
        results_steps.append(results[1])
        results_paths.append(results[2])

    #print(results_angles)
    #plt.plot(cons,results_steps,"-o")
    #plt.ylabel("Steps")
    #plt.xlabel(name+" vs. Steps")
    #plt.savefig("/Users/xwguo/Results/steps_vs"+str(name)+'.png',dpi=300)
    #plt.close()

    plt.plot(cons, results_angles, "-o")
    plt.ylabel("angles")
    plt.xlabel(name + " vs. angles")
    plt.show()
    #plt.savefig("/Users/xwguo/Results/angles_vs" + str(name) + '.png', dpi=300)
    plt.close()

    #plt.plot(cons, results_paths, "-o")
    #plt.ylabel("pathlength")
    #plt.xlabel(name + " vs. pathlength")
    #plt.savefig("/Users/xwguo/Results/pathlength_vs" + str(name) + '.png', dpi=300)
    #plt.close()