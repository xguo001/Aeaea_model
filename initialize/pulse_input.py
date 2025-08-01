import numpy as np
import initialize.set_parameters as set_parameters

def sample_launch_time():
    ti=np.random.normal(loc=set_parameters.get_material("pulse_peak_time"), scale=set_parameters.get_material("pulse_width") / 2.355)
    print("-------------------------Sample launch time-"+str(ti))
    return ti


