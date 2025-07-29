import initialize.results as results
import numpy as np


def absorption_monitor(position,width,length,height,absorbed_energy_matrix):
    #take the matrix of absorbed energy over the space (in an array of [e,x,y,z])
    #and sum total energy over the monitored space (given by cube center, w,l,h)

    total_energy = 0
    x,y,z=position

    x_bl=x-width/2
    y_bl=y-length/2
    z_bl=z-height/2

    x_br = x + width / 2
    y_br = y + length / 2
    z_br = z + height / 2

    for absorption_event in absorbed_energy_matrix:
        if x_bl <= absorption_event[1] <= x_br:
            if y_bl <= absorption_event[2] <= y_br:
                if z_bl <= absorption_event[3] <= z_br:
                    total_energy += absorption_event[0]

    return total_energy

def return_total_absoprtions():
    # Calculate energy values
    total_absorbed = results.return_absorption_matrix()[:, 0].sum()
    total_detected = results.return_absorption_matrix()[:, 0].sum()
    total_out_of_bound = results.return_out_of_bound_energy()[:, 0].sum()

    return (np.array([total_absorbed,total_detected,total_out_of_bound]))

