import gc_ar.results as results
from gc_ar.computations import energy_decay, compute_phi, rotation_matrix_phi, rotation_matrix_glucose, mie_scattering_matrix_rayleigh, cut_path_at_boundary
import gc_ar.set_parameters as set_parameters
from gc_ar.detectors import detect_boundary, detect_photon_v2, photon_roulette
from gc_ar.computations import mid_point
import numpy as np
    
def change_direction(g):
    """ Sample scattering direction using Henyey-Greenstein phase function """
    cos_theta = (1 / (2 * g)) * (1 + g ** 2 - ((1 - g ** 2) / (1 - g + 2 * g * np.random.rand())) ** 2)
    sin_theta = np.sqrt(1 - cos_theta ** 2)
    phi = 2 * np.pi * np.random.rand()
    dx = sin_theta * np.cos(phi)
    dy = sin_theta * np.sin(phi)
    dz = cos_theta
    return np.array([dx, dy, dz])

# -----------------------------
# PHOTON PROPAGATION AND MATRIX OPERATIONS (Eq. 5)
# -----------------------------
def simulate_one_photon():

    # Initialize variables

    pos, dir, stokes, energy, energy_m = set_parameters.initialize_photon()
    mu_t = set_parameters.get_material("mu_s") + set_parameters.get_material("mu_a")
    mu_a = set_parameters.get_material("mu_a")
    total_path_length = 0
    alive = True
    step_counter = 0
    gc = set_parameters.get_material("GC")

    while alive:
        pos_start = pos.copy()

        # Travel step
        s = -np.log(np.random.rand()) / mu_t
        total_path_length += s
        pos += dir * s
        step_counter += 1

        # Energy decay
        print ("I started with: ",energy)
        energy_0 = energy
        energy = energy_decay(energy,mu_t, s)

        #Branch #1: Look to see if unalived
        if energy <= set_parameters.get_material("energy_threshold"):
            print ("I'm in branch #1")
            energy = photon_roulette(energy,set_parameters.get_material("chance"))
            if energy == 0:
                alive = False
                #concatenante to results the energy level and midpoint of the travel
                results.conc_to_absorption_matrix(np.hstack(([energy_0 - energy],mid_point(pos_start,pos))))
                break

        #Branch #2: Look to see if detected
        detected_photon, t_value = detect_photon_v2(pos_start,pos,set_parameters.get_material("cone_axis"),set_parameters.get_material("cone_center"),set_parameters.get_material("r"))

        if detected_photon:
            print ("I'm in branch #2")
            #Recalculate path to only pre-detector path
            total_path_length -= (1-t_value)*s

            #scale energy and deposit energy at path midpoint
            results.conc_to_absorption_matrix(np.hstack(([(energy_0 - energy)*t_value],mid_point(pos_start,pos_start+t_value*(pos-pos_start)))))

            #deposit remaining energy to the detector
            results.conc_to_detected_energy([energy_0 - (energy_0 - energy)*t_value])

            # Rotate Angle and Polarize due to Glucose
            theta_glucose = set_parameters.get_material("alpha") * gc * (s/10)
            D = rotation_matrix_glucose(theta_glucose)
            stokes = D @ stokes

            return alive, step_counter, total_path_length, stokes

        #Branch #3: Look to see if out of bound
        # !!!!!!! <- boundary currently set to have the same radius as where detector is located
        #boundary detection must come after photon detection because we have same radius for both.
        if not detected_photon and detect_boundary(pos,set_parameters.get_material("r")):
            alive = False
            print ("I'm in branch #3")
            #Need to cut path off at the boundary to calculate energy level
            #t_value = length inside / total length
            t_value = cut_path_at_boundary(pos_start, pos, set_parameters.get_material("r"))

            #Scale energy and deposit at midpoint inside
            results.conc_to_absorption_matrix(np.hstack(([(energy_0 - energy)*t_value],mid_point(pos_start,pos_start+t_value*(pos-pos_start)))))

            #Deposit remainder of energy onto the boundary
            results.conc_to_out_of_bound_energy(energy_0 - (energy_0 - energy)*t_value)

            break

        #Branch 4: keep traveling
        print ("I'm in branch #4")

        # Deposit energy
        results.conc_to_absorption_matrix(np.hstack(([energy_0 - energy],mid_point(pos_start,pos))))

        # Rotate Angle and Polarize due to Glucose
        theta_glucose = set_parameters.get_material("alpha") * gc * (s/10)
        D = rotation_matrix_glucose(theta_glucose)
        stokes = D @ stokes

        # Scatter and Polarize due to scattering particle
        dir1= change_direction(g=0.9)
        phi= compute_phi(dir, dir1)
        dir=dir1
        M= mie_scattering_matrix_rayleigh(phi)
        D= rotation_matrix_phi(phi) @ M @ rotation_matrix_phi(phi)
        stokes = D @ stokes

    #If we reach here, alive = false, so the rest of return doesn't matter
    return alive, step_counter, total_path_length, stokes

