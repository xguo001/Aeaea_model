
from gc_ar.computations import energy_decay, compute_phi, rotation_matrix_phi, rotation_matrix_glucose, mie_scattering_matrix_rayleigh
import gc_ar.set_parameters as set_parameters
from gc_ar.detectors import detect_boundary, detect_photon_v2, photon_roulette
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
    pos_start = pos.copy()
    energy_m0= energy_m.copy()
    x0, y0, z0 = pos_start
    mu_t = set_parameters.get_material("mu_s") + set_parameters.get_material("mu_a")
    mu_a = set_parameters.get_material("mu_a")
    total_path_length = 0
    alive = True
    step_counter = 0
    gc = set_parameters.get_material("GC")
    abs_m=np.empty((0,4))


    while alive:
        pos_start = pos.copy()
        x0,y0,z0 = pos_start

        # Travel step
        s = -np.log(np.random.rand()) / mu_t
        total_path_length += s
        pos += dir * s
        step_counter += 1

        # Energy decay
        energy_0=energy
        energy = energy_decay(energy,mu_t, s)
        x,y,z = pos
        energy_m = [energy,x,y,z]
        abs_m = np.vstack([abs_m, [energy-energy_0,x,y,z]])
        #print ("abs"+str(abs_m))

        # Look to see if unalived
        if energy <= set_parameters.get_material("energy_threshold"):
            energy = photon_roulette(energy,set_parameters.get_material("chance"))
            if energy == 0:
                alive = False
                break

        # Look to see if detected
        detected_photon, t_value = detect_photon_v2(pos_start,pos,set_parameters.get_material("cone_axis"),set_parameters.get_material("cone_center"),set_parameters.get_material("r"))

        if detected_photon:
            #Recalculate path to only pre-detector path
            total_path_length -= (1-t_value)*s

        # Look to see if out of bound
        # !!!!!!! <- boundary currently set to have the same radius as where detector is located
        #boundary detection must come after photon detection because we have same radius for both.
        if not detected_photon and detect_boundary(pos,set_parameters.get_material("r")):
            alive = False
            break

        # Rotate Angle and Polarize due to Glucose
        theta_glucose = set_parameters.get_material("alpha") * gc * (s/10)
        D = rotation_matrix_glucose(theta_glucose)
        stokes = D @ stokes

        #If detected, we record. If not detected, we scatter and travel again
        if detected_photon:

            return alive, step_counter, total_path_length, stokes, abs_m

        # Scatter and Polarize due to scattering particle
        dir1= change_direction(g=0.9)
        phi= compute_phi(dir, dir1)
        dir=dir1
        M= mie_scattering_matrix_rayleigh(phi)
        D= rotation_matrix_phi(phi) @ M @ rotation_matrix_phi(phi)
        stokes = D @ stokes

    #If we reach here, alive = false, so the rest of return doesn't matter
    return alive, step_counter, total_path_length, stokes, abs_m

