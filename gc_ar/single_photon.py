import gc_ar.results as results
from gc_ar.computations import energy_decay, compute_phi, rotation_matrix_phi, rotation_matrix_glucose, mie_scattering_matrix_rayleigh, cut_path_at_boundary, change_direction
import gc_ar.set_parameters as set_parameters
from gc_ar.detectors import detect_boundary, detect_photon_v2, photon_roulette
from gc_ar.computations import mid_point
import numpy as np
from gc_ar.photon import Photon
from gc_ar.launch_photons import launch_a_photon


# -----------------------------
# PHOTON PROPAGATION AND MATRIX OPERATIONS (Eq. 5)
# -----------------------------
def simulate_one_photon():

    # Initialize variables
    this_photon= launch_a_photon(beam_radius=0.2,angle=0,z=0)
    mu_t = set_parameters.get_material("mu_s") + set_parameters.get_material("mu_a")
    alive = True

    while alive:
        pos_start = this_photon.position.copy()

        # Travel step
        s = -np.log(np.random.rand()) / mu_t
        this_photon.update_path_length(s)
        this_photon.update_position(s)
        this_photon.update_step_count()

        # Energy decay
        energy_0 = this_photon.energy
        this_photon.update_energy(s)

        #Branch #1: Look to see if unalived
        if this_photon.energy <= set_parameters.get_material("energy_threshold"):

            #this is the photon roulette
            this_photon.update_energy_without_decay(photon_roulette(this_photon.energy,set_parameters.get_material("chance")))

            if this_photon.energy == 0:
                alive = False
                #concatenante to results the energy level and midpoint of the travel
                results.conc_to_absorption_matrix(np.hstack(([energy_0 - this_photon.energy],mid_point(pos_start,this_photon.position))))

                return this_photon

        #Branch #2: Look to see if detected
        #t_value tells you what portion of the trajectory is inside vs. beyond the detector
        detected_photon, t_value = detect_photon_v2(pos_start,this_photon.position,set_parameters.get_material("cone_axis"),set_parameters.get_material("cone_center"),set_parameters.get_material("r"))

        if detected_photon:
            alive = False
            #Recalculate path to only pre-detector path
            this_photon.update_path_length(-(1-t_value)*s)

            #scale energy and deposit energy at path midpoint
            results.conc_to_absorption_matrix(np.hstack(([(energy_0 - this_photon.energy)*t_value],mid_point(pos_start,pos_start+t_value*(this_photon.position - pos_start)))))

            #deposit remaining energy to the detector
            results.conc_to_detected_energy([energy_0 - (energy_0 - this_photon.energy)*t_value])

            # Rotate Angle and Polarize due to Glucose
            this_photon.update_stokes_through_glucose_medium(s)

            this_photon.update_died_detected()

            return this_photon

        #Branch #3: Look to see if out of bound
        # !!!!!!! <- boundary currently set to have the same radius as where detector is located
        #boundary detection must come after photon detection because we have same radius for both.
        if not detected_photon and detect_boundary(this_photon.position,set_parameters.get_material("r")):
            alive = False
            #Need to cut path off at the boundary to calculate energy level
            #t_value = length inside / total length
            t_value = cut_path_at_boundary(pos_start, this_photon.position, set_parameters.get_material("r"))

            #Scale energy and deposit at midpoint inside
            results.conc_to_absorption_matrix(np.hstack(([(energy_0 - this_photon.energy)*t_value],mid_point(pos_start,pos_start+t_value*(this_photon.position-pos_start)))))

            #Deposit remainder of energy onto the boundary
            results.conc_to_out_of_bound_energy(energy_0 - (energy_0 - this_photon.energy)*t_value)

            return this_photon

        #Branch 4: keep traveling

        # Deposit energy
        results.conc_to_absorption_matrix(np.hstack(([energy_0 - this_photon.energy],mid_point(pos_start,this_photon.position))))

        # Rotate Angle and Polarize due to Glucose
        this_photon.update_stokes_through_glucose_medium(s)

        # Scatter and Polarize due to scattering particle
        this_photon.update_mie_scattering_and_rotate()

    #If we reach here, alive = false, and died detected = false, nothing to return since it's all in the photon class

    print ("THERE IS A PROBLEM")

    return  this_photon