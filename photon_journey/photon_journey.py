import initialize.results as results
from photon_journey.computations import cut_path_at_boundary
import initialize.set_parameters as set_parameters
from detectors_and_plots.detectors import detect_boundary, detect_photon_v2, photon_roulette
from photon_journey.computations import mid_point, rotation_angle_calculation
import numpy as np
from photon_journey.launch_photons import launch_a_photon



# -----------------------------
# PHOTON PROPAGATION AND MATRIX OPERATIONS (Eq. 5)
# -----------------------------
def simulate_one_photon():

    # Initialize variables
    #launch_time = np.random.normal(loc=5.5e-5, scale=1.5e-5)
    #launch_time = np.clip(launch_time, 1e-5, 1e-4)
    launch_time=0
    this_photon= launch_a_photon(beam_radius=0.2,angle=0,z=0,launch_time=launch_time,pulse=True)
    mu_t = set_parameters.get_material("mu_s") + set_parameters.get_material("mu_a")

    while this_photon.alive:
        pos_start = this_photon.position.copy()

        # Travel step
        s = -np.log(np.random.rand()) / mu_t
        this_photon.update_path_length(s)
        this_photon.update_position(s)
        this_photon.update_step_count()



        # Energy decay
        energy_0 = this_photon.energy
        this_photon.update_energy(s)

        #save results to energy matrix
        #results.conc_to_energy_matrix(np.hstack(([this_photon.energy], this_photon.position)))

        #Branch #1: Look to see if unalived

        if this_photon.energy <= set_parameters.get_material("energy_threshold"):

            #this is the photon roulette
            this_photon.update_energy_without_decay(photon_roulette(this_photon.energy,
                                                                    set_parameters.get_material("chance")))

            if this_photon.energy == 0:
                this_photon.alive = False
                #concatenante to results the energy level and midpoint of the travel
                results.conc_to_absorption_matrix(np.hstack(([energy_0 - this_photon.energy], mid_point(pos_start, this_photon.position), this_photon.time_alive)))
                #We have to break the while loop, so we don't have to check energy == 0 in later branches
                return this_photon

            #If photon survives photon roulette, then it keeps traveling like normal

        #Branch #2: Look to see if detected
        #t_value tells you what portion of the trajectory is inside vs. beyond the detector
        detected_photon, t_value = detect_photon_v2(pos_start, this_photon.position,
                                                    set_parameters.get_material("cone_axis"),
                                                    set_parameters.get_material("cone_center"),
                                                    set_parameters.get_material("r"))

        if detected_photon:

            this_photon.alive = False
            #Recalculate path to only pre-detector path
            this_photon.update_path_length(-(1-t_value)*s)

            #update position where detector is hit
            this_photon.update_position_hit_detector(pos_start,t_value)

            #scale energy and deposit energy at path midpoint
            results.conc_to_absorption_matrix(np.hstack(([(energy_0 - this_photon.energy) * t_value], mid_point(pos_start, this_photon.position_hit_detector), this_photon.time_alive)))

            #deposit remaining energy to the detector
            results.conc_to_detected_energy(np.hstack((energy_0 - (energy_0 - this_photon.energy) * t_value, this_photon.time_alive)))

            # Rotate Angle and Polarize due to Glucose
            this_photon.update_stokes_through_glucose_medium(t_value * s)

            #update this photon to detected
            this_photon.update_died_detected()

        #Branch #3: Look to see if out of bound
        # !!!!!!! <- boundary currently set to have the same radius as where detector is located
        #boundary detection must come after photon detection because we have same radius for both.

        hit_boundary = detect_boundary(this_photon.position, set_parameters.get_material("r"))

        if not detected_photon and hit_boundary:
            #Need to cut path off at the boundary to calculate energy level
            #t_value = length inside / total length
            t_value = cut_path_at_boundary(pos_start, this_photon.position, set_parameters.get_material("r"))


            #update position to where boundary is hit
            this_photon.update_position_hit_boundary(pos_start,t_value)

            #Scale energy and deposit at midpoint between start point and where boundary is
            results.conc_to_absorption_matrix(np.hstack(([(energy_0 - this_photon.energy) * t_value], mid_point(pos_start, this_photon.position),this_photon.time_alive)))
            #results.conc_to_energy_matrix(np.hstack(([this_photon.energy], this_photon.position)))
            reflected = this_photon.reflection_transmission()

            if reflected:
                #The direction of this photon has already been changed
                #We need to update this photon's position, energy and path length to the boundary point and let it go from there
                this_photon.update_energy_without_decay(energy_0 - (energy_0 - this_photon.energy)*t_value) #for energy formula, recall we have decay along the route and what would have been deposited to the boundary
                this_photon.update_path_length(-(1 - t_value) * s)

                #the energy matrix is not currently fully implemented. It was used for debugging.
                results.conc_to_energy_matrix(np.hstack(([this_photon.energy], this_photon.position)))

            else:
                #if not reflected, then deposit remaining energy and this photon's journey ends
                #Deposit remainder of energy onto the boundary (the energy for the journey to the boundary has already been deposited)
                this_photon.alive = False
                results.conc_to_out_of_bound_energy(np.hstack([energy_0 - (energy_0 - this_photon.energy) * t_value, this_photon.time_alive]))

        #Branch 4: keep traveling
        if not detected_photon and not hit_boundary:

            # Deposit energy
            results.conc_to_absorption_matrix(np.hstack(([energy_0 - this_photon.energy], mid_point(pos_start, this_photon.position), this_photon.time_alive)))

            # Rotate Angle and Polarize due to Glucose
            this_photon.update_stokes_through_glucose_medium(s)

            # Scatter and Polarize due to scattering particle
            this_photon.update_mie_scattering_and_rotate()

    #If we reach here, alive = false, and died detected = false, nothing to return since it's all in the photon class

    return  this_photon

# -----------------------------
# TAKE A GLUCOSE LEVEL AND SIMULATE N NUMBER OF PHOTONS AND RETURN THE MEAN
# -----------------------------
def simulate_multiple_photon(n_photons):
    results = []
    step_counters = []
    path_lengths_collector =[]
    death_counters = 0
    gc = set_parameters.get_material("GC")

    for _ in range(n_photons):
        this_photon = simulate_one_photon()

        if this_photon.died_detected:
            results.append(rotation_angle_calculation(gc,this_photon.total_path_length))
            step_counters.append(this_photon.total_path_length)
            path_lengths_collector.append(this_photon.total_path_length)
        else:
            death_counters += 1

    print ("This many died", death_counters)
    print ("Average path length", np.mean(path_lengths_collector))

    return np.mean(results), np.mean(step_counters), np.mean(path_lengths_collector), death_counters
