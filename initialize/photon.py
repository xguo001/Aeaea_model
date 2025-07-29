import numpy as np
import initialize.set_parameters as set_parameters
import initialize.results as results
from photon_journey.computations import rotation_matrix_glucose, change_direction, compute_phi, mie_scattering_matrix_rayleigh, rotation_matrix_phi, compute_ca1, \
    RFresnel,mu_a_circular_dichroism, compute_reflected_direction

class Photon:
    def __init__(self, position, direction, stokes,energy):
        self.position = position
        self.position_hit_detector = [0,0,0]
        self.direction = direction
        self.stokes = stokes
        self.energy = energy
        self.total_path_length = 0.0
        self.total_step_count = 0.0
        self.died_detected = False
        self.alive = True
        self.time_alive = 0.0

    def update_time_alive(self, added_distance):
        #take distance traveled and calculate time spent traveling
        #also add to global time as each photon travels (sequentially)
        self.time_alive += added_distance / set_parameters.get_material("light_speed")
        results.update_global_time(added_distance / set_parameters.get_material("light_speed"))

    def update_position(self, step_size):
        #In the single photon code, step_size could be + or - when we cut path at boundary
        #We will update time every time we update position.
        self.position += step_size * self.direction
        self.update_time_alive(step_size)

    def update_position_without_step(self, new_position):
        self.position = new_position

    def update_stokes_through_glucose_medium(self, s):
        # Rotate Angle and Polarize due to Glucose
        theta_glucose = set_parameters.get_material("alpha") * set_parameters.get_material("GC") * (s / 10)
        D = rotation_matrix_glucose(theta_glucose)
        self.stokes = D @ self.stokes

    def update_mie_scattering_and_rotate(self):
        # Scatter and Polarize due to scattering particle
        # Change direction
        dir1= change_direction() #change direction uses random number and constant g that it's calling directly
        phi= compute_phi(self.direction, dir1)
        self.direction = dir1

        # Change Stoke
        M= mie_scattering_matrix_rayleigh(phi)
        D= rotation_matrix_phi(phi) @ M @ rotation_matrix_phi(phi)
        self.stokes = D @ self.stokes

    def update_direction(self, matrix):
        self.direction = matrix @ self.direction

    def update_energy(self,s):
        stokes = self.stokes
        mu_at=mu_a_circular_dichroism(stokes)
        self.energy = self.energy * np.exp(-mu_at * s)

    def update_energy_without_decay(self,energy):
        self.energy = energy

    def update_path_length(self, step_length):
        self.total_path_length += step_length

    def update_step_count(self):
        self.total_step_count += 1

    def update_died_detected(self):
        self.died_detected = True

    def update_position_hit_detector(self,pos_start, t_value):
        self.position_hit_detector = pos_start + t_value * (self.position - pos_start)

    def update_position_hit_boundary(self,pos_start, t_value):
        self.position = pos_start + t_value * (self.position - pos_start)

    def reflection_transmission(self):
        #Will reflect or transmit based on recorded boundary hit position -- return value is whether photon is reflected
        #Will check whether the new direction is facing outside of the sphere or inside of the sphere. If outside, terminate the photon
        #Why this works:
        # The position vector P points from the origin (center) to the point on the sphere, so it always points radially outward. The dot product measures how much your direction vector D aligns with this outward radial direction. A positive dot product means they point in similar directions (outward), while a negative dot product means they point in opposite directions (inward).

        #calculate ca1
        ca1, normal_vector = compute_ca1(self.position, self.direction, set_parameters.get_material("r"))

        #calculate ca2
        r, ca2 = RFresnel(set_parameters.get_material("n"), set_parameters.get_material("n1"), ca1)

        #decide whether photon is reflected
        if np.random.rand() <= r:

            self.direction = compute_reflected_direction(self.direction, normal_vector, set_parameters.get_material("n"), set_parameters.get_material("n1"))

            #print ("I got reflected, ", self.direction)

            # decide whether the new direction points inwards
            # If P · D < 0: The direction points inward (toward the center)
            if np.dot(self.direction, self.position) < 0:

                #print ("my reflection was inwards")

                return True

            # If P · D > 0: The direction points outward (away from the center)
            # If P · D = 0: The direction is tangential to the sphere
            else:

                #print ("my reflection was outwards")

                return False

        else:
            #in this branch photon is transmitted.

            #print ("I'm transmitted")

            return False

            # # Photon is transmitted (only if ca2 is not None)
            # if ca2 is not None:
            #     self.direction = compute_transmitted_direction(self.direction, normal_vector, set_parameters.get_material("n"), set_parameters.get_material("n1"), ca2, ca1)
            #     print ("I'm turning this way: ", self.direction)
            #     print ("The dot plot is giving this answer: ", np.dot(self.direction, self.position))
            #
            # else:
            #     # Total internal reflection case - should not happen since r=1.0
            #     from photon_journey.computations import compute_reflected_direction
            #     self.direction = compute_reflected_direction(self.direction, normal_vector, set_parameters.get_material("n"), set_parameters.get_material("n1"))
            #     return False