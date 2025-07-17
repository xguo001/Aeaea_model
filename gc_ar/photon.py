import numpy as np
import gc_ar.set_parameters as set_parameters
from gc_ar.computations import rotation_matrix_glucose, change_direction, compute_phi, mie_scattering_matrix_rayleigh, rotation_matrix_phi, compute_ca1, compute_transmitted_direction, RFresnel,mu_a_circular_dichroism

class Photon:
    def __init__(self, position, direction, stokes,energy):
        self.position = position
        self.position_hit_detector = [0,0,0]
        self.position_hit_boundary = [0,0,0]
        self.direction = direction
        self.stokes = stokes
        self.energy = energy
        self.total_path_length = 0.0
        self.total_step_count = 0.0
        self.died_detected = False
        self.alive = True

    def update_position(self, step_size):
        self.position += step_size * self.direction

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
        self.position_hit_boundary = pos_start + t_value * (self.position - pos_start)

    def reflection_transmission(self):
        #Will reflect or transmit based on recorded boundary hit position -- will return error if boundary not hit
        if np.array_equal(self.position_hit_boundary, [0,0,0]):
            raise Exception("PROBLEM: REFLECTION FUNCTION CALLED WHEN BOUNDARY NOT HIT")

        #calculate ca1
        ca1, normal_vector = compute_ca1(self.position_hit_boundary, self.direction,set_parameters.get_material("r"))

        #calculate ca2
        r, ca2 = RFresnel(set_parameters.get_material("n"),set_parameters.get_material("n1"),ca1)

#        position= photon.position_hit_boundary
#        direction = photon.direction
#        radius = get_material('r')
#        ca1,normal=compute_ca1(position,direction,radius)
#        n=get_material('n')
#        n1=get_material('n1')
#        r,ca2=RFresnel(n,n1,ca1)

        #decide whether photon is reflected
        if r <= np.random.rand():

            self.direction = compute_transmitted_direction(self.direction, normal_vector, set_parameters.get_material("n"), set_parameters.get_material("n1"), ca2, ca1 )
            print("this photon is reflected")

            return True

#            direction_new= compute_transmitted_direction(direction,normal,n,n1,ca2,ca1)
#            photon.direction = direction_new
#            photon.died_detected=True

        else:

            print("this photon went out of bound")

            return False