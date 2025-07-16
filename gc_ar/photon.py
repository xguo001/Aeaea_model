import numpy as np
import gc_ar.set_parameters as set_parameters
from gc_ar.computations import rotation_matrix_glucose, change_direction, compute_phi, mie_scattering_matrix_rayleigh, rotation_matrix_phi

class Photon:
    def __init__(self, position, direction, stokes,energy):
        self.position = position
        self.direction = direction
        self.stokes = stokes
        self.energy = energy
        self.total_path_length = 0.0
        self.total_step_count = 0.0
        self.died_detected = False

    def update_position(self, step_size):
        self.position += step_size * self.direction

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
        self.energy = self.energy * np.exp(-(set_parameters.get_material("mu_a")+set_parameters.get_material("mu_s")) * s)

    def update_energy_without_decay(self,energy):
        self.energy = energy

    def update_path_length(self, step_length):
        self.total_path_length += step_length

    def update_step_count(self):
        self.total_step_count += 1

    def update_died_detected(self):
        self.died_detected = True
