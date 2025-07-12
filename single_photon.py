from computations import energy_decay, compute_phi, rotation_matrix_phi, rotation_matrix_glucose, mie_scattering_matrix_rayleigh
from set_parameters import initialize_photon, setup_detector, get_material
from detectors import detect_boundary, detect_photon_v2
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

    detector = setup_detector()
    pos, dir, stokes, energy = initialize_photon()
    mu_t = get_material("mu_s") + get_material("mu_a")
    total_path_length = 0
    alive = True
    step_counter = 0
    gc = get_material("GC")

    while alive:
        pos_start = pos.copy()

        # Travel step
        s = -np.log(np.random.rand()) / mu_t
        total_path_length += s
        pos += dir * s
        step_counter += 1

        # Energy decay
        energy = energy_decay(energy,mu_t, s)

        # Look to see if unalived
        if energy <= 1e-4:
            alive = False
            break

        # Look to see if detected
        detected_photon, t_value = detect_photon_v2(pos_start,pos,detector["cone_axis"],detector["alpha"],detector["r"])

        if detected_photon:
            #Recalculate path to only pre-detector path
            total_path_length -= (1-t_value)*s

        # Look to see if out of bound
        # !!!!!!! <- boundary currently set to have the same radius as where detector is located
        #boundary detection must come after photon detection because we have same radius for both.
        if not detected_photon and detect_boundary(pos,detector["r"]):
            alive = False
            break

        # Rotate Angle and Polarize due to Glucose
        theta_glucose = get_material("alpha") * gc * (s/10)
        D = rotation_matrix_glucose(theta_glucose)
        stokes = D @ stokes

        #If detected, we record. If not detected, we scatter and travel again
        if detected_photon:

            return alive, step_counter, total_path_length, stokes

        # Scatter and Polarize due to scattering particle
        dir1= change_direction(g=0.9)
        phi= compute_phi(dir, dir1)
        dir=dir1
        M= mie_scattering_matrix_rayleigh(phi)
        D= rotation_matrix_phi(phi) @ M @ rotation_matrix_phi(phi)
        stokes = D @ stokes

    #If we reach here, alive = false, so the rest of return doesn't matter
    return alive, step_counter, total_path_length, stokes

