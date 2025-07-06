import numpy as np

# --- Constants ---
MU_S = 10  # Scattering coefficient [1/cm]
MU_A = 1   # Absorption coefficient [1/cm]
G = 0.9    # Anisotropy
ALPHA = 0.01  # Optical rotation per concentration unit [rad/(cm*mmol/L)]
GLUCOSE_CONC = 5.0  # mmol/L

# --- Helper Functions ---

def sample_step(mu_t):
    """ Sample a step size from exponential distribution """
    return -np.log(np.random.rand()) / mu_t

def sample_scatter_direction(g):
    """ Sample scattering direction using Henyey-Greenstein phase function """
    cos_theta = (1/(2*g)) * (1 + g**2 - ((1 - g**2)/(1 - g + 2*g*np.random.rand()))**2)
    sin_theta = np.sqrt(1 - cos_theta**2)
    phi = 2 * np.pi * np.random.rand()
    dx = sin_theta * np.cos(phi)
    dy = sin_theta * np.sin(phi)
    dz = cos_theta
    return np.array([dx, dy, dz])

def rotate_stokes_vector(stokes, theta):
    """ Apply optical rotation to Stokes vector (rotate Q and U) """
    I, Q, U, V = stokes
    Q_rot = Q * np.cos(2*theta) - U * np.sin(2*theta)
    U_rot = Q * np.sin(2*theta) + U * np.cos(2*theta)
    return np.array([I, Q_rot, U_rot, V])

# --- Photon Simulation Function ---

def simulate_photon():
    position = np.array([0.0, 0.0, 0.0])  # Start at origin
    direction = np.array([0.0, 0.0, 1.0])  # Initial direction
    weight = 1.0
    stokes = np.array([1.0, 1.0, 0.0, 0.0])  # Linear polarized light along x

    mu_t = MU_S + MU_A
    path_length = 0.0

    for _ in range(1000):
        s = sample_step(mu_t)
        position += s * direction
        path_length += s
        direction = sample_scatter_direction(G)

        # Optical rotation by glucose (alpha * concentration * path_length)
        theta = ALPHA * GLUCOSE_CONC * s
        stokes = rotate_stokes_vector(stokes, theta)

        # Absorb a bit of energy
        weight *= np.exp(-MU_A * s)

        if weight < 1e-4:
            break  # Photon terminated

    return {
        "final_position": position,
        "final_direction": direction,
        "final_stokes": stokes,
        "total_path_length": path_length
    }

# --- Run Simulation ---
if __name__ == "__main__":
    result = simulate_photon()
    for key, value in result.items():
        print(f"{key}: {value}")

