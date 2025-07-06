import numpy as np
import miepython
from numba.np.arraymath import np_histogram

# --- Constants ---
MU_S = 10  # Scattering coefficient [1/cm]
MU_A = 1   # Absorption coefficient [1/cm]
G = 0.9    # Anisotropy
ALPHA = 0.01  # Optical rotation per concentration unit [rad/(cm*mmol/L)]
GLUCOSE_CONC = 5.0  # mmol/L
medium_index=1.33
sphere_index=1.57
wavelength=1.55
radius=0.05
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

def sample_mie_direction(radius, wavelength, n_particle, n_medium, n_angles=2000):
    """
    Sample a scattering direction based on Mie theory using unpolarized light intensity.
    """
    m = n_particle / n_medium                    # Relative refractive index
    x = 2 * np.pi * radius / wavelength          # Size parameter
    mu = np.linspace(-1, 1, n_angles)            # cos(theta)

    # Get unpolarized intensity from Mie theory
    intensity = miepython.i_unpolarized(m, x, mu)

    # Normalize to probability density
    pdf = intensity / np.trapezoid(intensity, mu)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]

    # Inverse transform sampling to get cos(theta)
    rnd = np.random.rand()
    idx = np.searchsorted(cdf, rnd)
    cos_theta = mu[idx]
    sin_theta = np.sqrt(1 - cos_theta**2)

    # Sample azimuthal angle phi
    phi = 2 * np.pi * np.random.rand()

    # Convert to Cartesian direction
    dx = sin_theta * np.cos(phi)
    dy = sin_theta * np.sin(phi)
    dz = cos_theta

    return np.array([dx, dy, dz])


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
        direction = sample_mie_direction(radius, wavelength, n_particle=1.57, n_medium=1.33, n_angles=2000)

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

def simulate_multiple_photons(n_photons):
    """ Run the simulation for multiple photons and collect results """
    results = []
    for _ in range(n_photons):
        result = simulate_photon()
        results.append(result)
    return results




# --- Run Simulation ---
if __name__ == "__main__":
    n_photons = 10
    results = simulate_multiple_photons(n_photons)

    for i, result in enumerate(results):
        print(f"\nPhoton {i + 1}:")
        for key, value in result.items():
            print(f"  {key}: {value}")

