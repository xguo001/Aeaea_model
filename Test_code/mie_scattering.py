"""
Mie scattering functions: computing scattering phase function and sampling scattering events.

This module may use the `miepython` library for accurate Mie scattering calculations.
If `miepython` is not installed, an approximate scattering behavior is used.
"""
import numpy as np
import math
try:
    import miepython
except ImportError:
    miepython = None

# Cache for scattering phase function data to avoid redundant calculations
_phase_cache = {}

def _compute_phase_function(radius, sphere_index, medium_index, wavelength):
    """
    Compute unpolarized scattering phase function (normalized differential scattering) for a sphere.
    Returns arrays (mu_vals, cdf_vals) where mu = cos(theta) and cdf_vals is the cumulative distribution function.
    """
    if miepython:
        # Use miepython to compute scattering intensities vs angle
        x = 2 * math.pi * radius / wavelength  # size parameter
        theta_vals = np.linspace(0, math.pi, 361)
        mu_vals = np.cos(theta_vals)
        m_rel = sphere_index / medium_index
        # Scattered intensity in planes parallel and perpendicular to incident field
        i_par = miepython.i_par(m_rel, x, mu_vals)
        i_per = miepython.i_per(m_rel, x, mu_vals)
        i_unpol = 0.5 * (i_par + i_per)
    else:
        # If miepython not available, approximate with Rayleigh phase function ~ (1 + cos^2θ)
        mu_vals = np.linspace(-1, 1, 201)
        i_unpol = 0.75 * (1 + mu_vals**2)  # not exact normalization but acceptable shape
    # Normalize to get PDF and CDF
    pdf = i_unpol / np.trapz(i_unpol, mu_vals)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    return mu_vals, cdf

def sample_scattering_event(dir_in, S_in, radius, sphere_index, medium_index, wavelength):
    """
    Sample a scattering event for a photon with incoming direction and Stokes state.
    Returns the new direction (unit vector) and the updated Stokes vector after scattering.

    This uses an approximate polarization scattering model:
    The scattering angle is sampled from the unpolarized phase function. The polarization state is partially depolarized 
    based on the scattering angle (larger angles cause more depolarization), and the polarization plane is rotated by the random azimuthal angle.
    For a more accurate simulation, a full Mie Mueller matrix should be applied [oai_citation:3‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=pure%20elements%20,as%20follows%3A%20H%20%3D%201) [oai_citation:4‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=The%20eigenvalue%20spectrum%20of%20the,29).
    """
    # Normalize input direction
    k_in = dir_in / np.linalg.norm(dir_in)
    # Retrieve or compute phase function for this radius
    key = (radius, sphere_index, medium_index, wavelength)
    if key not in _phase_cache:
        mu_vals, cdf_vals = _compute_phase_function(radius, sphere_index, medium_index, wavelength)
        _phase_cache[key] = (mu_vals, cdf_vals)
    else:
        mu_vals, cdf_vals = _phase_cache[key]
    # Sample scattering cosine (mu = cos θ) via inverse CDF
    r = np.random.random()
    idx = np.searchsorted(cdf_vals, r)
    if idx >= len(mu_vals):
        idx = len(mu_vals) - 1
    mu = mu_vals[idx]
    theta = math.acos(mu)
    # Sample random azimuth φ in [0, 2π)
    phi = 2 * math.pi * np.random.random()
    # Construct orthonormal basis (u, v, k_in) for polarization reference
    if abs(k_in[2]) < 0.999:
        temp = np.array([0.0, 0.0, 1.0])
    else:
        temp = np.array([1.0, 0.0, 0.0])
    u = np.cross(temp, k_in)
    if np.linalg.norm(u) < 1e-8:
        u = np.array([1.0, 0.0, 0.0])
    u = u / np.linalg.norm(u)
    v = np.cross(k_in, u)
    v = v / np.linalg.norm(v)
    # New direction unit vector
    new_dir = mu * k_in + math.sqrt(max(0.0, 1 - mu**2)) * (math.cos(phi) * u + math.sin(phi) * v)
    new_dir = new_dir / np.linalg.norm(new_dir)
    # Approximate polarization change:
    # Depolarization factor p = cos^2(theta) (retain polarization for small angles, reduce for large angles)
    cos_angle = np.dot(new_dir, k_in)
    p = cos_angle**2
    I, Q, U, V = S_in
    Q *= p; U *= p; V *= p
    S_temp = np.array([I, Q, U, V], dtype=float)
    # Rotate polarization plane by azimuth φ around k_in axis [oai_citation:5‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=the%20photon%20after%20the%20collision,%CF%86%29%20k%20%E2%88%8F%EF%B8%82%E2%80%B2) [oai_citation:6‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=Rk,%CF%86%29%20k%20%E2%88%8F%EF%B8%82%E2%80%B2%20k%3D0)
    angle_deg = math.degrees(phi)
    cos2a = math.cos(2 * math.radians(angle_deg))
    sin2a = math.sin(2 * math.radians(angle_deg))
    Qr = S_temp[1]*cos2a + S_temp[2]*sin2a
    Ur = -S_temp[1]*sin2a + S_temp[2]*cos2a
    S_out = np.array([S_temp[0], Qr, Ur, S_temp[3]], dtype=float)
    return new_dir, S_out













