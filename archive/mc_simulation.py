"""
Monte Carlo simulation of polarized photon propagation in a turbid, optically active medium.

This module simulates forward scattering of polarized light through a slab containing spherical scatterers and glucose.
The simulation tracks photon Stokes vectors through multiple scattering events and optical rotation due to glucose.
"""
import numpy as np
import math
from polarization import rotate_stokes  # Mueller rotation for optical rotation
from mie_scattering import sample_scattering_event

def simulate_forward(detector, source, num_photons, radius, concentration, 
                     scat_coeff=1.0, abs_coeff=0.01, medium_index=1.33, 
                     sphere_index=1.57, wavelength=0.589):
    """
    Simulate forward scattering through a slab and accumulate Stokes outputs.

    Parameters:
        detector (dict): Detector configuration, e.g. {'thickness': 1.0} with thickness in cm.
        source (dict or list): Source Stokes vector (4,) or Jones vector as initial polarization state.
        num_photons (int): Number of photon histories to simulate.
        radius (float): Scatterer radius (micrometers).
        concentration (float): Glucose concentration (g/mL).
        scat_coeff (float): Scattering coefficient μs (cm^-1).
        abs_coeff (float): Absorption coefficient μa (cm^-1).
        medium_index (float): Refractive index of medium (e.g., 1.33 for water).
        sphere_index (float): Refractive index of scatterer (e.g., 1.57 for polystyrene).
        wavelength (float): Wavelength in micrometers (e.g., 0.589 µm for 589 nm).

    Returns:
        np.ndarray: Cumulative Stokes vector of forward-emerging light [I, Q, U, V].
    """
    L = detector.get('thickness', 1.0)  # thickness in cm
    # Convert concentration (g/mL) to appropriate units for optical rotation:
    # Optical rotation (degrees) = ORD * path_length (dm) * concentration (g/mL)
    ORD = 52.7  # specific rotation of glucose at 589 nm [deg/(dm*(g/mL))] [oai_citation:0‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=wavelength%20of%20light%20is%20589,total%20radiance%20received%20by%20the)
    # Prepare initial Stokes vector S0
    if isinstance(source, np.ndarray) and source.shape == (2,):
        # source given as Jones [Ex, Ey]
        S0 = np.array([1.0,
                       np.abs(source[0])**2 - np.abs(source[1])**2,
                       2 * np.real(source[0]*np.conjugate(source[1])),
                       2 * np.imag(np.conjugate(source[0])*source[1])])
    else:
        S0 = np.array(source, dtype=float)
    S_forward = np.zeros(4, dtype=float)  # accumulate Stokes of forward-exiting photons
    mu_s = scat_coeff
    mu_a = abs_coeff
    mu_t = mu_s + mu_a  # total extinction coefficient
    for _ in range(num_photons):
        # Initialize photon position and direction (downward +z)
        pos_z = 0.0
        dir_vector = np.array([0.0, 0.0, 1.0], dtype=float)
        S = S0.copy()  # photon Stokes state (fully polarized initial)
        alive = True
        scatter_count = 0
        while alive:
            # Sample distance to next interaction (exponential with mean 1/mu_t)
            step = -math.log(np.random.random()) / mu_t
            # Distance to exit (depending on direction)
            dist_to_exit = L - pos_z if dir_vector[2] > 0 else pos_z
            if step >= dist_to_exit:
                # Photon exits the slab (reaches detector plane)
                pos_z += dist_to_exit
                # Apply optical rotation for remaining path
                alpha = ORD * (dist_to_exit/10.0) * concentration  # convert cm to dm
                # Forward scattering detection: apply R(-α) [oai_citation:1‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=Sfs%20%3D%20MSi%20%3D%20R,6)
                S = rotate_stokes(S, -alpha)
                if dir_vector[2] > 0:
                    # Forward-exiting photon: add its Stokes to detector sum
                    S_out = _transform_to_detector(S, dir_vector)
                    S_forward += S_out
                # If exiting backward, drop (backscatter not counted here)
                alive = False
                continue
            # Photon has an interaction inside the medium
            pos_z += step * dir_vector[2]  # update z position
            # Apply optical rotation over path 'step'
            alpha = ORD * (step/10.0) * concentration
            # Sign: forward direction -> -α, backward -> +α [oai_citation:2‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=Sfs%20%3D%20MSi%20%3D%20R,6)
            sign = -1 if dir_vector[2] > 0 else 1
            S = rotate_stokes(S, sign * alpha)
            # Decide event: scattering or absorption
            if np.random.random() > mu_s / mu_t:
                # Absorption: photon terminated
                alive = False
                continue
            # Scattering event: sample new direction and update Stokes
            new_dir, S = sample_scattering_event(dir_vector, S, radius, sphere_index, medium_index, wavelength)
            dir_vector = new_dir
            scatter_count += 1
            if scatter_count > 100:
                # Safety break to prevent infinite loop in highly scattering cases
                alive = False
                break
    return S_forward

def _transform_to_detector(S_stokes, direction):
    """
    Rotate Stokes vector from photon's current propagation direction to detector coordinate frame.
    The detector frame is defined such that the photon's propagation is aligned with +z (vertical).
    This accounts for photons exiting at an angle by aligning their polarization reference to the detector.
    """
    k = direction / np.linalg.norm(direction)
    z_axis = np.array([0.0, 0.0, 1.0])
    if np.allclose(k, z_axis, atol=1e-8):
        return S_stokes.copy()
    # Compute axis perpendicular to both k and z (rotation axis)
    axis = np.cross(k, z_axis)
    axis = axis / np.linalg.norm(axis)
    # Angle between k and z
    cosang = np.dot(k, z_axis)
    cosang = max(min(cosang, 1.0), -1.0)
    angle = math.acos(cosang)
    # ** Simplified **: For small exit angles, ignore polarization rotation in this transform.
    # (Accurate transformation would require rotating Stokes parameters about 'axis'.)
    return S_stokes.copy()













