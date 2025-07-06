"""
Polarization utility functions: Mueller matrices, rotation matrices, and conversions between polarization representations.
"""
import numpy as np
import math

def mueller_rotation(angle_deg):
    """
    Construct Mueller rotation matrix for rotation of polarization plane by `angle_deg` about propagation axis.
    This corresponds to rotating the Stokes Q-U parameters by 2*angle.
    """
    a = math.radians(angle_deg)
    cos2a = math.cos(2*a)
    sin2a = math.sin(2*a)
    M = np.array([[1,      0,       0,      0],
                  [0,   cos2a,   sin2a,     0],
                  [0,  -sin2a,   cos2a,     0],
                  [0,      0,       0,      1]], dtype=float)
    return M

def rotate_stokes(S, angle_deg):
    """
    Rotate Stokes vector `S` by `angle_deg` (degrees) about the propagation axis.
    Simulates optical rotation of the polarization plane by the given angle.
    """
    a = math.radians(angle_deg)
    cos2a = math.cos(2*a)
    sin2a = math.sin(2*a)
    I, Q, U, V = S
    Q_new = Q*cos2a + U*sin2a
    U_new = -Q*sin2a + U*cos2a
    # I and V remain the same
    return np.array([I, Q_new, U_new, V], dtype=float)

def jones_to_stokes(jones_vec):
    """
    Convert a Jones vector (Ex, Ey) to a normalized Stokes vector [I, Q, U, V].
    Assumes the Jones vector is normalized such that total intensity I=1.
    """
    Ex, Ey = jones_vec
    I = 1.0
    Q = np.abs(Ex)**2 - np.abs(Ey)**2
    U = 2 * np.real(Ex * np.conjugate(Ey))
    V = 2 * np.imag(np.conjugate(Ex) * Ey)
    return np.array([I, Q, U, V], dtype=float)

def stokes_to_jones(S):
    """
    Estimate a Jones vector for a (fully polarized) Stokes vector S.
    Assumes S is fully polarized (Degree of Polarization = 1).
    """
    I, Q, U, V = S
    if I == 0:
        return np.array([0+0j, 0+0j])
    # Degree of polarization
    dop = math.sqrt(Q**2 + U**2 + V**2) / I
    if dop < 1e-6:
        # Unpolarized: return arbitrary Jones (e.g., 45° polarization)
        return np.array([complex(math.sqrt(I/2)), complex(math.sqrt(I/2))])
    # Assume V ≈ 0 for simplicity (linear polarization).
    Ex_amp = math.sqrt((I + Q) / 2.0)
    Ey_amp = math.sqrt((I - Q) / 2.0)
    # Phase adjustment for U
    if U < 0:
        Ey_amp = -Ey_amp
    Ex = complex(Ex_amp)
    Ey = complex(Ey_amp)
    # If V is nonzero, introduce a phase difference (approximate)
    if abs(V) > 1e-6:
        sin_delta = V / (2*Ex_amp*Ey_amp) if Ex_amp*Ey_amp != 0 else 0
        sin_delta = max(min(sin_delta, 1.0), -1.0)
        delta = math.asin(sin_delta)
        Ey = Ey * complex(math.cos(delta), math.sin(delta))
    return np.array([Ex, Ey], dtype=complex)













