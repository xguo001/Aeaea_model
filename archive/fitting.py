"""
Curve fitting for P1 vs glucose concentration and coefficients vs particle radius, including inversion to retrieve concentration.
"""
import numpy as np

def fit_cubic_P1(GC_values, P1_values):
    """
    Fit a cubic polynomial P1(x) = a*x^3 + b*x^2 + c*x + d to the given data.
    Returns coefficients (a, b, c, d).
    """
    coeffs = np.polyfit(GC_values, P1_values, 3)
    a, b, c, d = coeffs
    return a, b, c, d

def fit_quadratic(x_vals, y_vals):
    """
    Fit a quadratic polynomial y(r) = A*r^2 + B*r + C.
    Returns coefficients (A, B, C).
    """
    coeffs = np.polyfit(x_vals, y_vals, 2)
    A, B, C = coeffs
    return A, B, C

def invert_P1_to_concentration(P1_meas, radius, coeff_funcs):
    """
    Invert the fitted polynomial model to retrieve glucose concentration from a measured P1 and known particle radius.
    `coeff_funcs` is a dict with keys 'a','b','c','d' that return polynomial coefficients as functions of radius.
    Returns the estimated glucose concentration (in the same units as the fit data).
    """
    # Coefficients at given radius
    a = coeff_funcs['a'](radius)
    b = coeff_funcs['b'](radius)
    c = coeff_funcs['c'](radius)
    d = coeff_funcs['d'](radius)
    # Solve a*x^3 + b*x^2 + c*x + (d - P1_meas) = 0 for x
    coeffs = [a, b, c, d - P1_meas]
    roots = np.roots(coeffs)
    real_roots = [r.real for r in roots if abs(r.imag) < 1e-6]
    if not real_roots:
        return None
    # Filter for non-negative roots (glucose concentration can't be negative)
    candidates = [r for r in real_roots if r >= 0]
    if not candidates:
        return real_roots[0]
    # Choose the smallest positive root as the physical solution
    return min(candidates)













