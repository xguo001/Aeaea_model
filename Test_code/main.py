"""
Main orchestration script for running simulations, performing analysis, and generating plots.
"""
import numpy as np
from mc_simulation import simulate_forward
from mueller_analysis import assemble_mueller_matrix, hermitian_covariance_matrix, compute_IPPs
from fitting import fit_cubic_P1, fit_quadratic, invert_P1_to_concentration
from visualization import plot_P1_DoP_vs_GC
import matplotlib.pyplot as plt

def main():
    # Particle radii (μm) and glucose concentration range (g/dL) [oai_citation:11‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=which%20the%20size%20of%20scattering,From)
    radii = [0.05, 0.125, 0.5]
    GC_values = np.linspace(0, 400, 5)  # e.g., 0,100,200,300,400 g/dL
    # Simulation parameters
    num_photons = 10000  # photons per input state (increase for smoother results)
    medium_index = 1.33
    sphere_index = 1.57
    thickness = 1.0  # cm slab thickness
    scat_coeff = 1.0  # scattering coefficient (cm^-1) [oai_citation:12‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=but%20the%20RI%20,pdf%20page%206%20of%2013)
    abs_coeff = 0.01  # absorption coefficient (cm^-1) [oai_citation:13‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=but%20the%20RI%20,pdf%20page%206%20of%2013)
    P1_results = {r: [] for r in radii}
    DoP_results = {r: [] for r in radii}
    for r in radii:
        for GC in GC_values:
            det = {'thickness': thickness}
            # Run simulation for multiple input polarizations
            outputs = {}
            S_out_H = simulate_forward(det, [1,1,0,0], num_photons, r, GC/100.0, scat_coeff, abs_coeff, medium_index, sphere_index, wavelength=0.589)
            outputs['H'] = S_out_H
            S_out_V = simulate_forward(det, [1,-1,0,0], num_photons, r, GC/100.0, scat_coeff, abs_coeff, medium_index, sphere_index, wavelength=0.589)
            outputs['V'] = S_out_V
            S_out_P45 = simulate_forward(det, [1,0,1,0], num_photons, r, GC/100.0, scat_coeff, abs_coeff, medium_index, sphere_index, wavelength=0.589)
            outputs['P45'] = S_out_P45
            S_out_R = simulate_forward(det, [1,0,0,1], num_photons, r, GC/100.0, scat_coeff, abs_coeff, medium_index, sphere_index, wavelength=0.589)
            outputs['R'] = S_out_R
            # Compute P1 from IPPs
            M = assemble_mueller_matrix(outputs)
            H = hermitian_covariance_matrix(M)
            P1, P2, P3 = compute_IPPs(H)
            P1_results[r].append(P1)
            # Compute DoP of emitted light for unpolarized input (average of H and V outputs) [oai_citation:14‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=GCs%20%28Fig,However) [oai_citation:15‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=Fig,we%20hope%20to%20obtain%20the)
            S_unpol = 0.5 * (S_out_H + S_out_V)
            I, Q, U, V = S_unpol
            DoP = 0.0
            if I > 0:
                DoP = np.sqrt(Q**2 + U**2 + V**2) / I
            DoP_results[r].append(DoP)
    # Convert to numpy arrays
    for r in radii:
        P1_results[r] = np.array(P1_results[r])
        DoP_results[r] = np.array(DoP_results[r])
    # Fit cubic P1(GC) for each radius and quadratic dependence of coefficients on radius [oai_citation:16‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=the%20data%20tracing%20points%20in,obtain%20the%20fitting%20function%20of) [oai_citation:17‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=P1%20%3D%20%E2%88%920,cubic%20polynomial%20is%20also%20fitted) [oai_citation:18‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=a%20%3D%20%E2%88%920,3895%20%2814)
    cubic_coeffs = {}
    for r in radii:
        a, b, c, d = fit_cubic_P1(GC_values, P1_results[r])
        cubic_coeffs[r] = (a, b, c, d)
        print(f"Radius {r} μm: P1(GC) ≈ {a:.4f}*x^3 + {b:.4f}*x^2 + {c:.4f}*x + {d:.4f}")
    r_arr = np.array(radii)
    a_vals = np.array([cubic_coeffs[r][0] for r in radii])
    b_vals = np.array([cubic_coeffs[r][1] for r in radii])
    c_vals = np.array([cubic_coeffs[r][2] for r in radii])
    d_vals = np.array([cubic_coeffs[r][3] for r in radii])
    A_a, B_a, C_a = fit_quadratic(r_arr, a_vals)
    A_b, B_b, C_b = fit_quadratic(r_arr, b_vals)
    A_c, B_c, C_c = fit_quadratic(r_arr, c_vals)
    A_d, B_d, C_d = fit_quadratic(r_arr, d_vals)
    print(f"Fitted a(r) = {A_a:.4f}*r^2 + {B_a:.4f}*r + {C_a:.4f}")
    print(f"Fitted b(r) = {A_b:.4f}*r^2 + {B_b:.4f}*r + {C_b:.4f}")
    print(f"Fitted c(r) = {A_c:.4f}*r^2 + {B_c:.4f}*r + {C_c:.4f}")
    print(f"Fitted d(r) = {A_d:.4f}*r^2 + {B_d:.4f}*r + {C_d:.4f}")
    # Create polynomial functions of r
    a_of_r = np.poly1d([A_a, B_a, C_a])
    b_of_r = np.poly1d([A_b, B_b, C_b])
    c_of_r = np.poly1d([A_c, B_c, C_c])
    d_of_r = np.poly1d([A_d, B_d, C_d])
    coeff_funcs = {'a': a_of_r, 'b': b_of_r, 'c': c_of_r, 'd': d_of_r}
    # Example inversion: given P1 at radius 0.05 μm and known radius, retrieve GC [oai_citation:19‡file-9cuejnth997c53y3ssoxyz](file://file-9cuEJntH997c53y3sSoXYZ#:~:text=values%20on%20each%20interval%20were,the%20bigger%20error%20occurs%20in)
    test_radius = radii[0]
    test_index = 2  # e.g., third point (around 200 g/dL)
    measured_P1 = P1_results[test_radius][test_index]
    est_GC = invert_P1_to_concentration(measured_P1, test_radius, coeff_funcs)
    print(f"Inversion test for r={test_radius} μm, P1≈{measured_P1:.4f}: retrieved GC ≈ {est_GC:.1f} g/dL (true GC≈{GC_values[test_index]} g/dL)")
    # Plot results
    fig = plot_P1_DoP_vs_GC(GC_values, P1_results, DoP_results, radii)
    plt.show()

if __name__ == "__main__":
    main()













