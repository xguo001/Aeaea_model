"""
Visualization of simulation results, generating plots analogous to Figures 3–5 from the paper.
"""
import matplotlib.pyplot as plt

def plot_P1_DoP_vs_GC(GC_values, P1_curves, DoP_curves, radii):
    """
    Plot P1 and DoP vs glucose concentration for multiple particle radii.
    radii: list of particle radii (in micrometers).
    P1_curves: dict mapping radius -> array of P1 values (same length as GC_values).
    DoP_curves: dict mapping radius -> array of DoP values.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10,5))
    ax1, ax2 = axes
    for r in radii:
        label = f"r = {r:.3f} μm"
        P1_vals = P1_curves[r]
        DoP_vals = DoP_curves[r]
        ax1.plot(GC_values, P1_vals, marker='o', label=label)
        ax2.plot(GC_values, DoP_vals, marker='s', label=label)
    ax1.set_xlabel("Glucose Concentration (g/dL)")
    ax1.set_ylabel("P1 (Index of Polarimetric Purity)")
    ax1.set_title("P1 vs Glucose Concentration")
    ax1.legend()
    ax1.grid(True)
    ax2.set_xlabel("Glucose Concentration (g/dL)")
    ax2.set_ylabel("Degree of Polarization (DoP)")
    ax2.set_title("DoP of Emitted Light vs Glucose Concentration")
    ax2.legend()
    ax2.grid(True)
    fig.suptitle("Forward Scattering Results for Different Particle Radii")
    plt.tight_layout()
    return fig













