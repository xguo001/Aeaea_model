import numpy as np
import matplotlib.pyplot as plt
import initialize.results as results
import initialize.set_parameters as set_parameters
from matplotlib import cm
from matplotlib.colors import Normalize
from detectors_and_plots.monitors import return_total_absoprtions
import initialize.results as results

def plot_energy_distribution(absorption_points):
    """Plots a 3D scatter of photon energy deposition"""
    data = np.array(absorption_points)  # shape [n_points, 4]

    energies = data[:, 0]
    xs = data[:, 1]
    ys = data[:, 2]
    zs = data[:, 3]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(xs, ys, zs, c=energies, cmap='viridis', s=10, alpha=0.8)
    plt.colorbar(sc, label='Energy')

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Photon Energy Deposition Map")
    plt.tight_layout()
    plt.show()



def plot_variable_vs_angle(name):
    #variable, np.mean(angular rotation), np.mean(step_counters), np.mean(path_lengths_collector), death_counters
    plt.plot(results.return_variable_vs_output()[:, 0], results.return_variable_vs_output()[:, 1], "-o")
    plt.ylabel("angles")
    plt.xlabel(name + " vs. angles")
    
    # Get all parameters and format them for display
    params = set_parameters.parameters
    param_text = []
    for key, value in params.items():
        if isinstance(value, np.ndarray):
            value_str = f"[{', '.join([f'{x:.2f}' if isinstance(x, (int, float)) else str(x) for x in value])}]"
        elif isinstance(value, (int, float)):
            value_str = f"{value:.2f}"
        else:
            value_str = str(value)
        param_text.append(f"{key}: {value_str}")
    
    # Split parameters into two columns
    mid_point = len(param_text) // 2
    left_column = param_text[:mid_point]
    right_column = param_text[mid_point:]
    
    # Calculate energy values
    total_absorption = return_total_absoprtions()

    # Create energy summary text
    energy_text = [
        "-------------------",
        f"total absorbed energy: {total_absorption[0]:.2f}",
        f"total detected energy: {total_absorption[1]:.2f}",
        f"total out of bound energy: {total_absorption[2]:.2f}"
    ]
    
    # Add whitespace below the plot and place parameters text outside plot area
    plt.subplots_adjust(bottom=0.5)
    plt.figtext(0.1, 0.01, '\n'.join(left_column + energy_text), fontsize=10, verticalalignment='bottom')
    plt.figtext(0.5, 0.01, '\n'.join(right_column), fontsize=10, verticalalignment='bottom')
    
    plt.show()

def spherical_to_cartesian(r, theta, phi=0):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x, y, z])

def plot_photon_paths(photon_paths, detector=None, sphere_radius=None):
    """
    Visualizes photon paths with a cone-based detector and optional sphere.

    Args:
        photon_paths: list of np.ndarray [N, 3] — list of photon paths
        detector: dict with keys 'cone_axis', 'cone_center', 'r'
        sphere_radius: float, radius of background sphere
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    cmap = cm.get_cmap('bone_r')  # darker = higher energy
    # Plot photon paths
    for path in photon_paths:
        path = np.atleast_2d(path)
        energies = path[:, 0]

        # Normalize energies for color mapping
        norm = Normalize(vmin=np.min(energies), vmax=np.max(energies))
        colors = cmap(norm(energies))

        # Plot each segment with its energy color
        for i in range(len(path) - 1):
            x = path[i:i + 2, 1]
            y = path[i:i + 2, 2]
            z = path[i:i + 2, 3]
            ax.plot(x, y, z, color=colors[i], linewidth=1.5)
        # Draw a circle where photon energy drops below threshold
        energy_threshold = 1e-4
        below_threshold = np.where(energies < energy_threshold)[0]
        if len(below_threshold) > 0:
            idx = below_threshold[0]
            ax.scatter(
                path[idx, 1], path[idx, 2], path[idx, 3],
                color='black', s=500, marker='o',
                label='Absorbed' if 'Absorbed' not in ax.get_legend_handles_labels()[1] else ""
            )

    # Plot sphere if needed
    if sphere_radius is not None:
        u = np.linspace(0, 2 * np.pi, 30)
        v = np.linspace(0, np.pi, 30)
        x = sphere_radius * np.outer(np.cos(u), np.sin(v))
        y = sphere_radius * np.outer(np.sin(u), np.sin(v))
        z = sphere_radius * np.outer(np.ones_like(u), np.cos(v))
        ax.plot_surface(x, y, z, color='lightblue', alpha=0.1, edgecolor='none')

    # Plot detector cone origin and direction
    if detector is not None:
        axis = detector["cone_axis"] / np.linalg.norm(detector["cone_axis"])
        r = detector["r"]
        origin = r * axis

        ax.scatter(*origin, color='red', s=100, label='Detector')

        # Draw axis line
        ax.quiver(*origin, *(0.2 * axis), color='red', arrow_length_ratio=0.1, linewidth=2)

        # Optionally: plot acceptance cone edges
        # Generate vectors perpendicular to axis
        perp1 = np.cross(axis, [1, 0, 0]) if abs(axis[0]) < 0.9 else np.cross(axis, [0, 1, 0])
        perp1 /= np.linalg.norm(perp1)
        perp2 = np.cross(axis, perp1)

        # Circle at cone opening
        theta = detector["cone_center"]
        radius = np.tan(theta) * 0.1  # visual size
        circle_pts = [
            origin + 0.15 * axis + radius * (np.cos(a) * perp1 + np.sin(a) * perp2)
            for a in np.linspace(0, 2 * np.pi, 50)
        ]
        circle_pts = np.array(circle_pts)
        ax.plot(circle_pts[:, 0], circle_pts[:, 1], circle_pts[:, 2], color='red', alpha=0.5)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Photon Paths with Directional Detector")
    ax.legend()
    plt.tight_layout()
    limit = detector["r"] + 0.1
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    plt.show()
    limit = detector["r"] + 0.1
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)


def plot_absorbed_energy_vs_time(n_bins):
    """Plots absorbed energy vs time using absorption matrix data with equal-sized bins"""
    absorption_data = results.return_absorption_matrix()

    if len(absorption_data) == 0:
        return

    energies = absorption_data[:, 0]  # First column: energy
    times = absorption_data[:, 4]  # Last column: global_time

    # Sort data by time to create equal-sized bins
    # Time should already be sorted -- this is written by Claude lol
    sorted_indices = np.argsort(times)
    sorted_times = times[sorted_indices]
    sorted_energies = energies[sorted_indices]

    # Calculate points per bin
    n_points = len(sorted_times)
    points_per_bin = n_points // n_bins

    # Create bins with equal number of points
    bin_centers = []
    binned_energies = []

    for i in range(n_bins):
        start_idx = i * points_per_bin
        if i == n_bins - 1:  # Last bin gets remaining points
            end_idx = n_points
        else:
            end_idx = (i + 1) * points_per_bin

        # Calculate bin center as median time in this bin
        bin_times = sorted_times[start_idx:end_idx]
        bin_center = np.median(bin_times)
        bin_centers.append(bin_center)

        # Sum energies in this bin
        bin_energy = np.sum(sorted_energies[start_idx:end_idx])
        binned_energies.append(bin_energy)

    plt.figure()
    plt.plot(bin_centers, binned_energies, 'o-', markersize=4, linewidth=1)
    plt.xlabel('Time')
    plt.ylabel('Absorbed Energy (Equal-sized bins)')
    plt.title('Absorbed Energy vs Time (Equal-sized bins)')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def plot_input():

    return None