import numpy as np
import matplotlib.pyplot as plt
import gc_ar.results as results
import gc_ar.set_parameters as set_parameters
from mpl_toolkits.mplot3d import Axes3D

def plot_photon_paths(photon_paths):
    """Plots all photon paths in one 3D figure"""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for path in photon_paths:
        ax.plot(path[:, 0], path[:, 1], path[:, 2], alpha=0.7)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'{len(photon_paths)} Photon Paths')
    plt.tight_layout()
    plt.show()

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
    plt.plot(results.return_variable_vs_output()[:,0],results.return_variable_vs_output()[:,1],"-o")
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
    total_absorbed = results.return_absorption_matrix()[:, 0].sum()
    total_detected = np.sum(results.return_detected_energy())
    total_out_of_bound = np.sum(results.return_out_of_bound_energy())
    
    # Create energy summary text
    energy_text = [
        "-------------------",
        f"total absorbed energy: {total_absorbed:.2f}",
        f"total detected energy: {total_detected:.2f}",
        f"total out of bound energy: {total_out_of_bound:.2f}"
    ]
    
    # Add whitespace below the plot and place parameters text outside plot area
    plt.subplots_adjust(bottom=0.5)
    plt.figtext(0.1, 0.01, '\n'.join(left_column + energy_text), fontsize=10, verticalalignment='bottom')
    plt.figtext(0.5, 0.01, '\n'.join(right_column), fontsize=10, verticalalignment='bottom')
    
    plt.show()
