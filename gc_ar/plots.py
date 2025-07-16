import numpy as np
import matplotlib.pyplot as plt
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