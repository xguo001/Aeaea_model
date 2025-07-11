import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

# -----------------------------
# TAKE A GLUCOSE LEVEL AND PLOT ALL ROTATION ANGLES INTO A HISTOGRAM
# -----------------------------
def plot_photon_histogram(results):
    data = np.array(results)
    sns.histplot(data, kde=True, stat="density", label="Data", bins=30)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, np.mean(data), np.std(data))
    plt.plot(x, p, 'k', linewidth=2, label="Gaussian fit")
    plt.legend()
    plt.title("Histogram with Gaussian Fit")
    plt.show()

# -----------------------------
# TAKE MULTIPLE GC LEVELS AND PLOT AGAINST THE AVERAGE ROTATION ANGLES
# -----------------------------
def plot_GC_vs_angles_plot(GC_a, angles):

    plt.plot(GC_a, angles)
    plt.xlabel("Glucose concentration")
    plt.ylabel("Rotation Angle")
    plt.title("GC vs. Angle")
    plt.show()
#    plt.savefig(r"/home/ubuntu/results/gc_angles.png",dpi=300)

def plot_GC_vs_angles_write(GC_a, angles):
    plt.plot(GC_a, angles)
    plt.xlabel("Glucose concentration")
    plt.ylabel("Rotation Angle")
    plt.title("GC vs. Angle")
    plt.savefig(r"/home/ubuntu/results/gc_angles.png",dpi=300)
