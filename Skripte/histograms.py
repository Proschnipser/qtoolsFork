"""
The dense_1 (output) layer kernel has shape (1891, 5): 1891 incoming
weights feeding into each of the 5 output neurons. This plots one
panel per output neuron (5 total), each showing its own 1891 weights
with x = weight index/identifier (0..1890) and y = weight value.

Usage:
    python dense1_per_neuron_histograms.py path/to/weights.h5
"""

import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_dense1_kernel(filepath):
    with h5py.File(filepath, "r") as f:
        kernel = np.asarray(f["dense_1/dense_1/kernel:0"][()])  # shape (1891, 5)
    return kernel


def plot_per_neuron_index_plots(kernel, out_path="dense1_per_neuron_histograms.png"):
    n_in, n_out = kernel.shape  # (1891, 5)
    fig, axes = plt.subplots(1, n_out, figsize=(4 * n_out, 4), sharey=True)

    for neuron_idx in range(n_out):
        weights = kernel[:, neuron_idx]      # the 1891 weights feeding this neuron
        idx = np.arange(weights.size)        # weight identifier: 0..1890
        ax = axes[neuron_idx]
        ax.bar(idx, weights, width=1.0, color="#4C72B0", linewidth=0)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_title(f"Output neuron {neuron_idx}\n(n={weights.size})")
        ax.set_xlabel("Weight index")
        if neuron_idx == 0:
            ax.set_ylabel("Weight value")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    filepath = sys.argv[1] if len(sys.argv) > 1 else "weights.h5"
    kernel = load_dense1_kernel(filepath)
    print(f"dense_1 kernel shape: {kernel.shape}")
    plot_per_neuron_index_plots(kernel)
