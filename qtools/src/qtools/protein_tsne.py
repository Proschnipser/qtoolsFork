"""
protein_tsne.py
───────────────
t-SNE visualisation of protein CIF structures based on pairwise Cα distance
distributions.

Usage
─────
    python protein_tsne.py --input /path/to/cif/folder

Or drop all your .cif files in the same folder as this script and run:
    python protein_tsne.py

Dependencies (install once)
───────────────────────────
    pip install biopython numpy scikit-learn matplotlib

How it works
────────────
1. Parse each .cif file → extract Cα (alpha-carbon) coordinates.
2. Compute the full N×N pairwise distance matrix for each structure.
3. Flatten the upper triangle into a histogram (100 distance bins, 0–150 Å).
   This gives every structure a fixed-length feature vector regardless of size.
4. Run t-SNE on the feature matrix.
5. Plot an interactive scatter where hovering shows the structure name, and
   points are coloured by protein length (number of Cα atoms).
"""

import argparse
import sys
import warnings
from pathlib import Path
import re
from natsort import  os_sorted, natsorted
import  data_prepper as dp

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from random import sample

# -- BioPython CIF parser ------------------------------------------------------
try:
    from Bio.PDB import MMCIFParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
except ImportError:
    sys.exit(
        "BioPython not found.\n"
        "Install it with:  pip install biopython"
    )

# -- scikit-learn t-SNE --------------------------------------------------------
try:
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler
except ImportError:
    sys.exit(
        "scikit-learn not found.\n"
        "Install it with:  pip install scikit-learn"
    )


# ------------------------------------------------------------------------------
# 1.  CIF -> Ca coordinates
# ------------------------------------------------------------------------------

def extract_ca_coords(cif_path):
    """Return (N, 3) array of Ca coordinates from a CIF file.
    Returns None if fewer than 5 Ca atoms are found.
    """
    parser = MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)
        try:
            structure = parser.get_structure(cif_path.stem, str(cif_path))
        except Exception as e:
            print(f"  [!] Could not parse {cif_path.name}: {e}")
            return None

    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    coords.append(residue["CA"].get_vector().get_array())
        break  # use only the first model

    if len(coords) < 5:
        print(f"  [!] Too few Ca atoms in {cif_path.name} ({len(coords)}), skipping.")
        return None

    return np.array(coords, dtype=np.float32)


# ------------------------------------------------------------------------------
# 2.  Distance matrix -> histogram feature vector
# ------------------------------------------------------------------------------

N_BINS   = 100
DIST_MIN = 0.0
DIST_MAX = 150.0   # Angstroms -- covers most intra-protein distances

def distance_histogram(coords):
    """Compute pairwise Ca distances, return a normalised histogram (N_BINS,)."""
    diff  = coords[:, None, :] - coords[None, :, :]   # (N,N,3)
    dists = np.sqrt((diff ** 2).sum(-1))               # (N,N)
    idx   = np.triu_indices(len(coords), k=1)
    upper = dists[idx]
    hist, _ = np.histogram(upper, bins=N_BINS, range=(DIST_MIN, DIST_MAX))
    total = hist.sum()
    return hist.astype(np.float32) / (total if total > 0 else 1.0)


# ------------------------------------------------------------------------------
# 3.  t-SNE
# ------------------------------------------------------------------------------

def run_tsne(feature_matrix, perplexity=None):
    """Scale features and run t-SNE. Perplexity is auto-tuned if not given."""
    n    = len(feature_matrix)
    perp = perplexity if perplexity else min(30, max(2, n // 3))
    print(f"\n  Running t-SNE  (n={n}, perplexity={perp}) ...")

    X = StandardScaler().fit_transform(feature_matrix)
    tsne = TSNE(
        n_components=2,
        perplexity=perp,
        n_iter=1000,
        random_state=42,
        init="pca",
        learning_rate="auto",
    )
    return tsne.fit_transform(X)


# ------------------------------------------------------------------------------
# 4.  Interactive plot
# ------------------------------------------------------------------------------

def plot_tsne(embedding, labels, lengths, output_path=None):
    """Scatter plot with hover labels and length-based colouring."""

    fig, ax = plt.subplots(figsize=(11, 8))
    fig.patch.set_facecolor("#0f0f1a")
    ax.set_facecolor("#0f0f1a")

    norm   = Normalize(vmin=min(lengths), vmax=max(lengths))
    cmap   = cm.plasma
    colors = cmap(norm(lengths))

    sc = ax.scatter(
        embedding[:, 0], embedding[:, 1],
        c=colors,
        s=90,
        alpha=0.85,
        linewidths=0.6,
        edgecolors="white",
        zorder=3,
    )

    for i, label in enumerate(labels):
        ax.annotate(
            label,
            (embedding[i, 0], embedding[i, 1]),
            fontsize=7,
            color="#e0e0e0",
            alpha=0.75,
            xytext=(5, 5),
            textcoords="offset points",
        )

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.03)
    cbar.set_label("Protein length (# Ca residues)", color="white", fontsize=10)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    ax.set_title(
        "Protein structure similarity — t-SNE  (Ca pairwise distance distributions)",
        color="white", fontsize=13, pad=14,
    )
    ax.set_xlabel("t-SNE dimension 1", color="#aaaaaa", fontsize=10)
    ax.set_ylabel("t-SNE dimension 2", color="#aaaaaa", fontsize=10)
    ax.tick_params(colors="#666666")
    for spine in ax.spines.values():
        spine.set_edgecolor("#333355")

    # Hover tooltip
    annot = ax.annotate(
        "", xy=(0, 0),
        xytext=(15, 15), textcoords="offset points",
        bbox=dict(boxstyle="round,pad=0.4", fc="#1e1e3a", ec="#9988ff", lw=1.2),
        fontsize=9, color="white",
        arrowprops=dict(arrowstyle="->", color="#9988ff"),
        zorder=10,
    )
    annot.set_visible(False)

    def on_motion(event):
        if event.inaxes != ax:
            annot.set_visible(False)
            fig.canvas.draw_idle()
            return
        cont, ind = sc.contains(event)
        if cont:
            i = ind["ind"][0]
            annot.xy = (embedding[i, 0], embedding[i, 1])
            annot.set_text(f"{labels[i]}\n{lengths[i]} residues")
            annot.set_visible(True)
        else:
            annot.set_visible(False)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", on_motion)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight",
                    facecolor=fig.get_facecolor())
        print(f"\n  Plot saved -> {output_path}")

    plt.show()


# ------------------------------------------------------------------------------
# 5.  Main
# ------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="t-SNE of protein CIF structures via Ca distance histograms."
    )
    ap.add_argument(
        "--input", "-i",
        type=Path, default=Path("."),
        help="Folder containing .cif files (default: current directory).",
    )
    ap.add_argument(
        "--perplexity", "-p",
        type=float, default=None,
        help="t-SNE perplexity (default: auto = min(30, n//3)).",
    )
    ap.add_argument(
        "--output", "-o",
        type=Path, default=None,
        help="Optional path to save the plot (e.g. tsne.png).",
    )
    args = ap.parse_args()

    #cif_files = sorted(args.input.glob("*.cif"))
    treefile = "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
    ref_tree = dp.read_tree(treefile)
    train_names = sample([leaf.name for leaf in ref_tree.get_terminals()],10)
    IDs= [re.search(r'_(\d+)_', name).group(1) for name in train_names]

    sorted_pairs = natsorted(zip(train_names, IDs), key=lambda pair: pair[1])
    train_names, IDs = zip(*sorted_pairs)
    print(train_names, IDs)

    pattern= '|'.join(IDs)
    regex = re.compile(rf'_({pattern})_\d+_')

    matches = os_sorted([p for p in Path(args.input).iterdir() if regex.search(p.name)])
    cif_files = [subdir.glob("*.cif") for subdir in matches]
    if not cif_files:
        sys.exit(f"No .cif files found in '{args.input}'.")

    print(f"\nFound {len(cif_files)} .cif file(s) in '{args.input}'\n")

    labels, features, lengths = [], [], []

    for cif in cif_files:
        print(f"  Parsing {cif.name} ...", end=" ", flush=True)
        coords = extract_ca_coords(cif)
        if coords is None:
            continue
        hist = distance_histogram(coords)
        print(f"{len(coords)} Ca atoms  OK")
        labels.append(cif.stem)
        features.append(hist)
        lengths.append(len(coords))

    if len(features) < 4:
        sys.exit(
            f"\nNeed at least 4 valid structures for t-SNE; "
            f"only {len(features)} parsed successfully."
        )

    embedding = run_tsne(np.vstack(features), args.perplexity)
    plot_tsne(embedding, labels, lengths, output_path=args.output)


if __name__ == "__main__":
    main()
