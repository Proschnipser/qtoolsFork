#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 08:46:15 2026

@author: -
"""
import pandas as pd
from pathlib import Path
import ast
import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
from math import sqrt

def shannon_entropy(weights, bins='auto', value_range=None):
    counts, _ = np.histogram(weights.flatten(), bins='auto', range=value_range)
    print("Count shape",counts.shape)
    probs = counts / counts.sum()
    probs = probs[probs > 0]  # drop empty bins (log2(0) undefined)
    return -np.sum(probs * np.log2(probs))

def plot_histogram(weights, title,subtitle, out_path, color="lightblue", bins='auto'):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(weights.flatten(), bins=bins, color=color, edgecolor="black",linewidth=0.6)
    ax.text(0.5, 1.13, title, transform=ax.transAxes,
        ha="center", fontsize=11)
    ax.text(0.5, 1.03, subtitle, transform=ax.transAxes,
        ha="center", fontsize=9, color="dimgray")

    ax.set_xlabel("weight value")
    ax.set_ylabel("count")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    
def makeheatmaps(data,n, typ):
    fig, ax = plt.subplots(1,data.shape[1],figsize=(15, 4.5))
    weight_mtx=np.zeros((n,n))
    z=0
    for k in range(data.shape[1]):
        z=0
        for i in range(62):
            j=i+1
            while j<62:
                print(i,j,z)
                weight_mtx[i,j]=data[z,k]
                z+=1
                j+=1
        print(weight_mtx)
        ax[k].imshow(weight_mtx); ax[k].set_title(str(k))
    plt.tight_layout()
    plt.savefig(out_path.replace("histogram","grid").replace(".png","_"+typ+".png"), dpi=150)
    plt.close(fig)

h5_file="/data/joscha/output/qtools/SRw3UZCwUl830gEIOhHRkw_newick/trained_models_test/2026_06_24__08_29_04/weights/m4_weights.h5"#sys.argv[1]
treefile= "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
out_path= h5_file.replace(".h5", "_histogram.png")
outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(treefile).stem).replace(".tree","")+"/"

vector_file = outfile_prefix+"vectors.csv"
data = pd.read_csv(vector_file)
data['mtxvector'] = data['mtxvector'].apply(ast.literal_eval)
mtxvectors = np.array(data['mtxvector'].tolist())
avg_vector = mtxvectors.mean(axis=0)#.tolist()

with h5py.File(h5_file, "r") as f:
    weights = f[f"dense_1/dense_1/kernel:0"][:]
    bias = f[f"dense_1/dense_1//bias:0"][:]
number_of_bins=80
print(f"avg_input shape: {avg_vector.shape}  ",f"(must match weights.shape[0] = {weights.shape[0]})")

title="35 Taxa of TANGO1 with 62 gap-free columns"
entropy=shannon_entropy(weights)
plot_histogram(weights, title,
f"weights\n"+
f"Entropy={entropy:.3f} bits", out_path)

title=title+""

norm_weights= weights/avg_vector[:, None]

entropy_norm = shannon_entropy(norm_weights)
print(f"Entropy of normalized weights: {entropy_norm:.4f} bits")
plot_histogram(
    norm_weights,title,
    f"weights, normalized\n"+
    f"Entropy={entropy_norm:.3f} bits",
    out_path.replace(".png","_norm.png"),
    color="#55A868",
)
n=int((sqrt(8*len(avg_vector)+1)+1)/2)
print("n=",n)
makeheatmaps(weights,n, "weights")
makeheatmaps(weights,n, "norm")
