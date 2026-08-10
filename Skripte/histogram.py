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
import matplotlib.colors as mcolors

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
    
def makeheatmaps(data,n, typ, title=False,saveto=False):#, gridspec_kw={'width_ratios': [1, 1, 1,1,1, 0.1]}
    fig, ax = plt.subplots(1,data.shape[1],figsize=(30,20), layout='constrained',squeeze=False)
    print(data.shape)
    ax = ax[0]
    weight_mtx=np.zeros((n,n))
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.9))
    if title:
        fig.text(0.5, 0.95, title, transform=fig.transFigure,
            ha="center", fontsize=14)
        fig.text(0.5, 0.91, typ+" weights for each distance pair", transform=fig.transFigure,
            ha="center", fontsize=12, color="dimgray")
    vmax = np.abs(data).max()
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = None
    for k in range(data.shape[1]):
        z=0
        print(k)
        for i in range(n-1):
            j=i+1
            while j<n:
                print(i,j,z)
                weight_mtx[i,j]=data[z,k]
                z+=1
                j+=1
        #weight_mtx[np.isclose(weight_mtx, 0, atol=1)] = 0
        new_arr = np.argwhere((weight_mtx > 1) | (weight_mtx<-1))
        new_set = set([])
        im = ax[k].imshow(weight_mtx, cmap="RdBu", norm=norm)
        # ax[k].text(0.5, -0.3, str(new_arr), 
        #    ha='center', va='top', transform=ax[k].transAxes, fontsize=10)
        ax[k].set_title("output "+str(k))
    if saveto:
        weight_mtx.save(saveto)
        
    cax = ax[-1].inset_axes([1.05, 0.05, 0.05, 0.9])
    fig.colorbar(im,cax=cax, label="weight value", location="right", orientation="vertical")
    #plt.tight_layout()
    plt.savefig(out_path.replace("histogram","grid").replace(".png","_"+typ+".png"), dpi=300)
    plt.close(fig)

weightdir="/data/joscha/output/qtools/SRw3UZCwUl830gEIOhHRkw_newick/trained_models_test/2026_07_14__13_03_44/weights/"
h5_file=weightdir+"m5_weights.h5"#sys.argv[1]
treefile= "/data/joscha/Downloads/SRw3UZCwUl830gEIOhHRkw_newick.tree"
out_path= h5_file.replace(".h5", "_histogram.png")
outfile_prefix = '/data/joscha/output/qtools/'+ str(Path(treefile).stem).replace(".tree","")+"/"


vector_file = outfile_prefix+"vectors.csv"
data = pd.read_csv(vector_file)
data['mtxvector'] = data['mtxvector'].apply(ast.literal_eval)
mtxvectors = np.array(data['mtxvector'].tolist())
avg_vector = mtxvectors.mean(axis=0)#.tolist()

with h5py.File(h5_file, "r") as f:
    weights = f[f"dense/dense/kernel:0"][:]
    bias = f[f"dense/dense//bias:0"][:]
number_of_bins=100
print(f"avg_input shape: {avg_vector.shape}  ",f"(must match weights.shape[0] = {weights.shape[0]})")

title="35 Taxa of TANGO1 with 62 gap-free columns"
entropy=shannon_entropy(weights)
plot_histogram(weights, title,
f"weights\n"+
f"Entropy={entropy:.3f} bits", out_path, bins = number_of_bins)


norm_weights= weights*avg_vector[:, None]
print(norm_weights.shape)
np.save(weightdir+"norm_weights.npy", norm_weights)
exit()
entropy_norm = shannon_entropy(norm_weights,bins=number_of_bins)
print(f"Entropy of normalized weights: {entropy_norm:.4f} bits")
plot_histogram(
    norm_weights,title,
    f"weights, normalized\n"+
    f"Entropy={entropy_norm:.3f} bits",
    out_path.replace(".png","_norm.png"),
    color="chocolate",bins = number_of_bins
)

n=int((sqrt(8*len(avg_vector)+1)+1)/2)
print("n=",n)

# makeheatmaps(weights,
#              n, 
#              "raw",
#              title=title)
# makeheatmaps(norm_weights,
#              n, 
#              "normalized",
#              title=title)
makeheatmaps(np.transpose(np.array([avg_vector])),
             n, 
             "averaged distances",
             title=title)
makeheatmaps(weights*np.transpose(np.array([avg_vector])),
              n, 
              "normalized experiment",
              title=title)
makeheatmaps(weights,
              n, 
              "raw weights",
              title=title)