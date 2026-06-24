#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 16:12:29 2026

@author: -
"""

import plotly.express as px
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.datasets import make_blobs
import plotly.io as pio
import webbrowser
import subprocess
webbrowser.register('firefox', None, webbrowser.GenericBrowser('/usr/bin/firefox'))
webbrowser.get('firefox')

# --- 1. Generate example data (120 points, 4 classes, 10 dimensions) ---
X, y = make_blobs(n_samples=120, centers=4, n_features=10,
                  cluster_std=0.8, random_state=42)

# --- 2. Fit t-SNE ---
tsne = TSNE(
    n_components=2,
    perplexity=20,        # 5–50; lower = tighter local clusters
    learning_rate=100,    # 10–300; 'auto' also works in sklearn ≥1.2
    max_iter=1000,
    random_state=42,
    init='pca'            # 'pca' is more stable than 'random'
)
X_2d = tsne.fit_transform(X)

colors  = ['#534AB7', '#1D9E75', '#D85A30', '#BA7517']
classes = ['Class A', 'Class B', 'Class C', 'Class D']


df = pd.DataFrame(X_2d, columns=['t-SNE 1', 't-SNE 2'])
df['Class'] = [classes[i] for i in y]


fig = px.scatter(df, x='t-SNE 1', y='t-SNE 2', color='Class',
                 color_discrete_sequence=colors,
                 title='Interactive t-SNE',
                 hover_data={'t-SNE 1': ':.2f', 't-SNE 2': ':.2f'})
fig.update_traces(marker=dict(size=9, opacity=0.85))
filepath="/data/joscha/output/qtools/plot.html"
fig.write_html(filepath)  # opens in browser; zoom, pan, hover all work out of the box
subprocess.Popen([
    "firefox",
    "--new-tab",
    filepath
])