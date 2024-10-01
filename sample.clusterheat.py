#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale

# Load tissue and species data
info = pd.read_csv('SRR.list', sep='\t', header=None)
tissue = info.pop(5)
species = info.pop(4)

# Define color palettes
tissue_colors = sns.hls_palette(15, s=0.9, l=0.6)
tissue_lut = dict(zip(tissue.unique(), tissue_colors))
row_colors = tissue.map(tissue_lut)

species_colors = sns.cubehelix_palette(4, start=.5, rot=-.75)
species_lut = dict(zip(species.unique(), species_colors))
row_colors2 = species.map(species_lut)

# Create the cluster heatmap using rOCR signal correlation
data = pd.read_csv('correlation.matrix', sep='\t', header=None)
g = sns.clustermap(data, row_colors=row_colors, col_colors=row_colors2,
                   xticklabels=species, yticklabels=tissue, vmin=0.4, vmax=1)


plt.savefig("all_rOCRs_clusterheat.svg", format="svg")
plt.close()
