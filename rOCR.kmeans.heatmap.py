#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale

f = 'SRR.filtered.list'
info = pd.read_csv(f, sep='\t', header=None)
tissue = info.pop(5)
tissue_palette = sns.hls_palette(14, s=0.9, l=0.6)
tissue_lut = dict(zip(tissue.unique(), tissue_palette))
col_colors = tissue.map(tissue_lut)

# Define color maps for heatmaps
cmap_zscore = sns.color_palette("RdBu_r", 7)
cmap_log = sns.cubehelix_palette(8, start=3, rot=0, dark=0.1, light=1)

# Read cpkm and cluster data
cpkm = pd.read_csv('cpkm..matrix', sep='\t', header=None)
cpkm = cpkm.rename(columns={0: 'rocr'})
cluster_list = pd.read_csv('kmeans.cluster.list', sep='\t', header=None)
cluster_list.columns = ['rocr', 'cluster']

df = pd.merge(cpkm, cluster_list, on='rocr')
df_sorted = df.sort_values(by='cluster')
rocr = df_sorted.pop('rocr')
clusters = df_sorted.pop('cluster')

# Define color palette for clusters
num_clusters = len(clusters.unique())
cluster_palette = sns.color_palette("cubehelix", num_clusters)
cluster_lut = dict(zip(clusters.unique(), cluster_palette))
row_colors = clusters.map(cluster_lut)

# Order of biosamples for heatmap
sample_order = [5, 12, 15, 21, 4, 2, 23, 1, 20, 7, 17, 14, 19, 16, 25, 9, 0, 10, 24, 6, 3, 13, 18, 8, 11, 22]

# Log transformation and normalization
log_df = np.log10(df_sorted.iloc[:, sample_order] + 0.1)
norm_df = log_df.apply(scale, axis=0)


name1 = 'rOCR.zscore.heatmap.png'
name2 = 'rOCR.log.cpkm.heatmap.png'
# Generate z-score normalized heatmap
plt.figure(figsize=(10, 10))
sns.clustermap(norm_df, col_colors=col_colors.values[sample_order], xticklabels=tissue[sample_order],
               row_colors=row_colors.values, yticklabels=False, row_cluster=False, col_cluster=False,
               cmap=cmap_zscore, vmin=-2, vmax=2).savefig(name1, dpi=100)
plt.close()

# Generate log-transformed CPKM heatmap
plt.figure(figsize=(10, 10))
sns.clustermap(log_df, col_colors=col_colors.values[sample_order], xticklabels=tissue[sample_order],
               row_colors=row_colors.values, yticklabels=False, row_cluster=False, col_cluster=False,
               cmap=cmap_log, vmax=2).savefig(name2, dpi=100)
plt.close()
