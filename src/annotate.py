import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import starcat
from starcat import starCAT
import mygene
import scipy.sparse as sp


# -------------------------------
# Step 0: Setup
# -------------------------------
data_dir = "/home/charlie/Desktop/tcell_activation_project/data/processed"
out_dir = "/home/charlie/Desktop/tcell_activation_project/data/annotated"

input_file = "adata_combined_combat.h5ad"
output_file = "adata_annotated.h5ad"

print("Step 0: Starting annotation workflow")

# -------------------------------
# Step 1: Load integrated dataset
# -------------------------------
print("Step 1: Loading integrated dataset...")
adata = sc.read_h5ad(os.path.join(data_dir, input_file))
print(f"Loaded dataset: {adata.n_obs} cells x {adata.n_vars} genes")

# -------------------------------
# Step 2: Clustering
# -------------------------------
print("Step 2: Clustering cells...")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden")
print("Clusters computed and stored in adata.obs['leiden']")

# -------------------------------
# Step 3: Compute marker genes
# -------------------------------
print("Step 3: Finding marker genes per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save="_markers.png")
print("Marker genes identified and plot saved")

# -------------------------------
# Step 4: Manual annotation (template)
# -------------------------------
# Map clusters to cell types after inspecting marker genes.
# Update the dictionary below with your annotations.

cluster_annotations = {
    '0': 'T cells',
    '1': 'B cells',
    '2': 'NK cells',
    '3': 'Monocytes',
    # Add/modify based on markers
}

adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations).fillna("Unknown")

#starcat tutorial
# Run starCAT to compute the usages and scores for the provided data 
tcat = starCAT(reference='TCAT.V1', cachedir='./cache')

raw_adata = adata.copy()
raw_adata.X = raw_adata.layers["counts"]

mg = mygene.MyGeneInfo()

ensembl_ids = [x.split('.')[0] for x in raw_adata.var_names]
res = mg.querymany(ensembl_ids, scopes="ensembl.gene", fields="symbol", species="human")
mapping = {r["query"]: r.get("symbol", r["query"]) for r in res}
raw_adata.var_names = [mapping[x] for x in ensembl_ids]
# Ensure var_names are strings
raw_adata.var_names = raw_adata.var_names.astype(str)

# Drop duplicate gene names (keep first)
raw_adata = raw_adata[:, ~raw_adata.var_names.duplicated()].copy()

# Remove cells with zero overlap genes
mask = raw_adata.X.sum(axis=1).A1 > 0 if sp.issparse(raw_adata.X) else raw_adata.X.sum(axis=1) > 0
raw_adata = raw_adata[mask].copy()

usage, scores = tcat.fit_transform(raw_adata)


# -------------------------------
# Step 5: UMAP on raw counts + starCAT labels
# -------------------------------

# Start with raw counts object
umap_adata = sc.AnnData(
    X=raw_adata.layers["counts"].copy(),
    obs=raw_adata.obs.copy(),
    var=raw_adata.var.copy()
)

# Normalize and log1p
sc.pp.normalize_total(umap_adata, target_sum=1e4)
sc.pp.log1p(umap_adata)

# HVGs
sc.pp.highly_variable_genes(umap_adata, n_top_genes=3000)
hvgs = umap_adata.var[umap_adata.var['highly_variable']].index
umap_adata = umap_adata[:, hvgs].copy()

# PCA
sc.pp.scale(umap_adata, max_value=10)
sc.tl.pca(umap_adata)

# Neighbors & UMAP
sc.pp.neighbors(umap_adata, n_neighbors=10, n_pcs=25)
sc.tl.umap(umap_adata)

# -------------------------------
# Step 6: Add starCAT outputs
# -------------------------------

# Add usage and scores to obs
umap_adata.obs = umap_adata.obs.join(usage)
umap_adata.obs = umap_adata.obs.join(scores)

# Ensure multinomial label is string (so scanpy can plot categorical colors)
if "Multinomial_Label" in umap_adata.obs:
    umap_adata.obs["Multinomial_Label"] = umap_adata.obs["Multinomial_Label"].astype(str)

# -------------------------------
# Step 7: Plotting
# -------------------------------
sc.pl.umap(umap_adata, color="Multinomial_Label", save="_starcat_multinomial.png")
sc.pl.umap(umap_adata, color=["ASA_binary", "Proliferation_binary"], save="_starcat_binaries.png")
