import pandas as pd
import numpy as np
import scanpy as sc
import os
import matplotlib
matplotlib.use("Agg")  # non-interactive backend, no GUI pop-ups
import matplotlib.pyplot as plt

import seaborn as sns
import starcat
from starcat import starCAT
import mygene
import scipy.sparse as sp

# -------------------------------
# Set global figure directory
# -------------------------------
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR

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
sc.pl.umap(adata, color="leiden", save="_leiden.png")
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, save="_markers.png")
print("Marker genes identified and plot saved")

# -------------------------------
# Step 4: Manual annotation (template)
# -------------------------------
cluster_annotations = {
    '0': 'T cells',
    '1': 'B cells',
    '2': 'NK cells',
    '3': 'Monocytes',
    # Add/modify based on markers
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations).fillna("Unknown")

# -------------------------------
# Step 5: Run starCAT on raw counts
# -------------------------------
print("Step 5: Running starCAT...")
tcat = starCAT(reference='TCAT.V1', cachedir='./cache')

raw_adata = adata.copy()
raw_adata.X = raw_adata.layers["counts"]

# Map Ensembl IDs → gene symbols
mg = mygene.MyGeneInfo()
ensembl_ids = [x.split('.')[0] for x in raw_adata.var_names]
res = mg.querymany(ensembl_ids, scopes="ensembl.gene", fields="symbol", species="human")
mapping = {r["query"]: r.get("symbol", r["query"]) for r in res}
raw_adata.var_names = [mapping[x] for x in ensembl_ids]
raw_adata.var_names = raw_adata.var_names.astype(str)
raw_adata = raw_adata[:, ~raw_adata.var_names.duplicated()].copy()

# Drop cells with zero counts
mask = raw_adata.X.sum(axis=1).A1 > 0 if sp.issparse(raw_adata.X) else raw_adata.X.sum(axis=1) > 0
raw_adata = raw_adata[mask].copy()

# Run starCAT
usage, scores = tcat.fit_transform(raw_adata)

# Join results back to the Harmony-integrated adata
adata.obs = adata.obs.join(usage)
adata.obs = adata.obs.join(scores)

if "Multinomial_Label" in adata.obs:
    adata.obs["Multinomial_Label"] = adata.obs["Multinomial_Label"].astype(str)

# -------------------------------
# Step 6: Plotting (all on Harmony UMAP)
# -------------------------------
sc.pl.umap(adata, color="Multinomial_Label", save="_starcat_multinomial.png")
sc.pl.umap(adata, color="batch", save="_batch.png")
sc.pl.umap(adata, color=["ASA_binary", "Proliferation_binary"], save="_starcat_binaries.png")
sc.pl.umap(adata, color=sorted(usage.columns), ncols=3, vmin=0, vmax='p99', save="_starcat_continuous.png")


