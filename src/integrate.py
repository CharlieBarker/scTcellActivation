import scanpy as sc
import scipy.sparse as sp
import anndata as ad
import os
import numpy as np
import harmonypy as hm
import matplotlib.pyplot as plt
import mygene
import pandas as pd


# -------------------------------
# Step 0: Setup
# -------------------------------
data_dir = "/home/charlie/Desktop/tcell_activation_project/data/processed"


file_list = [
   "GSM5008737.h5ad",
   "GSM5008740.h5ad",
   "GSM3589414.h5ad"
]

file_list = [
   "GSM3589406.h5ad",
   "GSM3589407.h5ad",
   "GSM3589408.h5ad",
   "GSM3589409.h5ad",
   "GSM3589410.h5ad",
   "GSM3589411.h5ad",
   "GSM3589412.h5ad",
   "GSM3589413.h5ad",
   "GSM3589414.h5ad"
]


print("Step 0: Starting workflow")

# -------------------------------
# Step 1: Load datasets
# -------------------------------
print("Step 1: Loading datasets...")
adatas = [sc.read_h5ad(os.path.join(data_dir, f)) for f in file_list]
print(f"Loaded {len(adatas)} datasets")

# -------------------------------
# Step 2: Preprocess individually
# -------------------------------
for i, adata in enumerate(adatas):
    print(f"Step 2.{i+1}: Preprocessing {file_list[i]}")
    if sp.issparse(adata.X):
        # Replace NaNs in the stored data with 0
        adata.X.data = np.nan_to_num(adata.X.data)

    # 🔥 Save raw counts before any normalization/logging
    adata.layers["counts"] = adata.X.copy()

    sc.pp.filter_genes(adata, min_cells=3)  # remove genes with 0 counts
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

mg = mygene.MyGeneInfo()

for i, adata in enumerate(adatas):
    print(f"Processing AnnData object {i+1}/{len(adatas)}")

    # Only convert if var_names look like Ensembl IDs (start with ENSG)
    if all(x.startswith("ENSG") for x in adata.var_names):
        ensembl_ids = [x.split(".")[0] for x in adata.var_names]  # strip version numbers

        # Query MyGene.info for gene symbols
        res = mg.querymany(ensembl_ids, scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=True)
        
        # Map Ensembl -> Symbol; fallback to original Ensembl if symbol not found
        mapping = {ens: row['symbol'] if not pd.isna(row['symbol']) else ens 
                   for ens, row in res.iterrows()}
        
        # Apply mapping
        adata.var_names = [mapping[x.split(".")[0]] for x in adata.var_names]
        
        # Drop duplicate gene names (keep first)
        adata = adata[:, ~adata.var_names.duplicated()].copy()
        print(f"Converted Ensembl IDs to gene symbols. Now {adata.n_vars} genes.")
    else:
        print("Var names do not appear to be Ensembl IDs, skipping.")

    # Optional: update the layer names as well (they follow var_names automatically)
    if "counts" in adata.layers:
        adata.layers["counts"] = adata.layers["counts"][:, ~adata.var_names.duplicated()]

    # Save back into the list
    adatas[i] = adata

# -------------------------------
# Step 3: Concatenate datasets
# -------------------------------
print("Step 3: Concatenating datasets...")

# Use file names without .h5ad as batch keys
keys = [f.replace(".h5ad", "") for f in file_list]

adata_combined = ad.concat(
    adatas,
    join='inner',
    label='batch',
    keys=keys
)

print("Layers in combined object:", adata_combined.layers.keys())


# Make cell names unique to avoid warnings
adata_combined.obs_names_make_unique()

print(f"Combined dataset: {adata_combined.n_obs} cells x {adata_combined.n_vars} genes")
print("Batch distribution:")
print(adata_combined.obs['batch'].value_counts())


# -------------------------------
# Step 4: HVG selection
# -------------------------------
print("Step 4: Selecting highly variable genes...")
sc.pp.highly_variable_genes(adata_combined, flavor='seurat', n_top_genes=5000)
adata_combined = adata_combined[:, adata_combined.var.highly_variable]
print(f"Selected {adata_combined.n_vars} highly variable genes")

# -------------------------------
# Step 5: Subsample cells (optional for memory)
# -------------------------------
n_cells = 30000
if adata_combined.n_obs > n_cells:
    print(f"Step 5: Subsampling {n_cells} cells for memory efficiency...")
    idx = np.random.choice(adata_combined.n_obs, n_cells, replace=False)
    adata_combined = adata_combined[idx, :]
    print(f"Subsampled dataset: {adata_combined.n_obs} cells x {adata_combined.n_vars} genes")

# Step 6: Dimensionality reduction
print("Step 6: PCA...")
sc.pp.scale(adata_combined, max_value=10)
sc.tl.pca(adata_combined, svd_solver='arpack', n_comps=50)
print(np.cumsum(adata_combined.uns['pca']['variance_ratio'])[:10])

print("Step 7: Running Harmony...")

adata_combined.obs['batch'] = adata_combined.obs['batch'].astype(str)
ho = hm.run_harmony(adata_combined.obsm['X_pca'], adata_combined.obs, 'batch')
adata_combined.obsm['X_pca_harmony'] = ho.Z_corr.T

sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony', n_neighbors=15)

# -------------------------------
# Step 8: UMAP
# -------------------------------
print("Step 8: Computing UMAP...")
sc.tl.umap(adata_combined)
print("UMAP finished")

# -------------------------------
# Step 9: Plot UMAP colored by batch
# -------------------------------
print("Step 9: Plotting UMAP...")
sc.pl.umap(adata_combined, color='batch', save="_combat_umap.png")
print("UMAP plot saved as *_combat_umap.png")

# -------------------------------
# Step 10: Save integrated dataset
# -------------------------------
save_path = os.path.join(data_dir, "adata_combined_combat.h5ad")
adata_combined.write(save_path)
print(f"Integrated dataset saved to {save_path}")
