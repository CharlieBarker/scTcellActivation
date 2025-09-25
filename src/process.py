import GEOparse
import scanpy as sc
import pandas as pd
import urllib.request
from pathlib import Path
from scipy import sparse
import gzip
import shutil
import time
from scipy import io
import anndata as ad

# ---------------------------
# Paths
# ---------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DATA_RAW.mkdir(parents=True, exist_ok=True)
DATA_PROCESSED.mkdir(parents=True, exist_ok=True)


# ---------------------------
# Robust loader for dense matrix.txt.gz
# ---------------------------
def load_dense_gsm(gse_id):
    print(f"Processing GSE {gse_id} ...")
    gse_folder = DATA_RAW / gse_id
    gse_folder.mkdir(parents=True, exist_ok=True)

    # Load GEO metadata
    gse = GEOparse.get_GEO(geo=gse_id, destdir=gse_folder)
    gse.download_supplementary_files(directory=gse_folder)

    saved_files = []
    records = []  # collect metadata here

    # Loop over GSMs
    for gsm_id, gsm in gse.gsms.items():
        print(f"Processing {gsm_id}: {gsm.metadata.get('title', '')}")

        # Find the supplementary folder that matches this GSM
        gsm_folder = next((f for f in gse_folder.iterdir() if gsm_id in f.name and f.is_dir()), None)
        if gsm_folder is None:
            print(f"No folder found for {gsm_id}, skipping.")
            continue

        # Initialize AnnData
        adata = None

        # ----------------
        # dense matrix
        # ----------------
        matrix_file = next(gsm_folder.glob("*matrix.txt.gz"), None)
        if matrix_file:
            # Read the table, force numeric values
            with gzip.open(matrix_file, 'rt') as f:
                df = pd.read_csv(f, sep='\t', index_col=0)

            # Convert all values to numeric, coercing errors (turn non-numeric to 0)
            df = df.apply(pd.to_numeric, errors='coerce').fillna(0)

            # Check orientation: genes as rows, cells as columns
            X_sparse = sparse.csr_matrix(df.T.values)  # transpose to cells x genes

            # Create AnnData
            adata = ad.AnnData(X=X_sparse)
            adata.var_names = df.index  # genes
            adata.obs_names = df.columns  # cells


        if adata is None:
            print(f"No Matrix files found for {gsm_id}, skipping.")
            continue

        # ----------------
        # Save h5ad
        # ----------------
        processed_dir = DATA_PROCESSED
        processed_dir.mkdir(parents=True, exist_ok=True)
        save_path = processed_dir / f"{gsm_id}.h5ad"
        adata.write(save_path)
        print(f"Saved {save_path}")
        saved_files.append(save_path)

        # ----------------
        # Collect metadata
        # ----------------
        meta = gsm.metadata.copy()
        meta["gsm_id"] = gsm_id
        meta["gse_id"] = gse_id
        meta = {k: v[0] if isinstance(v, list) and len(v) == 1 else v for k, v in meta.items()}  # flatten
        records.append(meta)

    # ----------------
    # Save metadata CSV
    # ----------------
    if records:
        meta_df = pd.DataFrame(records)
        meta_dir = PROJECT_ROOT / "data" / "meta"
        meta_dir.mkdir(parents=True, exist_ok=True)
        meta_path = meta_dir / f"{gse_id}_meta.csv"
        meta_df.to_csv(meta_path, index=False)
        print(f"Metadata saved to {meta_path}")

    if len(saved_files) == 0:
        print(f"No valid GSMs found for {gse_id}.")
        return None

    return saved_files


# ---------------------------
# Memory-efficient process_gse function
# ---------------------------

def process_gse(gse_id):
    print(f"Processing GSE {gse_id} ...")
    gse_folder = DATA_RAW / gse_id
    gse_folder.mkdir(parents=True, exist_ok=True)

    # Load GEO metadata
    gse = GEOparse.get_GEO(geo=gse_id, destdir=gse_folder)
    gse.download_supplementary_files(directory=gse_folder)

    saved_files = []
    records = []  # collect metadata here

    # Loop over GSMs
    for gsm_id, gsm in gse.gsms.items():
        print(f"Processing {gsm_id}: {gsm.metadata.get('title', '')}")

        # Find the supplementary folder that matches this GSM
        gsm_folder = next((f for f in gse_folder.iterdir() if gsm_id in f.name and f.is_dir()), None)
        if gsm_folder is None:
            print(f"No folder found for {gsm_id}, skipping.")
            continue

        # Initialize AnnData
        adata = None

        # ----------------
        # RNA
        # ----------------
        rna_file = next(gsm_folder.glob("*RNA*matrix.mtx.gz"), None)
        if rna_file:
            features_file = rna_file.parent / rna_file.name.replace("-matrix.mtx.gz","-features.tsv.gz")
            barcodes_file = rna_file.parent / rna_file.name.replace("-matrix.mtx.gz","-barcodes.tsv.gz")

            matrix = io.mmread(rna_file).T.tocsr()
            features = pd.read_csv(features_file, header=None, sep='\t')
            barcodes = pd.read_csv(barcodes_file, header=None, sep='\t')

            adata = ad.AnnData(X=matrix)
            adata.var['gene_symbols'] = features[1].values
            adata.var_names = features[0]
            adata.obs_names = barcodes[0].values

        # ----------------
        # ADT
        # ----------------
        adt_file = next(gsm_folder.glob("*ADT*matrix.mtx.gz"), None)
        if adt_file and adata is not None:
            features_file = adt_file.parent / adt_file.name.replace("-matrix.mtx.gz","-features.tsv.gz")
            barcodes_file = adt_file.parent / adt_file.name.replace("-matrix.mtx.gz","-barcodes.tsv.gz")
            adt_matrix = io.mmread(adt_file).T.tocsr()
            adt_features = pd.read_csv(features_file, header=None, sep='\t')
            adata.obsm["protein_counts"] = adt_matrix.toarray()
            adata.uns["protein_names"] = adt_features[1].values

        # ----------------
        # HTO
        # ----------------
        hto_file = next(gsm_folder.glob("*HTO*matrix.mtx.gz"), None)
        if hto_file and adata is not None:
            features_file = hto_file.parent / hto_file.name.replace("-matrix.mtx.gz","-features.tsv.gz")
            barcodes_file = hto_file.parent / hto_file.name.replace("-matrix.mtx.gz","-barcodes.tsv.gz")
            hto_matrix = io.mmread(hto_file).T.tocsr()
            hto_features = pd.read_csv(features_file, header=None, sep='\t')
            adata.obsm["hto_counts"] = hto_matrix.toarray()
            adata.uns["hto_names"] = hto_features[1].values

        if adata is None:
            print(f"No RNA/ADT/HTO files found for {gsm_id}, skipping.")
            continue

        # ----------------
        # Save h5ad
        # ----------------
        processed_dir = DATA_PROCESSED
        processed_dir.mkdir(parents=True, exist_ok=True)
        save_path = processed_dir / f"{gsm_id}.h5ad"
        adata.write(save_path)
        print(f"Saved {save_path}")
        saved_files.append(save_path)

        # ----------------
        # Collect metadata
        # ----------------
        meta = gsm.metadata.copy()
        meta["gsm_id"] = gsm_id
        meta["gse_id"] = gse_id
        meta = {k: v[0] if isinstance(v, list) and len(v) == 1 else v for k, v in meta.items()}  # flatten
        records.append(meta)

    # ----------------
    # Save metadata CSV
    # ----------------
    if records:
        meta_df = pd.DataFrame(records)
        meta_dir = PROJECT_ROOT / "data" / "meta"
        meta_dir.mkdir(parents=True, exist_ok=True)
        meta_path = meta_dir / f"{gse_id}_meta.csv"
        meta_df.to_csv(meta_path, index=False)
        print(f"Metadata saved to {meta_path}")

    if len(saved_files) == 0:
        print(f"No valid GSMs found for {gse_id}.")
        return None

    return saved_files




# ---------------------------
# 3. Download & process each GSE
# ---------------------------
#h5ads_gse2 = process_gse("GSE164378")
h5ads_gse1 = load_dense_gsm("GSE126030")
h5ads_gse1 = load_dense_gsm("GSE197067")


