import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import boto3
import os

## VIASH START
par = {
    "output": "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
    "brain_section_label": "Zhuang-ABCA-1.080",
}
meta = {
    "temp_dir": "cache"
}
## VIASH END

cache_path = meta["temp_dir"] + "/Zhuang_2023_merfish_brain"

# create dir if it does not exist
import os
os.makedirs(cache_path, exist_ok=True)

s3 = boto3.client("s3")

# loading raw expression matrix via s3
expression_matrix_raw = f"{cache_path}/Zhuang-ABCA-1-raw.h5ad"
print(f"Loading raw expression matrix to '{expression_matrix_raw}'", flush=True)
s3.download_file("allen-brain-cell-atlas", "expression_matrices/Zhuang-ABCA-1/20230830/Zhuang-ABCA-1-raw.h5ad", expression_matrix_raw)
adata = ad.read_h5ad(expression_matrix_raw)

# loading annotation with pooch
cell_metadata = f"{cache_path}/cell_metadata_with_cluster_annotation.csv"
print(f"Loading cell annotation to '{cache_path}'", flush=True)
s3.download_file("allen-brain-cell-atlas", "metadata/Zhuang-ABCA-1/20231215/views/cell_metadata_with_cluster_annotation.csv", cell_metadata)
df_cell = pd.read_csv(cell_metadata, index_col=0)

# loading gene metadata with pooch
gene_metadata = f"{cache_path}/gene_metadata.csv"
print(f"Loading gene metadata to '{gene_metadata}'", flush=True)
s3.download_file("allen-brain-cell-atlas", "metadata/Zhuang-ABCA-1/20231215/gene.csv", gene_metadata)
df_gene = pd.read_csv(f"{cache_path}/gene_metadata.csv", index_col=0)


print("Processing raw data", flush=True)
adata.var = df_gene

print("Adding annotation to adata.obs", flush=True)
adata.obs = adata.obs.join(df_cell, rsuffix="_right")
adata = adata[adata.obs.brain_section_label_right.notna()]

if par["brain_section_label"] is not None:
    print(f"Filtering for brain section '{par['brain_section_label']}'", flush=True)
    adata = adata[adata.obs.brain_section_label_right == par["brain_section_label"]].copy()

print("Adding data collection schema to adata.obs", flush=True)
adata.obs["dataset_id"] = "Zhuang_2023"
adata.obs["condition_id"] = "normal"
adata.obs["tissue"] = "brain"
adata.obs["organism"] = "Mus musculus"
adata.obs["assay_ontology"] = "MERFISH"
adata.obs["assay"] = "MERFISH"
adata.obs["celltype"] = adata.obs["class"].astype("category")
adata.obs["fov"] = adata.obs["brain_section_label"].astype("category")
adata.obsm["spatial"] = np.array(adata.obs[["x", "y"]])
del adata.obs["x"]
del adata.obs["y"]

print("Convert object columns to category", flush=True)
for obs_col in adata.obs.columns:
    if adata.obs[obs_col].dtype == "object":
        adata.obs[obs_col] = adata.obs[obs_col].astype("category")

print("Normalize counts", flush=True)
adata.layers["counts"] = adata.X.copy()
adata.layers["normalized"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4, layer="normalized")
sc.pp.log1p(adata, layer="normalized")
del adata.X

print(f"Writing processed anndata object to '{par['output']}'", flush=True)
adata.write_h5ad(par["output"], compression="gzip")

print("Done. All files processed.")