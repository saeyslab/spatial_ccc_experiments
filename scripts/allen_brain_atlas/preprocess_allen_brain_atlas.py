import scanpy as sc
import numpy as np
import pandas as pd
import sys

# Set global variables.
print("sys.argv", sys.argv)

path = sys.argv[1]

out_path = sys.argv[1]

print("Processing raw data")
adata = sc.read_h5ad(f"{path}raw/Zhuang-ABCA-1-raw.h5ad")
df_cell = pd.read_csv(f"{path}raw/cell_metadata.csv", index_col=0)
df_gene = pd.read_csv(f"{path}raw/gene_metadata.csv", index_col=0)

adata.var = df_gene

print("...adding annotation to adata.obs")
adata.obs = adata.obs.join(df_cell, rsuffix="_right")
adata = adata[adata.obs.brain_section_label_right.notna()]

df_cell_sub = df_cell[df_cell['brain_section_label'] == 'Zhuang-ABCA-1.080']

adata = adata[adata.obs.index.isin(list(df_cell_sub.index))]

print("...adding data collection schema to adata.obs")
adata.obs["dataset_id"] = "Zhuang_2023"
adata.obs["dataset_id"] = adata.obs["dataset_id"].astype("category")
adata.obs["sample_id"] = adata.obs["feature_matrix_label"].astype("category")
adata.obs["donor_id"] = adata.obs["donor_label"].astype("category")

adata.obs["condition_id"] = "normal"
adata.obs["condition_id"] = adata.obs["condition_id"].astype("category")

adata.obs["tissue"] = "brain"
adata.obs["tissue"] = adata.obs["tissue"].astype("category")

adata.obs["organism"] = "Mus musculus"
adata.obs["organism"] = adata.obs["organism"].astype("category")

adata.obs["assay_ontology"] = "MERFISH"
adata.obs["assay_ontology"] = adata.obs["assay_ontology"].astype("category")

adata.obs["assay"] = "MERFISH"
adata.obs["assay"] = adata.obs["assay"].astype("category")

adata.obs["celltype"] = adata.obs["class"].astype("category")
adata.obs["fov"] = adata.obs["brain_section_label"].astype("category")

adata.obsm["spatial"] = np.array(adata.obs[["x", "y"]])

print(f"Writing processed anndata object to '{out_path}/processed/'")
adata.write(f"{path}processed/yao_2023_merfish_brain.h5ad")

print("Done. All files processed.")