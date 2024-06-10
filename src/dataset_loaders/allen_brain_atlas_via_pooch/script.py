import pooch
import anndata as ad
import pandas as pd
import numpy as np

## VIASH START
par = {
    "output": "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
    "brain_section_label": "Zhuang-ABCA-1.080",
}
meta = {
    "temp_dir": "/tmp"
}
## VIASH END

cache_path = meta["temp_dir"] + "/Zhuang_2023_merfish_brain"

# loading raw expression matirx with pooch
print(f"Loading raw expression matrix to '{cache_path}'")
expression_matrix_raw = pooch.retrieve(
    url="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/Zhuang-ABCA-1/20230830/Zhuang-ABCA-1-raw.h5ad",
    fname="Zhuang-ABCA-1-raw.h5ad",
    path=cache_path,
    known_hash=None,
)

# loading annotation with pooch
print(f"Loading annotation to '{cache_path}'")
cell_metadata = pooch.retrieve(
    url="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/Zhuang-ABCA-1/20231215/cell_metadata.csv",
    fname="cell_metadata_with_cluster_annotation.csv",
    path=cache_path,
    known_hash=None,
)

gene_metadata = pooch.retrieve(
    url="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/Zhuang-ABCA-1/20231215/gene.csv",
    fname="gene_metadata.csv",
    path=cache_path,
    known_hash=None,
)

print("Processing raw data")
adata = ad.read_h5ad(f"{cache_path}/Zhuang-ABCA-1-raw.h5ad")
df_cell = pd.read_csv(f"{cache_path}/cell_metadata_with_cluster_annotation.csv", index_col=0)
df_gene = pd.read_csv(f"{cache_path}/gene_metadata.csv", index_col=0)

adata.var = df_gene

print("...adding annotation to adata.obs")
adata.obs = adata.obs.join(df_cell, rsuffix="_right")
adata = adata[adata.obs.brain_section_label_right.notna()]

if par["brain_section_label"] is not None:
    df_cell_sub = df_cell[df_cell['brain_section_label'] == par["brain_section_label"]]
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

print(f"Writing processed anndata object to '{par['output']}'")
adata.write_h5ad(par["output"], compression="gzip")

print("Done. All files processed.")