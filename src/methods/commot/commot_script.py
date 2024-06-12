import anndata as ad
import pandas as pd
import commot as ct
import scanpy as sc
import numpy as np

## VIASH START
par = {
  "input": "/lustre/groups/ml01/datasets/2024_spatial_ccc_experiments_francesca.drummer/dataset.h5ad",
  "output": "/lustre/groups/ml01/datasets/2024_spatial_ccc_experiments_francesca.drummer/experiments/prediction_commot.h5ad"
}
## VIASH END

adata = ad.read_h5ad(par["input"])

print("Set log1p data...")
adata.X = adata.layers['normalized']
adata.var_names = adata.var["gene_symbol"]

print("Get LR database...")
df_cellchat = ct.pp.ligand_receptor_database(species='mouse', signaling_type='Secreted Signaling', database='CellChat')
df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata, min_cell_pct=0.01)
print(f"Nr. of LR pairs {df_cellchat_filtered.shape}...")

print("Run COMMOT...")
ct.tl.spatial_communication(adata,
    database_name='user_database', df_ligrec=df_cellchat_filtered, dis_thr=200, heteromeric=True, pathway_sum=True)

print("Save prediction adata...")
adata.write(par["output"])
