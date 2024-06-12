import os
import pandas as pd
import anndata as ad

import spatialdm as sdm
import liana as li

## VIASH START
par = {
  "input": "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
  "bandwidth": 130,
  "cutoff": 0.01,
  "species": "mouse",
  "output": "resources/results/example_python/prediction.h5ad"
}
## VIASH END

adata = ad.read_h5ad(par["input"])

# set X data to log1p
adata.X = adata.layers["normalized"]

# calculate spatial neighbors
li.ut.spatial_neighbors(adata, bandwidth=par["bandwidth"], cutoff=par["cutoff"], kernel='gaussian', set_diag=True)

# TODO: add resource_name to param and only run this if it is not in LIANA
# creating cellchatDB (since cellchatDB mouse is not in LIANA) 
min_cell = 0.01*adata.shape[0] # setting min cells to 1%
sdm.extract_lr(adata, species=par["species"], min_cell=min_cell)

ligand_receptor = adata.uns["receptor"].index.tolist()
ligands = [s.split("_")[0] for s in ligand_receptor]
receptors = [s.split("_")[1] for s in ligand_receptor]

resource = pd.DataFrame()
resource["ligand"] = ligands
resource["receptor"] = receptors

# generate prediction
# run global and local interactions
# TODO: only add resource if it is not in LIANA - otherwise, set resource_name
li.mt.bivariate(
    adata,
    #resource_name='', # NOTE: uses HUMAN gene symbols!
    resource=resource,
    local_name='morans', # Name of the function
    global_name="morans", # Name global function
    n_perms=100, # Number of permutations to calculate a p-value
    mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
    add_categories=True, # Whether to add local categories to the results
    nz_prop=0.2, # Minimum expr. proportion for ligands/receptors and their subunits
    use_raw=False,
    verbose=True
)


# create output
ccc_pred = adata.obsm['local_scores'] # NOTE: ccc_pred here is an AnnData object


output = ad.AnnData(
  shape=(0, 0),
  uns={
    "ccc_pred": ccc_pred
  }
)

output.write_h5ad(par["output"], compression = "gzip")