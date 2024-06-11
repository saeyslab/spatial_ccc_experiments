from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import anndata as ad

## VIASH START
par = {
    "dataset_label": "WMB-10Xv3",
    "file_name": "WMB-10Xv3-CB/log2",
    "output": "resources/datasets/WMB-10Xv3-CB-log2.h5ad",
}
meta = {
    "temp_dir": "/tmp"
}
## VIASH END

cache_path = meta["temp_dir"] + "/ABA_data"

abc_cache = AbcProjectCache.from_s3_cache(cache_path)

# list all data files in the cache
data_files = abc_cache.list_data_files(par["dataset_label"])

# ['WMB-10Xv3-CB/log2', 'WMB-10Xv3-CB/raw', 'WMB-10Xv3-CTXsp/log2', 'WMB-10Xv3-CTXsp/raw', 'WMB-10Xv3-HPF/log2', 'WMB-10Xv3-HPF/raw', 'WMB-10Xv3-HY/log2', 'WMB-10Xv3-HY/raw', 'WMB-10Xv3-Isocortex-1/log2', 'WMB-10Xv3-Isocortex-1/raw', 'WMB-10Xv3-Isocortex-2/log2', 'WMB-10Xv3-Isocortex-2/raw', 'WMB-10Xv3-MB/log2', 'WMB-10Xv3-MB/raw', 'WMB-10Xv3-MY/log2', 'WMB-10Xv3-MY/raw', 'WMB-10Xv3-OLF/log2', 'WMB-10Xv3-OLF/raw', 'WMB-10Xv3-P/log2', 'WMB-10Xv3-P/raw', 'WMB-10Xv3-PAL/log2', 'WMB-10Xv3-PAL/raw', 'WMB-10Xv3-STR/log2', 'WMB-10Xv3-STR/raw', 'WMB-10Xv3-TH/log2', 'WMB-10Xv3-TH/raw']

assert len(data_files) == 26

# download the data to the cache and get the path
path = abc_cache.get_data_path(directory=par["dataset_label"], file_name=par["file_name"])

# load the data
adata = ad.read_h5ad(path)

adata

# AnnData object with n_obs × n_vars = 182026 × 32285
#     obs: 'cell_barcode', 'library_label', 'anatomical_division_label'
#     var: 'gene_symbol'
#     uns: 'normalization', 'parent', 'parent_layer', 'parent_rows'