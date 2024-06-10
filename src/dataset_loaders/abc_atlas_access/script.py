from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

## VIASH START
par = {
    "dataset_label": "WMB-10Xv3",
    "file_name": "WMB-10Xv3-CB-log2.h5ad",
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

assert len(data_files) == 26

abc_cache.get_data_path(directory=par["dataset_label"], file_name=par["file_name"])