from pathlib import Path
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
import os

# Get argument from the command line
import sys
if len(sys.argv) != 2:
    print("Usage: python download_ABA_sc.py <task_id>")
    sys.exit(1)

task_id = int(sys.argv[1])

# Combine VSC_DATA_VO_USER and "ABA_data" to get the download base
download_base = 'resources/ABA_data'
abc_cache = AbcProjectCache.from_s3_cache(download_base)

data_files = abc_cache.list_data_files('WMB-10Xv3')

# Check if there are 26 files in the cache
assert len(data_files) == 26

data_file = data_files[task_id-1]

# Download data
abc_cache.get_data_path(directory='WMB-10Xv3', file_name=data_file)

