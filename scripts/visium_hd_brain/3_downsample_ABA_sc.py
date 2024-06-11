import os
import sys
import pandas as pd
import anndata as ad

if len(sys.argv) != 2:
    print("Usage: python downsample_ABA_sc.py <task_id>")
    sys.exit(1)

task_id = int(sys.argv[1])

# Define the paths
input_folder = os.environ['VSC_DATA_VO_USER'] + '/ABA_data/expression_matrices/WMB-10Xv3/20230630/'
downsampled_csv = os.environ['VSC_DATA_VO_USER'] + '/ABA_data/cell_metadata_with_cluster_annotation_downsampled.csv'
output_folder = os.environ['VSC_DATA_VO_USER'] + '/ABA_data/subset/'

# Ensure the output directory exists
os.makedirs(output_folder, exist_ok=True)

# Get the list of .h5ad files
h5ad_files = sorted([f for f in os.listdir(input_folder) if f.endswith('-raw.h5ad')])

# Ensure TASK_ID is within the correct range
if task_id < 1 or task_id > len(h5ad_files):
    raise ValueError(f"TASK_ID must be between 1 and {len(h5ad_files)}")

# Read the appropriate .h5ad file based on TASK_ID
input_file = os.path.join(input_folder, h5ad_files[task_id-1])
adata = ad.read_h5ad(input_file)

# Read the downsampled dataframe
downsampled_df = pd.read_csv(downsampled_csv)

# Subset the AnnData object based on the 'cell_label' column
cell_labels_to_keep = downsampled_df['cell_label'].values
subset_adata = adata[adata.obs_names.isin(cell_labels_to_keep)]

# Get basename of the input file
input_basename = os.path.basename(input_file)

# Save the subsetted AnnData object
output_file = os.path.join(output_folder, input_basename)
subset_adata.write(output_file)

print(f"Subset data saved to {output_file}")