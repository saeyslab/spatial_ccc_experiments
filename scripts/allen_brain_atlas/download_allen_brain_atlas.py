import pooch
import sys

# Set global variables.
print("sys.argv", sys.argv)

out_path = sys.argv[1]

# loading raw expression matirx with pooch
print(f"Loading raw expression matrix to '{out_path}raw/'")
expression_matrix_raw = pooch.retrieve(
    url="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/expression_matrices/Zhuang-ABCA-1/20230830/Zhuang-ABCA-1-raw.h5ad",
    fname="Zhuang-ABCA-1-raw.h5ad",
    path=f"{out_path}raw/",
    known_hash=None,
)

# loading annotation with pooch
print(f"Loading annotation to '{out_path}raw/'")
cell_metadata = pooch.retrieve(
    url="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/Zhuang-ABCA-1/20231215/views/cell_metadata_with_cluster_annotation.csv",
    fname="cell_metadata_with_cluster_annotation.csv",
    path=f"{out_path}raw/",
    known_hash=None,
)

gene_metadata = pooch.retrieve(
    url="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/Zhuang-ABCA-1/20231215/gene.csv",
    fname="gene_metadata.csv",
    path=f"{out_path}raw/",
    known_hash=None,
)

print("Done. All files loaded.")
