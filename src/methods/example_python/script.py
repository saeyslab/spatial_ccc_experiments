import anndata as ad
import pandas as pd

## VIASH START
par = {
  "input": "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
  "output": "resources/results/example_python/prediction.h5ad"
}
## VIASH END

adata = ad.read_h5ad(par["input"])

# generate prediction
ccc_pred = pd.DataFrame(
  # ...
)

# create output
output = ad.AnnData(
  shape=(0, 0),
  uns={
    "ccc_pred": ccc_pred
  }
)
output.write_h5ad(par["output"], compression = "gzip")