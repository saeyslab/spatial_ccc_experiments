## VIASH START
par <- list(
  input = "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
  output = "resources/results/example_python/prediction.h5ad"
)
## VIASH END

adata <- anndata::read_h5ad(par$input)

# generate prediction
ccc_pred <- data.frame(
  ...
)

# create output
output <- ad.AnnData(
  shape = c(0L, 0L),
  uns = list(
    ccc_pred = ccc_pred
  )
)
output$write_h5ad(par$output, compression = "gzip")