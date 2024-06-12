## VIASH START
par <- list(
  input = "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
  output = "resources/results/spatalk/prediction.h5ad"
)
## VIASH END
library(tidyverse)
library(SpaTalk)
library(anndata)

adata <- read_h5ad(par$input)
count_matrix <- adata$layers[["counts"]] %>% t %>%
  `rownames<-`(adata$var$gene_symbol)
coords <- adata$obsm %>% data.frame %>% 
  mutate(cell = colnames(count_matrix), .before=1) %>% 
  rename(x = spatial.1, y = spatial.2)

pyad <- adata$.__enclos_env__$private$.anndata
pyad$obs <- reticulate::py_get_item(pyad$obs, c("celltype", "fov"))

# create SpaTalk data
obj <- createSpaTalk(st_data = count_matrix,
                     st_meta = coords,
                     species = "Mouse",
                     if_st_is_sc = T,
                     spot_max_cell = 1,
                     celltype = as.character(adata$obs[["celltype"]]))

obj <- find_lr_path(object = obj, lrpairs = lrpairs, pathways = pathways)
obj <- dec_cci_all(object = obj)

# generate prediction
ccc_pred <- data.frame(
  obj@lrpair 
) %>% arrange(desc(score))

# create output
output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    ccc_pred = ccc_pred
  )
)
output$write_h5ad(par$output, compression = "gzip")
