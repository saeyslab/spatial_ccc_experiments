## VIASH START
par <- list(
  input = "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
  output = "resources/results/stMLnet/prediction.h5ad"
)

## VIASH END
source("src/methods/stMLnet/Select_Lig_Rec_TGs.R")
library(hdf5r)
library(tidyverse)

# Load data
h5_file <- H5File$new(par$input)

# Get some basic information
cell_ids <- h5_file[["obs"]][[h5attr(h5_file[["obs"]], "_index")]][]
# Usually this code, but it's ENSEMBL ID
#gene_ids <- h5_file[["var"]][[h5attr(h5_file[["var"]], "_index")]][]
gene_ids <- h5_file[["var"]][["gene_symbol"]][]
ncells <- h5_file[["layers"]][["counts"]]$maxdims[2]
ngenes <- h5_file[["layers"]][["counts"]]$maxdims[1]

# Get cell types
celltypes <-  h5_file[["obs"]][["celltype"]][["categories"]][]
celltype_names <- h5_file[["obs"]][["cellname"]][["categories"]][]
celltype_mapping <- data.frame(celltypes = celltypes,
           index = 0:(length(celltypes)-1))
celltype_list <- h5_file[["obs"]][["celltype"]][["codes"]][] %>%
  data.frame(index = .) %>% 
  # Map this with celltype_mapping
  left_join(celltype_mapping, by = 'index') %>%
  `rownames<-`(cell_ids) %>% 
  select(-index) %>%
  # Hardcoded as Cluster in the tool
  rename(Cluster=celltypes) 

# Get coordinates
coords <- t(h5_file[["obsm"]][["spatial"]][, 1:ncells]) %>% 
  `rownames<-`(cell_ids) %>% 
  as.data.frame()

# Get count matrix
count_matrix <- h5_file[["layers"]][["counts"]][1:ngenes, 1:ncells] %>% 
  `rownames<-`(gene_ids) %>%
  `colnames<-`(cell_ids)

# They only have human database
# load(url("https://zenodo.org/records/10901775/files/ex_databases.rda?download=1"))

# So, instead use NicheNet database
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds")) %>% 
  group_by(from, to) %>% tally() %>%  rename(source=from, target=to, score=n)
weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

# Put it into format of stMLNet prior database
nichenet_db <- list(
  LigRec.DB = lr_network,
  RecTF.DB = weighted_networks$lr_sig %>% rename(source=from, target=to, score=weight),
  TFTG.DB = weighted_networks$gr %>% rename(source=from, target=to, score=weight)
)

# Get LR-TF-TG of interest
# THIS CODE USES TOO MUCH MEMORY!
ex_inputs <- Select_Lig_Rec_TGs(ExprMat = count_matrix, 
                                AnnoMat = celltype_list, 
                                LocaMat = coords,
                                Databases = nichenet_db,
                                # python.exe path for Giotto
                                python_path = "/home/chananchidas/.local/share/r-miniconda/envs/giotto_env/bin/python", 
                                min.pct = 0.05, 
                                expr.ct = 0.1, 
                                pct.ct = 0.05)

# If you manage to get the above to run, this part should work
library(stMLnet)
library(ggalluvial)

outputDir <- getwd()

DistMat <- as.matrix(dist(ex_inputs$locaMat))
colnames(DistMat) <- colnames(ex_inputs$exprMat)
rownames(DistMat) <- colnames(ex_inputs$exprMat)

resList <- runstMLnet(ExprMat = ex_inputs$exprMat,
                      AnnoMat = ex_inputs$annoMat,
                      DistMat,
                      LigClus = NULL,
                      RecClus = celltype_names,
                      Normalize = F, 
                      OutputDir = outputDir,
                      TGList=ex_inputs$tgs_of_inter, 
                      LigList=ex_inputs$ligs_of_inter, 
                      RecList=ex_inputs$recs_of_inter)
