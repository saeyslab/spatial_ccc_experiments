library(CellChat)
library(patchwork)

## VIASH START
par <- list(
  input = "resources/datasets/Zhuang_2023_merfish_brain/dataset.h5ad",
  output = "resources/results/Zhuang_2023_merfish_brain/prediction.h5ad"
)
## VIASH END

adata <- anndata::read_h5ad(par$input)

# Create CellChat object
data.input <- t(adata$layers[['counts']])
rownames(data.input) <- adata$var[['gene_symbol']]

meta = data.frame(labels = adata$obs[['celltype']], samples = unique(adata$obs[['brain_section_label']]), row.names = rownames(adata$obs))
spatial.locs <- as.data.frame(adata$obsm[['spatial']])
colnames(spatial.locs) <- c("imagerow", "imagecol")
conversion.factor = 1
spot.size = 10 # use the typical human cell size
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB) # Show the structure of the database 
dplyr::glimpse(CellChatDB$interaction) 
# use a subset of CellChatDB for cell-cell communication analysis 
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling 
CellChatDB.use <- CellChatDB # simply use the default CellChatDB to use all CellChatDB for cell-cell communication analysis 
# set the used database in the object 
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
cellchat <- updateCellChat(cellchat)
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database 
future::plan("multisession", workers = 4) # do parallel 
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE) 
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.mouse) ---> function not found

# Inference of cell-cell communication network
cellchat <- updateCellChat(cellchat)
# adjust the parameter scale.distance (default = 0.01) when working on data from spatial transcriptomics technologies (other than 10X Visium)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, distance.use = TRUE, interaction.range = 250, scale.distance = 100, contact.dependent = TRUE, contact.range = 100)

cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# Extract the inferred cellular communication network as a data frame
ccc_pred <- subsetCommunication(cellchat)

# Visualisation

# Visualise the aggregated network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(2,1))
aggr_net_count <- netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
aggr_net_weight <- netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Heatmaps
heatmap_count <- netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
heatmap_weight <- netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

# Visualization of cell-cell communication network
pathways.show <- c("TGFb") 

# Circle plot
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
tgfb_circle <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Spatial plot
par(mfrow=c(1,1))
tgfb_spatial <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

# Compute and visualize the network centrality scores:
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
tgfb_centrality_heatmap <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
tgfb_centrality_spatial <- netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 4, vertex.label.cex = 3.5)

# Visualize gene expression distribution on tissue
spatialFeaturePlot(cellchat, features = c("Tgfb2","ACVR1C_TGFbR2"), point.size = 0.8, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, pairLR.use = "NTS_NTSR1", point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)

# create output
output <- anndata::AnnData(
  shape = c(0L, 0L),
  uns = list(
    ccc_pred = ccc_pred
  )
)
output$write_h5ad(par$output, compression = "gzip")
