# spatial_ccc_experiments

This repository contains experiments with performing spatial ccc. Our current aims are to:

* List a few candidate datasets
* Decide on a common input format
* Decide on a common output format
* Implement a selection of spatial CCC methods
* Run all methods on all datasets
* Compare outputs (visually, or with a metric)

## Candidate datasets

## Input format

Related issue: [#2](https://github.com/saeyslab/spatial_ccc_experiments/issues/2)

```
AnnData object
 obs: 'celltype'
 obsm: 'spatial'
 layers: 'counts'
```

How to read in Python:

```python
import anndata as ad

adata = ad.read_h5ad('path/to/file.h5ad')

# Access celltype
adata.obs['celltype']

# Access spatial coordinates
adata.obsm['spatial']

# Access counts
adata.layers['counts']
```

How to read in R:

```R
library(anndata)

adata <- read_h5ad('path/to/file.h5ad')

# Access celltype
adata$obs[["celltype"]]

# Access spatial coordinates
adata$obsm[["spatial"]]

# Access counts
adata$layers[["counts"]]
```

## Output format

Related issue: [#1](https://github.com/saeyslab/spatial_ccc_experiments/issues/1)

Spatial CCC methods typically output one of three formats.

### Per-spot prediction scores

Examples: SpatialDM, LIANA+, NICHES

### Per slide

That is, one summary of interactions for the whole slide / sample

Examples: COMMOT / MISTy / LIANA+

### Per source-target

That is, a summary for the whole slide but with respect to the cell type annotations

Examples: Giotto / CellChatv2

## Implement a selection of spatial CCC methods

The goal of the group was to run multiple spatial CCC methods, compare evaluations/visualizations and results. We selected the methods from Armingol et al., 2024. 

| Tool           | Data | Input | Output | Synthetic data | Experimental | Visualizations | Benchmarked against | Comments |
|----------------|------|-------|--------|----------------|--------------|----------------|----------------------|----------|
| COMMOT         | Mouse MERFISH, mouse placenta & brain STARMap, mouse hippocampus Slide-seqV2, Visium | LR pairs from CellChat DB (Heteromeric/homomeric binding) or manually defined genes, Log1p counts | Sparse matrix (in adata.obsp) with interaction scores for cell to cell | LR binding | Immunostaining of proteins; RNAScope imaging | Spatial signaling direction between spots/cells, heatmaps of DGE in CCC | CellChat, Giotto, CellPhoneDB v.3 | Got output; Easy to set up and run, but took a long time to run; Weird output because every interaction is in separate files |
| stMLNet        | Breast cancer ST, seqFISH+, MERFISH, Slide-seq v2, Visium glioma & COVID-19 | Counts, location, and cell type annotations; method only provides human DB | Ligand-receptors scores for each sender-receiver-target gene combination | Data following ligand diffusion model | N/A | Sankey plots, CCC network, chord diagram | NicheNet, CytoTalk, MISTy | Ran on the example dataset from their documentation; Tried running on chosen dataset, but the preprocessing is quite slow; Returns folders per cell type pair |
| SpatialDM      | Visium intestine, melanoma, seqFISH subventricular zone | CellChat DB | SVCA (Gaussian process model) simulated ST data | Chord diagram, spatial plots for interaction pattern clustering | CellChat, Giotto, SpaTalk | Ran the global part; The local part, but it's currently error | Apparently LIANA+ already has SpatialDM included so we're running this – ran through and much faster, will check the results | Initial implementation is pretty slow |
| SpaTalk        | Many untargeted and targeted datasets | Counts, location, cell type annotation | Prioritized ligand-receptors for each sender-receiver pair; enriched pathways | Simulated spots for deconvolution | N/A | Chord diagrams, snakey plots, heatmaps | Giotto, SpaOTsc, NicheNet, CytoTalk, CellCall, CellPhoneDB, and CellChat | No Seurat v5 compatibility |
| HoloNet        | Breast cancer & liver cancer Visium, Slide-seq v2 | Strategy of Tang et al., generated spatial patterning of two cell types | N/A | CCC network overlaid on spatial data | NicheNet, SpaTalk | N/A | N/A | N/A |
| CellNeighborEX | seqFISH mouse embryo, slide-seq from mouse embryo, mouse liver, mouse brain | Artificial heterotypic spots | Cell isolation and sorting with specific markers | Spatial visualization of cell neighbor-dependent gene expression | N/A | N/A | N/A | N/A |
| LIANA+ (cellphonedb) | Visium, Slide-tags, Spatial metabolite-transcriptome (MALDI-MSI x Visium) | - | Cell type specific score of communication events (lr_means) | Intracellular signalling network, ugly dotplots, heatmaps | Slide-tags spatially-informed vs not, classification of malignant spots, classification of conditions | Received dotplots and heatmaps for interactions cell-type specific, not spot specific | Need to investigate the resulting matrix, didn't yet find the specific ligand-receptor pairs |
| MEBOCOST       | scRNAseq | Gene expression of enzymes and knowledge repository of enzyme-metabolite-sensor partners (scFEA, COMPASS, HMDB, literature mining) | Number of communication events and score (sum of -log10(FDR)) of detected metabolite-sensor communications between a pair of cell types| - | - | - | - | Got metabolite-abundance estimation, is not comparable to other outputs | N/A |
| NicheCompass   | seqFISH mouse organogenesis; Slide-seq mouse hippocampus | - | Circular plot | CellCharter; STACI; GraphST | Not very easy to run, because we encountered a bug | Not very responsive on GitHub | N/A |
| Giotto         | seqFISH+, merFISH, osmFISH, STARmap, Visium 10X brain & kidney, Slide-seq, t-cyCIF, MIBI, CODEX | a gene-by-cell count matrix and the spatial coordinates for the centroid position of each cell | score based on adjusted p-value and log2 fold change was used to rank a ligand-receptor pair across all identified cells of interacting cell types | Simulated seqFISH+ and spatial pattern | N/A | All types of plots, GiottoVisuals, GiottoViewer | N/A | Not very easy to run; Not user friendly; Working running the method |
| CellChat v2           | Visium 10X mouse brain and human fetal intestine | Normalized and log-transformed gene-by-cell matrix, spatial coordinates,Metadata (cell type labels and sample names), CellChatDB | Ligand, receptor, source, target, interaction probabilities, pathways | - | - | Circular plot, heatmap, spatial plot | - | Gene expression data should use gene symbols |

## Run all methods on all datasets

## Compare outputs

Either visually or with a metric
