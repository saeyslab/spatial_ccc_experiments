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

## Run all methods on all datasets

## Compare outputs

Either visually or with a metric
