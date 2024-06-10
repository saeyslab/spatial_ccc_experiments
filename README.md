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

## Implement a selection of spatial CCC methods

## Run all methods on all datasets

## Compare outputs

