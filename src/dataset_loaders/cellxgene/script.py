import anndata as ad
import scanpy as sc
import urllib.request

## VIASH START
par = {
    "url": "https://datasets.cellxgene.cziscience.com/2ae703d3-e98f-4db0-8477-43b75f7c543c.h5ad",
    "output": "resources/datasets/Puck_200104_15.h5ad",
}
meta = {
    "temp_dir": "/tmp"
}
## VIASH END

# download the file
path = meta["temp_dir"] + "/Puck_200104_15.h5ad"

urllib.request.urlretrieve(par["url"], path)

# load the data
adata = ad.read_h5ad(path)

# move counts
adata.layers["counts"] = adata.X.copy()
adata.layers["normalized"] = adata.X.copy()
del adata.X

# compute log1p
sc.pp.normalize_total(adata, target_sum=1e4, layer="normalized")
sc.pp.log1p(adata, layer="normalized")

# AnnData object with n_obs × n_vars = 25741 × 25391
#     obs: 'mapped_reference_annotation', 'donor_id', 'self_reported_ethnicity_ontology_term_id', 'donor_living_at_sample_collection', 'organism_ontology_term_id', 'sample_uuid', 'sample_preservation_method', 'tissue_ontology_term_id', 'development_stage_ontology_term_id', 'sample_derivation_process', 'sample_source', 'donor_BMI_at_collection', 'tissue_type', 'tissue_section_uuid', 'tissue_section_thickness', 'library_uuid', 'assay_ontology_term_id', 'is_primary_data', 'cell_type_ontology_term_id', 'author_cell_type', 'disease_ontology_term_id', 'reported_diseases', 'sex_ontology_term_id', 'suspension_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
#     var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length'
#     uns: 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'spatial', 'title'
#     obsm: 'X_spatial', 'spatial'
#     layers: 'counts', 'normalized'

adata.write_h5ad(par["output"], compression="gzip")

