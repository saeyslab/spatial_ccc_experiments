# install the R anndata package
install.packages("anndata")

# to install miniconda
reticulate::install_miniconda()

# to install the python anndata package
anndata::install_anndata()

# download the example dataset
# experiment/download_test_dataset.sh

# to read the example dataset
library(anndata)

adata <- anndata::read_h5ad("resources/datasets/Puck_200104_15.h5ad")

# view structure
adata

# AnnData object with n_obs × n_vars = 25741 × 25391
#     obs: 'mapped_reference_annotation', 'donor_id', 'self_reported_ethnicity_ontology_term_id', 'donor_living_at_sample_collection', 'organism_ontology_term_id', 'sample_uuid', 'sample_preservation_method', 'tissue_ontology_term_id', 'development_stage_ontology_term_id', 'sample_derivation_process', 'sample_source', 'donor_BMI_at_collection', 'tissue_type', 'tissue_section_uuid', 'tissue_section_thickness', 'library_uuid', 'assay_ontology_term_id', 'is_primary_data', 'cell_type_ontology_term_id', 'author_cell_type', 'disease_ontology_term_id', 'reported_diseases', 'sex_ontology_term_id', 'suspension_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
#     var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length'
#     uns: 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'spatial', 'title'
#     obsm: 'X_spatial', 'spatial'

# view obs
head(adata$obs)

# view var
head(adata$var)

# view obsm['X_spatial']
head(adata$obsm$X_spatial)
