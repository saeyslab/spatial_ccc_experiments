# read example dataset
import anndata as ad

adata = ad.read_h5ad("resources/datasets/Puck_200104_15.h5ad")

# view structure
adata

# AnnData object with n_obs × n_vars = 25741 × 25391
#     obs: 'mapped_reference_annotation', 'donor_id', 'self_reported_ethnicity_ontology_term_id', 'donor_living_at_sample_collection', 'organism_ontology_term_id', 'sample_uuid', 'sample_preservation_method', 'tissue_ontology_term_id', 'development_stage_ontology_term_id', 'sample_derivation_process', 'sample_source', 'donor_BMI_at_collection', 'tissue_type', 'tissue_section_uuid', 'tissue_section_thickness', 'library_uuid', 'assay_ontology_term_id', 'is_primary_data', 'cell_type_ontology_term_id', 'author_cell_type', 'disease_ontology_term_id', 'reported_diseases', 'sex_ontology_term_id', 'suspension_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
#     var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length'
#     uns: 'citation', 'default_embedding', 'schema_reference', 'schema_version', 'spatial', 'title'
#     obsm: 'X_spatial', 'spatial'

# view obs
adata.obs

# view var
adata.var

# view obsm['X_spatial']
adata.obsm['X_spatial']