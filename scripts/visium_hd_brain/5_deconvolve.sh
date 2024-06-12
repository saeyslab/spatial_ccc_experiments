nextflow run main.nf -profile hpc \
-params-file $VSC_DATA_VO_USER/spatial_ccc_experiments/scripts/visium_hd_brain/cell2location_params.yaml \
--mode run_dataset \
--sc_input $VSC_DATA_VO_USER/ABA_data/WMB-10Xv3_subset.h5ad \
--sp_input $VSC_DATA_VO_USER/rds/Visium_HD_Mouse_Brain_008um.h5ad \
--annot class --methods cell2location --skip_metrics --gpu
