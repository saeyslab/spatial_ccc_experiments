#/bin/bash

viash run src/dataset_loaders/allen_brain_atlas_via_pooch/config.vsh.yaml \
  --output resources/datasets/allen_brain_atlas_via_pooch/dataset.h5ad \
  --brain_section_label "Zhuang-ABCA-1.080"