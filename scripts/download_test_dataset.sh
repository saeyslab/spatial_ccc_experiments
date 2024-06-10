#!/bin/bash

[ ! -d resources/datasets ] && mkdir -p resources/datasets

if [ ! -f resources/datasets/Puck_200104_15.h5ad ]; then
    wget https://datasets.cellxgene.cziscience.com/2ae703d3-e98f-4db0-8477-43b75f7c543c.h5ad \
      -O resources/datasets/Puck_200104_15.h5ad
fi