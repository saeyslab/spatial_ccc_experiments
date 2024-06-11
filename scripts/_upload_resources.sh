#!/bin/bash

aws s3 sync \
  s3://openproblems-data/resources/spatial_ccc_experiments/datasets/ \
  resources/datasets/ \
  --delete --dryrun

aws s3 sync \
  resources/datasets/ \
  s3://openproblems-data/resources/spatial_ccc_experiments/datasets/ \
  --delete --dryrun