#!/bin/bash

set -e

echo ">> Downloading resources"

# OPTION 1
viash run src/sync_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/spatial_ccc_experiments/" \
  --output "resources" \
  --delete

# OPTION 2
aws s3 sync \
  s3://openproblems-data/resources/spatial_ccc_experiments/ \
  resources/ \
  --delete --dryrun