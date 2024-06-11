#!/bin/bash

set -e

echo ">> Downloading resources"

viash run src/sync_resources/config.vsh.yaml -- \
  --input "s3://openproblems-data/resources/spatial_ccc_experiments/" \
  --output "resources" \
  --delete