#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake-8.9.0
ROOT=$(realpath .)
snakemake \
  --directory .tests/integration \
  --configfile .tests/integration/config.yaml \
  --profile .tests/integration/profile \
  --snakefile workflow/Snakefile \
  --all-temp \
  --singularity-args="--bind $ROOT --bind $HOME" \
  --dry-run
