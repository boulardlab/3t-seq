#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake-8.9.0
ROOT=$(realpath .)
# Temporarily override validation using a custom flag or simply testing `--rulegraph` or disabling the check in common.smk
# Actually, the simplest way is to just run `snakemake -n` on the unit test config since it expects dummy files, or use the `test_snakemake_dag.yaml` with some dummy files created.
# Let's just create the missing dummy files so the validation passes.
mkdir -p .tests/integration/GSE130735-subset
touch .tests/integration/GSE130735-subset/SRX5795113_SRR9016959_1.fastq.gz
touch .tests/integration/GSE130735-subset/SRX5795113_SRR9016959_2.fastq.gz
touch .tests/integration/GSE130735-subset/SRX5795117_SRR9016963_1.fastq.gz
touch .tests/integration/GSE130735-subset/SRX5795117_SRR9016963_2.fastq.gz
touch .tests/integration/GSE130735-subset/SRX5795118_SRR9016964_1.fastq.gz
touch .tests/integration/GSE130735-subset/SRX5795118_SRR9016964_2.fastq.gz

snakemake \
  --directory .tests/integration \
  --configfile .tests/integration/config.yaml \
  --profile .tests/integration/profile \
  --snakefile workflow/Snakefile \
  --all-temp \
  --singularity-args="--bind $ROOT --bind $HOME" \
  --dry-run
