#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate refgenie
export REFGENIE=/Users/francesco/scratch/refgenie/genome_config.yaml
echo "Testing Fasta remote:"
refgenie seekr mm10/fasta || echo "failed"
echo "Testing ensembl_gtf remote:"
refgenie seekr mm10/ensembl_gtf || echo "failed"
echo "Testing gencode_gtf remote:"
refgenie seekr mm10/gencode_gtf || echo "failed"
