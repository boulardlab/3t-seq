#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate refgenie
export REFGENIE=/Users/francesco/scratch/refgenie/genome_config.yaml
echo "Testing Fasta:"
refgenie seek mm10/fasta || echo "failed"
echo "Testing GTF / annotation:"
refgenie seek mm10/ensembl_gtf || echo "failed"
refgenie seek mm10/gencode_gtf || echo "failed"
refgenie seek mm10/gtf || echo "failed"
