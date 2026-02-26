#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate refgenie
export REFGENIE=/Users/francesco/scratch/refgenie/genome_config.yaml
refgenie listr
