#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake-8.9.0
python validate_schema.py
