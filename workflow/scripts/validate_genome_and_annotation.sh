#!/usr/bin/env bash

# Disable strict mode
set +euo pipefail

# Enable debug printing to stdout
set -x

# Define input and log variables
fasta_file="${snakemake_input[genome_fasta_file]}"
gtf_file="${snakemake_input[genome_annotation_file]}"
log_file="${snakemake_log}"

# Check if input files exist
if [ ! -f "$fasta_file" ] || [ ! -f "$gtf_file" ]; then
    echo "Error: Input files do not exist or are inaccessible." >&2
    exit 1
fi

# Check fasta
if head -n1000 "$fasta_file" | grep -q '^>chr'; then
    echo "Fasta has chr" >> "$log_file"
    FASTA_HAS_CHR=0
else
    echo "Fasta does not have chr" >> "$log_file"
    FASTA_HAS_CHR=1
fi
echo "FASTA_HAS_CHR=$FASTA_HAS_CHR" > "$log_file"

# Check gtf
if head -n1000 "$gtf_file" | grep -v '#' | head -n1 | grep -q '^chr'; then
    echo "GTF has chr" >> "$log_file"
    GTF_HAS_CHR=0
else
    echo "GTF does not have chr" >> "$log_file"
    GTF_HAS_CHR=1
fi
echo "GTF_HAS_CHR=$GTF_HAS_CHR" >> "$log_file"

# If fasta has chr but gtf has not, add chr to gtf
if [ "$FASTA_HAS_CHR" -eq 0 ] && [ "$GTF_HAS_CHR" -eq 1 ]; then
    echo "Fixing GTF" >> "$log_file"
    sed -i -r 's/^([^#]+)$/chr\1/g' "$gtf_file"
fi

# If gtf has chr but fasta has not, add chr to fasta
if [ "$FASTA_HAS_CHR" -eq 1 ] && [ "$GTF_HAS_CHR" -eq 0 ]; then
    echo "Fixing fasta" >> "$log_file"
    sed -i -r 's/^>(.+)$/>chr\1/g' "$fasta_file"
fi
