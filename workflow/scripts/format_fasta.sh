#!/usr/bin/env bash
set -euo pipefail

# Redirect all output and stderr to the snakemake log file
exec &> "${snakemake_log[0]}"

# Inputs
INPUT_FASTA="${snakemake_input[fasta]}"
OUTPUT_FASTA="${snakemake_output[fasta]}"
SELECTED_CHRS="${snakemake_params[selected_chrs]:-}"

# Decompress if needed
if [[ "$INPUT_FASTA" == *.gz ]]; then
    gunzip -c "$INPUT_FASTA" > "${OUTPUT_FASTA}.tmp"
else
    cp "$INPUT_FASTA" "${OUTPUT_FASTA}.tmp"
fi

# Check if chromosomes have 'chr' prefix
set +o pipefail
HAS_CHR=$(head -n 1000 "${OUTPUT_FASTA}.tmp" | awk '/^>/ {if ($0 ~ /^>chr/) {print "yes"; exit} else {print "no"; exit}}' || echo "no")
set -o pipefail

# Subset chromosomes if requested
if [ -n "$SELECTED_CHRS" ]; then
    if [ "$HAS_CHR" == "no" ]; then
        CHRS_TO_EXTRACT=$(echo "$SELECTED_CHRS" | sed 's/chr//g')
        samtools faidx "${OUTPUT_FASTA}.tmp" $CHRS_TO_EXTRACT > "${OUTPUT_FASTA}.tmp2"
        sed '/^>/s/^>/>chr/' "${OUTPUT_FASTA}.tmp2" > "$OUTPUT_FASTA"
        rm "${OUTPUT_FASTA}.tmp2"
    else
        samtools faidx "${OUTPUT_FASTA}.tmp" $SELECTED_CHRS > "$OUTPUT_FASTA"
    fi
    rm "${OUTPUT_FASTA}.tmp"
else
    if [ "$HAS_CHR" == "no" ]; then
        sed '/^>/s/^>/>chr/' "${OUTPUT_FASTA}.tmp" > "$OUTPUT_FASTA"
        rm "${OUTPUT_FASTA}.tmp"
    else
        mv "${OUTPUT_FASTA}.tmp" "$OUTPUT_FASTA"
    fi
fi
