#!/bin/bash
set -euo pipefail

# Inputs
INPUT_FASTA=$1
OUTPUT_FASTA=$2
SELECTED_CHRS=$3

# Decompress if needed
if [[ "$INPUT_FASTA" == *.gz ]]; then
    gunzip -c "$INPUT_FASTA" > "${OUTPUT_FASTA}.tmp"
else
    cp "$INPUT_FASTA" "${OUTPUT_FASTA}.tmp"
fi

HAS_CHR=$(head -n 1000 "${OUTPUT_FASTA}.tmp" | awk '/^>/ {if ($0 ~ /^>chr/) {print "yes"; exit} else {print "no"; exit}}' || echo "no")

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
