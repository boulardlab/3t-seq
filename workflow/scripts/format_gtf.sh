#!/bin/bash
set -euo pipefail

# Inputs
INPUT_GTF=$1
OUTPUT_GTF=$2
SELECTED_CHRS=$3

# Decompress if needed
if [[ "$INPUT_GTF" == *.gz ]]; then
    gunzip -c "$INPUT_GTF" > "${OUTPUT_GTF}.tmp"
else
    cp "$INPUT_GTF" "${OUTPUT_GTF}.tmp"
fi

HAS_CHR=$(head -n 1000 "${OUTPUT_GTF}.tmp" | grep -v '^#' | head -n1 | grep -q '^chr' && echo "yes" || echo "no")

# Subset chromosomes if requested
if [ -n "$SELECTED_CHRS" ]; then
    if [ "$HAS_CHR" == "no" ]; then
        CHRS_TO_EXTRACT=$(echo "$SELECTED_CHRS" | sed 's/chr//g')
        awk -v chrs="$CHRS_TO_EXTRACT" 'BEGIN{split(chrs,a,","); for(i in a) keep[a[i]]} /^#/ || $1 in keep' "${OUTPUT_GTF}.tmp" > "${OUTPUT_GTF}.tmp2"
        sed '/^[^#]/s/^/chr/' "${OUTPUT_GTF}.tmp2" > "$OUTPUT_GTF"
        rm "${OUTPUT_GTF}.tmp2"
    else
        awk -v chrs="$SELECTED_CHRS" 'BEGIN{split(chrs,a,","); for(i in a) keep[a[i]]} /^#/ || $1 in keep' "${OUTPUT_GTF}.tmp" > "$OUTPUT_GTF"
    fi
    rm "${OUTPUT_GTF}.tmp"
else
    if [ "$HAS_CHR" == "no" ]; then
        sed '/^[^#]/s/^/chr/' "${OUTPUT_GTF}.tmp" > "$OUTPUT_GTF"
        rm "${OUTPUT_GTF}.tmp"
    else
        mv "${OUTPUT_GTF}.tmp" "$OUTPUT_GTF"
    fi
fi
