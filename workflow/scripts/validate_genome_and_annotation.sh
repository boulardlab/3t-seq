#!/usr/bin/env bash

set +e

head -n1000 "${snakemake_input[genome_fasta_file]}" | grep -q '^>chr'
FASTA_HAS_CHR=$?

head -n1000 "${snakemake_input[genome_annotation_file]}" | grep -v '#' | head -n1 | grep -q '^chr'
GTF_HAS_CHR=$?

if [ $FASTA_HAS_CHR -eq 0 ] & [ $GTF_HAS_CHR -eq 1 ]; then
    sed -i -r 's/^([^#]+)$/chr\1/g' "${snakemake_input[genome_annotation_file]}"
fi

if [ $FASTA_HAS_CHR -eq 1 ] & [ $GTF_HAS_CHR -eq 0 ]; then
    sed -i -r 's/^>(.+)$/>chr\1/g' "${snakemake_input[genome_fasta_file]}"
fi


