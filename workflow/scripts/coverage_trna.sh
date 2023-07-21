#!/usr/bin/env bash

set -e
set -x 

SORTED=$(mktemp)


CHROMOSOMES=$(awk '{print $1}' "${snakemake_input[genome]}")

awk -v j="$CHROMOSOMES" '$1~j' "${snakemake_input[annotation]}" | \
    grep -v 'chrUn' | \
    bedtools sort -faidx "${snakemake_input[genome]}" -i - \
    > $SORTED 2> ${snakemake_log}

bedtools coverage \
    -g "${snakemake_input[genome]}" \
    -sorted \
    -a $SORTED \
    -b "${snakemake_input[bam]}" | \
    sort -k1,1 -k2,2n > "${snakemake_output}" 2>> ${snakemake_log}

rm $SORTED
