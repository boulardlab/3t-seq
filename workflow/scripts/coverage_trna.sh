#!/usr/bin/env bash

set -e
set -x

T=$(mktemp)
SORTED=$(mktemp -p ${snakemake_resources[tmpdir]})

# Keep tRNA lines whose chromosome is in the .genome file (exact match, #52)
awk 'NR == FNR { keep[$1]; next } $1 in keep' \
    "${snakemake_input[genome]}" "${snakemake_input[annotation]}" > "$T"

bedtools sort -faidx "${snakemake_input[genome]}" -i $T > $SORTED

bedtools coverage \
    -g "${snakemake_input[genome]}" \
    -sorted \
    -a $SORTED \
    -b "${snakemake_input[bam]}" | \
    sort -k1,1 -k2,2n > "${snakemake_output}"

rm $SORTED
