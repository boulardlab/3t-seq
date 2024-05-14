#!/usr/bin/env bash

set -e
set -x

T=$(mktemp)
SORTED=$(mktemp -p ${snakemake_resources[tmpdir]})

CHROMOSOMES=$(awk '{print $1}' "${snakemake_input[genome]}")

for C in $CHROMOSOMES; do
    awk -v j="$C" '$1~j' "${snakemake_input[annotation]}" >> $T
done 

bedtools sort -faidx "${snakemake_input[genome]}" -i $T > $SORTED 2> ${snakemake_log}

bedtools coverage \
    -g "${snakemake_input[genome]}" \
    -sorted \
    -a $SORTED \
    -b "${snakemake_input[bam]}" | \
    sort -k1,1 -k2,2n > "${snakemake_output}" 2>> ${snakemake_log}

rm $SORTED
