#!/usr/bin/env bash

set -e

if [[ -n "$TEMPDIR" ]]; then
    TMPDIR="$TEMPDIR"
elif [[ -n "$TMP" ]]; then
    TMPDIR="$TMP"
elif [[ -z "$TMPDIR" ]]; then
    TMPDIR="/tmp"
fi

T=$(mktemp -d -p $TMPDIR)


I=""
for F in ${snakemake_input}; do
    BN=$(basename $F)
    if [[ $BN == *.gz ]]; then
        O=$T/${BN%.gz}
        O=${O/txt/fq}
        gunzip -c $F > $O
        I="$I $O"
    else
        I="$I $F"
    fi
done

python /opt/SalmonTE/SalmonTE.py quant \
--reference=${snakemake_params[reference_genome]} \
--outpath=${snakemake_output[outfolder]} \
--num_threads=${snakemake[threads]} $I |& \
tee ${snakemake_log}