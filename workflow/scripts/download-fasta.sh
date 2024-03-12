#!/usr/bin/env bash

set -e  

pwd

sleep $(( $RANDOM % 10 + 2 ))

URL="${snakemake_params[url]}"
OUTPUT="${snakemake_output}"

if [[ "$URL" = *.gz && ! "$OUTPUT" = *.gz ]]; then
    T=$(mktemp)
    wget --quiet -O $T "$URL" 
    gunzip -c $T > $OUTPUT
    rm $T
else
    wget --quiet -O "$OUTPUT" "$URL"
fi