#!/usr/bin/env bash

set -e
set -e  

pwd

sleep $(( $RANDOM % 10 + 2 ))

URL="${snakemake_params[url]}"
OUTPUT="${snakemake_output}"

if [[ "$URL" = *.gz && ! "$OUTPUT" = *.gz ]]; then
    wget --quiet -O - "$URL" | gunzip -c > $OUTPUT
else
    wget --quiet -O "$OUTPUT" "$URL"
fi