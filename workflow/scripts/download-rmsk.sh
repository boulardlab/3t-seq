#!/usr/bin/env bash

set -e
set -x 

pwd

sleep $(( $RANDOM % 10 + 2 ))

URL="${snakemake_params[url]}"
OUTPUT="${snakemake_output}"
if [[ "$URL" = *.gz && ! "$OUTPUT" = *.gz ]]; then
  wget -q -O - "$URL" | gunzip -c > $OUTPUT |& tee "${snakemake_log}"
else
  wget -q -O "$OUTPUT" "$URL" |& tee "${snakemake_log}"
fi
