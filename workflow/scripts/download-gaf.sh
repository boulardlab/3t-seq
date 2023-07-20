#!/usr/bin/env bash

set -e
set -x

pwd

sleep $(( $RANDOM % 10 + 2 ))

URL="${snakemake_params[url]}"
OUTPUT="${snakemake_output}"

if [[ "$URL" = *.gz ]]; then
    if [["$OUTPUT" = *.gz ]]; then
        wget -q -O - "$URL" | zgrep -v ! | pigz -c -p "${snakemake_threads}" > "$OUTPUT "
    else
        wget -q -O - "$URL" | zgrep -v ! >" $OUTPUT"
    fi
else
    if [["$OUTPUT" = *.gz ]]; then
        wget -q -O - "$URL" | grep -v ! | pigz -c -p "${snakemake_threads}" > "$OUTPUT" 
    else
        wget -q -O - "$URL" | grep -v ! > "$OUTPUT" 
    fi
fi 