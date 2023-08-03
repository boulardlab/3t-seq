#!/usr/bin/env bash

set -e
set -e  

pwd

sleep $(( $RANDOM % 10 + 2 ))

URL="${snakemake_params[url]}"
TMP=$(mktemp -p ${snakemake_params[tmp]})

OUTPUT=${snakemake_output}
mkdir -p $(dirname $OUTPUT)

if [ ${URL: -3} == ".gz" ] && [ ! ${OUTPUT: -3} == ".gz" ]; then
    wget --quiet -O - "$URL" | gunzip -c > $TMP 
else
    wget --quiet -O $TMP "$URL"
fi

if ! grep -v '#' $TMP | head -n 1 | grep -q 'chr'; then
    echo "Adding \"chr\" to first column, then mv'ing to $OUTPUT"
    awk '{ if ("#"~$0) {print $0} else {print "chr"$0}}' $TMP > $OUTPUT 
else
    echo "Mv'ing to $OUTPUT"
    mv $TMP $OUTPUT
fi
