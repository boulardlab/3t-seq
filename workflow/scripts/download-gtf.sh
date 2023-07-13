#!/usr/bin/env bash

set -e
set -x

URL="${snakemake_params[url]}"
TMP=$(mktemp -p ${snakemake_params[tmp]})

OUTPUT=${snakemake_output}
if [ ${URL: -3} == ".gz" ] && [ ! ${OUTPUT: -3} == ".gz" ]; then
    wget -q -O - "$URL" | gunzip -c > $TMP 2> ${snakemake_log}
else
    wget -q -O $TMP "$URL" 2> ${snakemake_log}
fi

if ! grep -v '#' $TMP | head -n 1 | grep -q 'chr' 2>> ${snakemake_log}; then
    echo "Adding \"chr\" to first column, then mv'ing to $OUTPUT"
    awk '{ if ("#"~$0) {print $0} else {print "chr"$0}}' $TMP > $OUTPUT 2>> ${snakemake_log}
else
    echo "Mv'ing to $OUTPUT"
    mv $TMP $OUTPUT 2>> ${snakemake_log}
fi
