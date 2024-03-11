#!/usr/bin/env bash

set -e

URL="${snakemake_params[url]}"

if [[ $URL == *.gz ]]; then 
  TMP=$(mktemp -u --suffix .gz)
else
  TMP=$(mktemp -u)
fi

echo "Downloading to $TMP" | tee -a ${snakemake_log}

OUTPUT=${snakemake_output}
mkdir -pv $(dirname $OUTPUT)

# wget -O $TMP "$URL"
curl "$URL" \
-H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/118.0' \
-H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' \
-H 'Accept-Language: en-US,en;q=0.5' \
-H 'Accept-Encoding: gzip, deflate' \
-H 'Connection: keep-alive' \
-H 'Upgrade-Insecure-Requests: 1' \
-H 'DNT: 1' \
-H 'Sec-GPC: 1' \
-H 'Pragma: no-cache' \
-H 'Cache-Control: no-cache' \
--silent \
--output $TMP

sleep $(( $RANDOM % 10 + 2 ))

if [[ $URL == *.gz ]] && [[ ! $OUTPUT == *.gz ]]; then
    gunzip $TMP 
    sleep $(( $RANDOM % 10 + 2 ))
fi

if grep -v '#' "${TMP%.gz}" | head -n 1 | grep -q '^chr' | tee -a ${snakemake_log}; then
    echo "Mv'ing to $OUTPUT" | tee -a ${snakemake_log}
    mv $TMP $OUTPUT
else
    echo "Adding \"chr\" to first column, then move to $OUTPUT" | tee -a ${snakemake_log}
    awk -F "\t" -v OFS="\t" '!/^#/{print "chr"$0}/#/{print}' ${TMP%.gz} > $OUTPUT 
fi
