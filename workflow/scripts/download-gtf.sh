#!/usr/bin/env bash

set -e

URL="${snakemake_params[url]}"

mkdir -p $(dirname ${snakemake_log})
touch ${snakemake_log}

if [[ $URL == *.gz ]]; then 
  TMP=$(mktemp -u --suffix .gz)
  echo "$URL refers to gzipped file. Temp file: $TMP" >> ${snakemake_log}
else
  TMP=$(mktemp -u)
  echo "$URL does not refer to gzipped file. Temp file: $TMP" >> ${snakemake_log}
fi

echo "Downloading to $TMP" >> ${snakemake_log}

OUTPUT=${snakemake_output}
mkdir -pv $(dirname $OUTPUT)

wget "$URL" \
--user-agent=' Mozilla/5.0 (X11; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/118.0' \
--header='Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' \
--header='Accept-Language: en-US,en;q=0.5' \
--header='Accept-Encoding: gzip, deflate' \
--header='Connection: keep-alive' \
--header='Upgrade-Insecure-Requests: 1' \
--header='DNT: 1' \
--header='Sec-GPC: 1' \
--header='Pragma: no-cache' \
--header='Cache-Control: no-cache' \
--quiet \
--output-document="$TMP"

sleep $(( $RANDOM % 10 + 2 ))

if [[ $URL == *.gz ]] && [[ ! $OUTPUT == *.gz ]]; then
    gunzip $TMP 
    sleep $(( $RANDOM % 10 + 2 ))
fi

if grep -v '#' "${TMP%.gz}" | head -n 1 | grep -q '^chr' >> ${snakemake_log}; then
    echo "Mv'ing to $OUTPUT" >> ${snakemake_log}
    mv $TMP $OUTPUT
else
    echo "Adding \"chr\" to first column, then move to $OUTPUT" >> ${snakemake_log}
    awk -F "\t" -v OFS="\t" '!/^#/{print "chr"$0}/#/{print}' ${TMP%.gz} > $OUTPUT 
fi
