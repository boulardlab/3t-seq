#!/usr/bin/env bash

set -e
set -x 

sleep $(( $RANDOM % 10 + 2 ))

echo "Download GtRNAdb" 2>&1 > "${snakemake_log}"

mkdir -p "${snakemake_output}"
cd "${snakemake_output}"

F=$(basename "${snakemake_params[url]}")

wget --quiet "${snakemake_params[url]}"

tar xvf "$F"