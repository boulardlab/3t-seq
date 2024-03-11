#!/usr/bin/env bash

set -euo pipefail

echo "Download GtRNAdb" 2>&1 > "${snakemake_log}"

mkdir -pv "${snakemake_params[output_dir]}"
cd "${snakemake_params[output_dir]}"

F=$(basename "${snakemake_params[url]}")

wget --quiet "${snakemake_params[url]}"

tar xvf "$F"
sleep $(( $RANDOM % 10 + 2 ))