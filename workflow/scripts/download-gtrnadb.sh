#!/usr/bin/env bash

set -euo pipefail
set -x 

# Redirect all output and stderr to the snakemake log file
exec &> "${snakemake_log[0]}"

echo "Starting GtRNAdb download script..."

OUTPUT_DIR="${snakemake_params[output_dir]}"
URL="${snakemake_params[url]}"

echo "Creating output directory: $OUTPUT_DIR"
mkdir -pv "$OUTPUT_DIR"

F=$(basename "$URL")
TARGET_FILE="$OUTPUT_DIR/$F"

echo "Downloading GtRNAdb from $URL to $TARGET_FILE..."
wget -O "$TARGET_FILE" "$URL"

echo "Extracting $TARGET_FILE into $OUTPUT_DIR..."
tar -xvf "$TARGET_FILE" -C "$OUTPUT_DIR"

echo "Cleaning up downloaded archive $TARGET_FILE..."
rm "$TARGET_FILE"

echo "GtRNAdb download and extraction completed successfully."