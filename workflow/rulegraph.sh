#!/bin/bash

snakemake -j1 --configfile config.yaml --rulegraph | dot -Tsvg > rulegraph.svg
convert rulegraph.svg rulegraph.png