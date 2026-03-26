#!/usr/bin/env bash

python /opt/SalmonTE/SalmonTE.py test \
	--inpath=${snakemake_input[infolder]} \
	--outpath=${snakemake_output} \
	--tabletype=csv \
	--figtype=png \
	--analysis_type=DE \
	--conditions=control,treatment |& \
	tee ${snakemake_log}
