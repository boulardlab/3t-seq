#!/bin/bash
set -e

N_RUNNING_JOBS=$(squeue -u $USER | grep RUNNING | wc -l)
[[ $N_RUNNING_JOBS -gt 0 ]] && ( echo -n "STOPPING JOBS..." && scancel -u tabaro && echo -e "\tDONE" && sleep 3 ) || echo "No running jobs"
[[ -d ../../tmp ]] && rm -r ../../tmp/  || echo "No tmp folder"

git checkout tabaro/hpc
git pull --no-edit origin tabaro/main
snakemake -j1 --configfile config.yaml --unlock

echo "Run pipeline?" && read
./run.sh


