#!/bin/bash

WD="../.."

LOG_FOLDER="$WD/var/log/slurm/$(date -I)"
mkdir -p "$LOG_FOLDER"

SBATCH="sbatch -J {cluster.jobname} --partition={cluster.partition} --mem={cluster.mem} --time={cluster.time} --cpus-per-task={cluster.cpus} -e ${LOG_FOLDER%/}/{cluster.log} -o ${LOG_FOLDER%/}/{cluster.log}"

CONFIGFILE="config.yaml"
CLUSTERCONFIG="$WD/config/cluster_config.json"

SINGULARITY_FOLDER="$WD/env/container"
SINGULARITY_MOUNTS="$(grep mounts config.yaml | sed 's/ //g' | awk -F":" '{print $2}')"
echo "$SINGULARITY_MOUNTS"

snakemake --jobs 100 \
  --reason \
  --printshellcmds \
  --show-failed-log \
  --configfile $CONFIGFILE \
  --use-singularity \
  --singularity-prefix $SINGULARITY_FOLDER \
  --singularity-args "--bind $SINGULARITY_MOUNTS" \
  --latency-wait 120 \
  --cluster "$SBATCH" \
  --cluster-config $CLUSTERCONFIG \
  --rerun-incomplete \
  |& tee "$LOG_FOLDER/snakemake-$(date '+%H%M%S').log"

