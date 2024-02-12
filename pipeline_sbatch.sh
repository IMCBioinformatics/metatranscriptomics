#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=10G
#SBATCH --error=record_sbatch_run.%J.err
#SBATCH --output=record_sbatch_run.%J.out

log_dir="$(pwd)"
log_file="logs/pipeline-analysis.log.txt"
num_jobs=60

snakemake --unlock

echo "started at: `date`"

snakemake --latency-wait 100 --rerun-incomplete --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --use-conda &>> $log_dir/$log_file

echo "finished with exit code $? at: `date`"


