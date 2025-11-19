#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=7-00:00:00
#SBATCH --mem=10G
#SBATCH --error=record_sbatch_run.%J.err
#SBATCH --output=record_sbatch_run.%J.out



time_in_seconds=$(date +%s)

log_dir="$(pwd)"
log_file="logs/pipeline-analysis_${time_in_seconds}.log.txt"

num_jobs=60

snakemake --unlock

conda activate snakemake

echo "started at: `date`"

snakemake --latency-wait 100 --rerun-incomplete --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --use-conda &>> $log_dir/$log_file

output_dir=$(grep "output_dir" < config.yaml | cut -d ' ' -f2 | sed 's/"//g')
list_files=$(grep "sampletable" < config.yaml | cut -d ' ' -f2 | sed 's/"//g')

#Copying all snakemake/log files of the run in the output folder
snakemake_file_dir="${output_dir}/snakemake_files"
mkdir -p $snakemake_file_dir
rm -rf "$snakemake_file_dir"/*

cp $list_files $snakemake_file_dir
cp Snakefile $snakemake_file_dir
cp config.yaml $snakemake_file_dir
cp cluster.json $snakemake_file_dir
cp pipeline_sbatch.sh $snakemake_file_dir 
cp -rf logs $snakemake_file_dir
cp -rf utils $snakemake_file_dir


echo "finished with exit code $? at: `date`"



