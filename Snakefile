# ********************************
# * Snakefile for metqc pipeline *
# ********************************

# **** Variables ****
configfile: "config.yaml"

# **** Imports ****
import pandas as pd
import os

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# The forward read number for output files.
forward_read_num = config["forward_read_suffix"].split(".",1)[0]

# The reverse read number for output files.
reverse_read_num = config["reverse_read_suffix"].split(".",1)[0]


myoutput = list()

if config['Sortmerna_run'] == True:
    myoutput.append(expand(os.path.join(config["output_dir"], "sortmerna", "output", "{sample}.fq"), sample=SAMPLES))


# **** Rules ****

rule all:
    input:
        os.path.join(config["output_dir"],"metqc/multiqc","multiqc_report_raw.html"),
        os.path.join(config["output_dir"],"metqc/multiqc","multiqc_report_prinseq_filtered.html"),
        os.path.join(config["output_dir"],"metqc/multiqc","multiqc_report_bmtagger_filtered.html"),
	## comment out the line below if host_contamination rule fails.
        os.path.join(config["output_dir"],"metqc/seqkit","qc_seqkit.csv"),
        expand(config["output_dir"] + "/metaphlan/SGB/{sample}_bowtie2.bz2",sample=SAMPLES),
        expand(config["output_dir"] + "/metaphlan/SGB/{sample}_profile.txt",sample=SAMPLES),
        expand(config["output_dir"] + "/metaphlan/GTDB/{sample}_profile.txt",sample=SAMPLES),
        config["output_dir"] + "/metaphlan/results/merged_abundance_table_SGB.txt",
        config["output_dir"] + "/metaphlan/results/merged_abundance_table_SGB_species.txt",
        config["output_dir"] + "/metaphlan/results/merged_abundance_table_GTDB.txt",
        config["output_dir"] + "/metaphlan/results/merged_abundance_table_GTDB_species.txt",
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_uniref_renamed_relab_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_ko_renamed_relab_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_rxn_renamed_relab_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_eggnog_renamed_relab_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_pathabundance_relab_stratified.tsv", 
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_uniref_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_ko_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_rxn_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_eggnog_renamed_relab_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_relab/merged_pathabundance_relab_unstratified.tsv", 
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_uniref_renamed_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_ko_renamed_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_rxn_renamed_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_eggnog_renamed_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_pathabundance_stratified.tsv", 
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_uniref_renamed_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_ko_renamed_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_rxn_renamed_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_eggnog_renamed_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_raw/merged_pathabundance_unstratified.tsv", 
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_ko_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_stratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_pathabundance_cpm_stratified.tsv", 
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_ko_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_unstratified.tsv",
        config["output_dir"] + "/metannotate/final_results_cpm/merged_pathabundance_cpm_unstratified.tsv" 




include: "utils/rules/metqc.smk"
include: "utils/rules/sortmerna.smk"
include: "utils/rules/metaphlan.smk"
include: "utils/rules/metannotate.smk"
