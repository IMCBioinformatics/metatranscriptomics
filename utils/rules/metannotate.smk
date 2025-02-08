rule humann3:
    input:
        rules.merge_reads.output,
        rules.metaphlan.output
    output:
        genefam = config["output_dir"] + "/metannotate/raw/{sample}_genefamilies.tsv",
        pathcov = config["output_dir"] + "/metannotate/raw/{sample}_pathcoverage.tsv",
        pathabun = config["output_dir"] + "/metannotate/raw/{sample}_pathabundance.tsv"
    params:
        threads = config["metannotate_cpus"],
        db2 = config["output_dir"] + "/metannotate/databases/{sample}_database",
        metaphlan= config["path"] + "/output/metaphlan/SGB/{sample}_profile.txt",
        s = "{sample}",
        output=config["output_dir"] + "/metannotate/raw"
    conda: config["humann_version"]
    shell:
            """
            humann3 --input {input[0]} --threads {params.threads} --output {params.output} --output-basename {params.s} --nucleotide-database {config[nuc_db]} --protein-database {config[prot_db]} --taxonomic-profile {params.metaphlan}
# This has been taken out becaue we can use the metaphlan outputs from the metaphlan run and don't need to re-do them.
#--metaphlan-options="{params.db1}"
            rm -rf {params.db2}
            """

rule regroup:
     input: genefam=config["output_dir"] + "/metannotate/raw/{sample}_genefamilies.tsv"
     output:genefam_eggnog=config["output_dir"] + "/metannotate/regrouped/{sample}_genefamilies_eggnog.tsv",
            genefam_ko=config["output_dir"] + "/metannotate/regrouped/{sample}_genefamilies_ko.tsv",
            genefam_rxn=config["output_dir"] + "/metannotate/regrouped/{sample}_genefamilies_rxn.tsv"
     params: eggnog=config["eggnog_uniref90"],
             ko=config["ko_uniref90"]
     conda: config["humann_version"]
     shell:
       """ humann_regroup_table --input {input.genefam}     --output {output.genefam_ko}  -c {params.ko} ; 
           humann_regroup_table --input {input.genefam}     --output {output.genefam_eggnog}  -c {params.eggnog};
           humann_regroup_table --input {input.genefam}     --output {output.genefam_rxn}  -g uniref90_rxn;"""


rule rename:
      input:genefam=config["output_dir"] + "/metannotate/raw/{sample}_genefamilies.tsv",
            genefam_eggnog=config["output_dir"] + "/metannotate/regrouped/{sample}_genefamilies_eggnog.tsv",
            genefam_ko=config["output_dir"] + "/metannotate/regrouped/{sample}_genefamilies_ko.tsv",
            genefam_rxn=config["output_dir"] + "/metannotate/regrouped/{sample}_genefamilies_rxn.tsv"

      output: 
            genefam_name=config["output_dir"] + "/metannotate/renamed/uniref/{sample}_genefamilies_renamed.tsv",
            genefam_ko_name=config["output_dir"] + "/metannotate/renamed/ko/{sample}_genefamilies_ko_renamed.tsv",
            genefam_rxn_name=config["output_dir"] + "/metannotate/renamed/rxn/{sample}_genefamilies_rxn_renamed.tsv",
            genefam_eggnog_name=config["output_dir"] + "/metannotate/renamed/eggnog/{sample}_genefamilies_eggnog_renamed.tsv"
      conda: config["humann_version"]
      params: uniref=config["uniref90_name"], 
              eggnog=config["eggnog_name"],
              ko=config["ko_name"]
      shell:
       """ 
	humann_rename_table -i {input.genefam} -c {params.uniref} -o {output.genefam_name};
	humann_rename_table -i {input.genefam_ko} -c {params.ko} -o {output.genefam_ko_name};
	humann_rename_table -i {input.genefam_eggnog} -c {params.eggnog} -o {output.genefam_eggnog_name};
	humann_rename_table -i {input.genefam_rxn} -n metacyc-rxn -o {output.genefam_rxn_name} """
     


rule merge_output:
    input:
            genefam_name=expand(config["output_dir"] + "/metannotate/renamed/uniref/{sample}_genefamilies_renamed.tsv",sample=SAMPLES),
            genefam_eggnog_name=expand(config["output_dir"] + "/metannotate/renamed/eggnog/{sample}_genefamilies_eggnog_renamed.tsv",sample=SAMPLES),
            genefam_ko_name=expand(config["output_dir"] + "/metannotate/renamed/ko/{sample}_genefamilies_ko_renamed.tsv",sample=SAMPLES),
            genefam_rxn_name=expand(config["output_dir"] + "/metannotate/renamed/rxn/{sample}_genefamilies_rxn_renamed.tsv",sample=SAMPLES),
            pathabun = expand(config["output_dir"] + "/metannotate/raw/{sample}_pathabundance.tsv", sample=SAMPLES),
            pathcov = expand(config["output_dir"] + "/metannotate/raw/{sample}_pathcoverage.tsv", sample=SAMPLES)

    output:
        genefam_uniref = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_eggnog = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed.tsv",
        genefam_ko = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed.tsv",       
        pathabun = config["output_dir"] + "/metannotate/merged/merged_pathabundance.tsv",
        pathcov = config["output_dir"] + "/metannotate/merged/merged_pathcoverage.tsv"
    conda: config["humann_version"]
    params:
        dir= config["output_dir"] 
    shell:
         """
            humann_join_tables --input {params.dir}/metannotate/renamed/uniref/ --output {output.genefam_uniref} --file_name genefamilies; 
            humann_join_tables --input {params.dir}/metannotate/renamed/ko/ --output {output.genefam_ko} --file_name genefamilies; 
            humann_join_tables --input {params.dir}/metannotate/renamed/rxn/ --output {output.genefam_rxn} --file_name genefamilies; 
            humann_join_tables --input {params.dir}/metannotate/renamed/eggnog/ --output {output.genefam_eggnog} --file_name genefamilies; 
            humann_join_tables --input {params.dir}/metannotate/raw/ --output {output.pathcov} --file_name pathcoverage; 
            humann_join_tables --input {params.dir}/metannotate/raw/  --output {output.pathabun} --file_name _pathabundance; """


rule relab:
    input:
        genefam_uniref = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_ko = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed.tsv",
        genefam_eggnog = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed.tsv",
        pathabun = config["output_dir"] + "/metannotate/merged/merged_pathabundance.tsv"
    output:
        genefam_uniref_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed_relab.tsv",
        genefam_ko_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed_relab.tsv",
        genefam_rxn_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed_relab.tsv",
        genefam_eggnog_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed_relab.tsv",
        pathabun_relab = config["output_dir"] + "/metannotate/merged/merged_pathabundance_relab.tsv" 
    conda: config["humann_version"]
    shell:
         """
         humann_renorm_table --input {input.genefam_uniref} --output {output.genefam_uniref_relab} -s n --units relab;
         humann_renorm_table --input {input.genefam_ko} --output {output.genefam_ko_relab} -s n --units relab;
         humann_renorm_table --input {input.genefam_rxn} --output {output.genefam_rxn_relab} -s n --units relab;
         humann_renorm_table --input {input.genefam_eggnog} --output {output.genefam_eggnog_relab} -s n --units relab;
         humann_renorm_table --input {input.pathabun} --output {output.pathabun_relab} -s n --units relab; """


rule cpm:
    input:
        genefam_uniref = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_ko = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed.tsv",
        genefam_eggnog = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed.tsv",
        pathabun = config["output_dir"] + "/metannotate/merged/merged_pathabundance.tsv"
    output:
        genefam_uniref_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed_cpm.tsv",
        genefam_ko_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed_cpm.tsv",
        genefam_rxn_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed_cpm.tsv",
        genefam_eggnog_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed_cpm.tsv",
        pathabun_cpm = config["output_dir"] + "/metannotate/merged/merged_pathabundance_cpm.tsv" 
    conda: config["humann_version"] 
    shell:
         """
         humann_renorm_table --input {input.genefam_uniref} --output {output.genefam_uniref_cpm} -s n --units cpm;
         humann_renorm_table --input {input.genefam_ko} --output {output.genefam_ko_cpm} -s n --units cpm;
         humann_renorm_table --input {input.genefam_rxn} --output {output.genefam_rxn_cpm} -s n --units cpm;
         humann_renorm_table --input {input.genefam_eggnog} --output {output.genefam_eggnog_cpm} -s n --units cpm;
         humann_renorm_table --input {input.pathabun} --output {output.pathabun_cpm} -s n --units cpm; """


rule final_results_raw:
    input: 
        genefam_uniref_raw = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed.tsv",
        genefam_ko_raw = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed.tsv",
        genefam_rxn_raw = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed.tsv",
        genefam_eggnog_raw = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed.tsv",
        pathabun_raw = config["output_dir"] + "/metannotate/merged/merged_pathabundance.tsv" 
    params: 
        output_dir=config["output_dir"] + "/metannotate/final_results_raw"
    output:
        genefam_uniref_raw_s = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_uniref_renamed_stratified.tsv",
        genefam_ko_raw_s = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_ko_renamed_stratified.tsv",
        genefam_rxn_raw_s = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_rxn_renamed_stratified.tsv",
        genefam_eggnog_raw_s = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_eggnog_renamed_stratified.tsv",
        pathabun_raw_s = config["output_dir"] + "/metannotate/final_results_raw/merged_pathabundance_stratified.tsv", 
        genefam_uniref_raw_us = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_uniref_renamed_unstratified.tsv",
        genefam_ko_raw_us = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_ko_renamed_unstratified.tsv",
        genefam_rxn_raw_us = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_rxn_renamed_unstratified.tsv",
        genefam_eggnog_raw_us = config["output_dir"] + "/metannotate/final_results_raw/merged_genefamilies_eggnog_renamed_unstratified.tsv",
        pathabun_raw_us = config["output_dir"] + "/metannotate/final_results_raw/merged_pathabundance_unstratified.tsv" 

    conda: config["humann_version"]
    shell:
        """ 
            humann_split_stratified_table --input {input.genefam_uniref_raw} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_ko_raw} --output  {params.output_dir};
            humann_split_stratified_table --input {input.genefam_rxn_raw} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_eggnog_raw} --output {params.output_dir};
            humann_split_stratified_table --input {input.pathabun_raw} --output {params.output_dir};

         """

rule final_results_relab:
    input: 
        genefam_uniref_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed_relab.tsv",
        genefam_ko_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed_relab.tsv",
        genefam_rxn_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed_relab.tsv",
        genefam_eggnog_relab = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed_relab.tsv",
        pathabun_relab = config["output_dir"] + "/metannotate/merged/merged_pathabundance_relab.tsv" 
    
    output:
        genefam_uniref_relab_s = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_uniref_renamed_relab_stratified.tsv",
        genefam_ko_relab_s = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_ko_renamed_relab_stratified.tsv",
        genefam_rxn_relab_s = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_rxn_renamed_relab_stratified.tsv",
        genefam_eggnog_relab_s = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_eggnog_renamed_relab_stratified.tsv",
        pathabun_relab_s = config["output_dir"] + "/metannotate/final_results_relab/merged_pathabundance_relab_stratified.tsv", 
        genefam_uniref_relab_us = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_uniref_renamed_relab_unstratified.tsv",
        genefam_ko_relab_us = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_ko_renamed_relab_unstratified.tsv",
        genefam_rxn_relab_us = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_rxn_renamed_relab_unstratified.tsv",
        genefam_eggnog_relab_us = config["output_dir"] + "/metannotate/final_results_relab/merged_genefamilies_eggnog_renamed_relab_unstratified.tsv",
        pathabun_relab_us = config["output_dir"] + "/metannotate/final_results_relab/merged_pathabundance_relab_unstratified.tsv" 

    params: 
        output_dir=config["output_dir"] + "/metannotate/final_results_relab"

    conda: config["humann_version"]
    shell:
        """ 
            humann_split_stratified_table --input {input.genefam_uniref_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_ko_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_rxn_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_eggnog_relab} --output {params.output_dir};
            humann_split_stratified_table --input {input.pathabun_relab} --output {params.output_dir};
 
         """


rule final_results_cpm:
    input: 
        genefam_uniref_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_uniref_renamed_cpm.tsv",
        genefam_ko_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_ko_renamed_cpm.tsv",
        genefam_rxn_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_rxn_renamed_cpm.tsv",
        genefam_eggnog_cpm = config["output_dir"] + "/metannotate/merged/merged_genefamilies_eggnog_renamed_cpm.tsv",
        pathabun_cpm = config["output_dir"] + "/metannotate/merged/merged_pathabundance_cpm.tsv" 
    
    output:
        genefam_uniref_cpm_s = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_stratified.tsv",
        genefam_ko_cpm_s = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_ko_renamed_cpm_stratified.tsv",
        genefam_rxn_cpm_s = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_stratified.tsv",
        genefam_eggnog_cpm_s = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_stratified.tsv",
        pathabun_cpm_s = config["output_dir"] + "/metannotate/final_results_cpm/merged_pathabundance_cpm_stratified.tsv", 
        genefam_uniref_cpm_us = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_uniref_renamed_cpm_unstratified.tsv",
        genefam_ko_cpm_us = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_ko_renamed_cpm_unstratified.tsv",
        genefam_rxn_cpm_us = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_rxn_renamed_cpm_unstratified.tsv",
        genefam_eggnog_cpm_us = config["output_dir"] + "/metannotate/final_results_cpm/merged_genefamilies_eggnog_renamed_cpm_unstratified.tsv",
        pathabun_cpm_us = config["output_dir"] + "/metannotate/final_results_cpm/merged_pathabundance_cpm_unstratified.tsv" 
    params: 
        output_dir=config["output_dir"] + "/metannotate/final_results_cpm"

    conda: config["humann_version"]
    shell:
        """ 
            humann_split_stratified_table --input {input.genefam_uniref_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_ko_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_rxn_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.genefam_eggnog_cpm} --output {params.output_dir};
            humann_split_stratified_table --input {input.pathabun_cpm} --output {params.output_dir};
         """

