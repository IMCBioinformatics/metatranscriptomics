rule metaphlan:
    input:
        reads = config["output_dir"] + "/sortmerna/output/{sample}.fq" if config.get("Sortmerna_run", True) else rules.merge_reads.output
    output:
        bt = config["output_dir"] + "/metaphlan/SGB/{sample}_bowtie2.bz2",
        pr = config["output_dir"] + "/metaphlan/SGB/{sample}_profile.txt",
        sam = config["output_dir"] + "/metaphlan/SGB/{sample}_sam.bz2"
    params:
        metaphlan_database = config["metaphlan_database"],
        index = config["index"],
        threads = config["metaphlan_cpus"]
    conda: config["metaphlan_version"]
    shell:
        "metaphlan -t rel_ab_w_read_stats --unclassified_estimation {input.reads} --input_type fastq -s {output.sam} "
        "--bowtie2db {params.metaphlan_database} --index {params.index} --bowtie2out {output.bt} --nproc {params.threads} -o {output.pr}"



#giving resuls in count
rule sgb_to_GTDB_count:
    input:
         sg=config["output_dir"] + "/metaphlan/SGB/{sample}_profile.txt",
    params:
         sgb_to_gtdb_tsv_file=config["sgb_to_gtdb_tsv_file"],
         gtdb_output_dir=config["output_dir"] + "/metaphlan/GTDB"
    output:
         gtdb=config["output_dir"] + "/metaphlan/GTDB/{sample}_profile.txt"
    conda: config["metaphlan_version"]
    shell:
        """
        mkdir -p {params.gtdb_output_dir}
        python utils/scripts/sgb_to_gtdb_profile_abundances.py --metaphlan_SGB_profile_infile {input.sg} --sgb_to_gtdb_tsv_file {params.sgb_to_gtdb_tsv_file} --output_dir {params.gtdb_output_dir}
        """


rule sgb_to_GTDB_relab:
    input:
        sg=config["output_dir"] + "/metaphlan/SGB/{sample}_profile.txt"
    params:
        gtdb_output_dir=config["output_dir"] + "/metaphlan/GTDB_relab"
    output:
        gtdb=config["output_dir"] + "/metaphlan/GTDB_relab/{sample}_profile.txt"
    conda: "utils/envs/metaphlan4.yaml"
    shell:
        "sgb_to_gtdb_profile.py  -i {input.sg} -o {output.gtdb}"


rule mergeprofiles:
    input: expand(config["output_dir"] + "/metaphlan/SGB/{sample}_profile.txt", sample=SAMPLES),
           expand(config["output_dir"] + "/metaphlan/GTDB/{sample}_profile.txt", sample=SAMPLES)
    params:
         metaphlan_SGB_profile_dir=config["output_dir"] + "/metaphlan/SGB",
         metaphlan_GTDB_profile_dir=config["output_dir"] + "/metaphlan/GTDB",
         output_dir=config["output_dir"]
    output: o1=config["output_dir"] + "/metaphlan/results/merged_abundance_table_SGB.txt",
            o2=config["output_dir"] + "/metaphlan/results/merged_abundance_table_SGB_species.txt",
            o3=config["output_dir"] + "/metaphlan/results/merged_abundance_table_GTDB.txt",
            o4=config["output_dir"] + "/metaphlan/results/merged_abundance_table_GTDB_species.txt"
    conda: config["metaphlan_version"]
    shell: """
           python utils/scripts/merge_metaphlan_profiles_to_tables.py --metaphlan_SGB_profile_dir {params.metaphlan_SGB_profile_dir} --metaphlan_GTDB_profile_dir {params.metaphlan_GTDB_profile_dir} --output_dir {params.output_dir} 
           grep -E "(s__)|(^ID)|(clade_name)|(UNKNOWN)|(UNCLASSIFIED)" {output.o1} | grep -v "t__"  > {output.o2}
           grep -E "(s__)|(^ID)|(clade_name)|(UNKNOWN)|(UNCLASSIFIED)" {output.o3} | grep -v "t__"  > {output.o4}
           """

rule mergeprofiles_SGB_relab:
    input: expand(config["output_dir"] + "/metaphlan/SGB/{sample}_profile.txt", sample=SAMPLES)
    output: o1=config["output_dir"] + "/metaphlan/results/merged_abundance_table_SGB_relab.txt",
            o2=config["output_dir"] + "/metaphlan/results/merged_abundance_table_SGB_species_relab.txt"
    params: profiles=config["output_dir"]+"/metaphlan/SGB/*_profile.txt"
    conda: config["metaphlan_version"]
    shell: """
           python utils/scripts/merge_metaphlan_tables_relab.py {params.profiles} > {output.o1}
           grep -E "(s__)|(^ID)|(clade_name)|(UNKNOWN)|(UNCLASSIFIED)" {output.o1} | grep -v "t__"  > {output.o2}
           """


rule mergeprofiles_GTDB_relab:
    input: expand(config["output_dir"] + "/metaphlan/GTDB_relab/{sample}_profile.txt", sample=SAMPLES)
    output: o1=config["output_dir"] + "/metaphlan/results/merged_abundance_table_GTDB_relab.txt",
            o2=config["output_dir"] + "/metaphlan/results/merged_abundance_table_GTDB_species_relab.txt"
    params: profiles=config["output_dir"]+"/metaphlan/GTDB_relab/*_profile.txt"
    conda: config["metaphlan_version"]
    shell: """
           python utils/scripts/merge_metaphlan_tables_relab.py {params.profiles} > {output.o1}
           grep -E "(s__)|(^ID)|(clade_name)|(UNKNOWN)|(UNCLASSIFIED)" {output.o1} | grep -v "t__"  > {output.o2}
           """


