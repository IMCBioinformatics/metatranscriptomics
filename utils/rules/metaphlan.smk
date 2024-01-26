rule metaphlan:
    input:
        reads = config["output_dir"] + "/sortmerna/output/{sample}.fq"
    output:
        bt = config["output_dir"] + "/metaphlan/{sample}_bowtie2.bz2",
        pr = config["output_dir"] + "/metaphlan/{sample}_profile.txt"
    params: 
        metaphlan_database = config["metaphlan_database"],
	threads = config["threads"]
    conda: "metaphlan4.0.6"
    shell:
            "metaphlan -t rel_ab_w_read_stats --unclassified_estimation {input.reads} --input_type fastq "
            "--bowtie2db {params.metaphlan_database} --bowtie2out {output.bt} --nproc {params.threads} -o {output.pr}"



#giving resuls in count
rule sgb_to_GTDB:
    input:
         sg=config["output_dir"] + "/metaphlan/{sample}_profile.txt",
    params:
         sgb_to_gtdb_tsv_file=config["sgb_to_gtdb_tsv_file"],
         gtdb_output_dir=config["output_dir"] + "/metaphlan/GTDB"
    output:
         gtdb=config["output_dir"] + "/metaphlan/GTDB/{sample}_profile.txt"
    conda: "metaphlan4.0.6"
    shell:
        """
        mkdir -p {params.gtdb_output_dir}
        python utils/scripts/sgb_to_gtdb_profile_abundances.py --metaphlan_SGB_profile_infile {input.sg} --sgb_to_gtdb_tsv_file {params.sgb_to_gtdb_tsv_file} --output_dir {params.gtdb_output_dir}
        """


rule mergeprofiles:
    input: expand(config["output_dir"] + "/metaphlan/{sample}_profile.txt", sample=SAMPLES),
           expand(config["output_dir"] + "/metaphlan/GTDB/{sample}_profile.txt", sample=SAMPLES)
    params:
         metaphlan_SGB_profile_dir=config["output_dir"] + "/metaphlan",
         metaphlan_GTDB_profile_dir=config["output_dir"] + "/metaphlan/GTDB",
         output_dir=config["output_dir"]
    output: o1=config["output_dir"] + "/metaphlan/results/merged_abundance_table.txt",
            o2=config["output_dir"] + "/metaphlan/results/merged_abundance_table_species.txt",
            o3=config["output_dir"] + "/metaphlan/results/merged_abundance_table_GTDB.txt",
            o4=config["output_dir"] + "/metaphlan/results/merged_abundance_table_species_GTDB.txt"
    conda: "metaphlan4.0.6"
    shell: """
           python utils/scripts/merge_metaphlan_profiles_to_tables.py --metaphlan_SGB_profile_dir {params.metaphlan_SGB_profile_dir} --metaphlan_GTDB_profile_dir {params.metaphlan_GTDB_profile_dir} --output_dir {params.output_dir} 
           grep -E "(s__)|(^ID)|(clade_name)|(UNKNOWN)|(UNCLASSIFIED)" {output.o1} | grep -v "t__"  > {output.o2}
           grep -E "(s__)|(^ID)|(clade_name)|(UNKNOWN)|(UNCLASSIFIED)" {output.o3} | grep -v "t__"  > {output.o4}
           """
