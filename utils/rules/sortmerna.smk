rule sortmerna:
    input:
        reads = config["output_dir"] + "/sortmerna/merged_data/{sample}.fastq"
    output:
        bt = config["output_dir"] + "/sortmerna/output/{sample}.fq",
    params: 
        ref1=config["ref1"],
        ref2=config["ref2"],
        ref3=config["ref3"],
        ref4=config["ref4"],
        ref5=config["ref5"],
        ref6=config["ref6"],
	threads = config["threads"],
        dir=config["output_dir"] + "/sortmerna/{sample}",
        unaligned=config["output_dir"] + "/sortmerna/output/{sample}",
    conda: "sortmerna"
    shell:
         """
        if [[ "{config[Sortmerna_run]}" == "True" ]]; then
            sortmerna --ref {params.ref1} --ref {params.ref2} --ref {params.ref3} --ref {params.ref4} --ref {params.ref5} --ref {params.ref6} \
                --reads {input.reads} --threads {params.threads}  --workdir {params.dir} --other {params.unaligned} --fastx
        else
            echo "Rule 'Sortmerna' is not executed because this is a metagenomics run and 'Sortmerna_run' is set to 'false' in the config file."
        fi
        """
