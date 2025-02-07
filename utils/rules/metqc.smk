rule fastqc_raw:
    input:
        r1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        r1 = os.path.join(config["output_dir"],"metqc/fastqc_raw","{sample}"+forward_read_num+"_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"metqc/fastqc_raw","{sample}"+reverse_read_num+"_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"metqc/fastqc_raw")
    conda: "QC"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_raw:
    input:
        r1 = expand(os.path.join(config["output_dir"],"metqc/fastqc_raw","{sample}"+forward_read_num+"_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"metqc/fastqc_raw","{sample}"+reverse_read_num+"_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"metqc/multiqc","multiqc_report_raw.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"metqc/fastqc_raw/"),
	multiqc_dir = os.path.join(config["output_dir"],"metqc/multiqc/")
    conda: "QC"
    shell: "multiqc -c utils/envs/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_raw.html"

rule cutadapt:
    input:
        r1 = os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    output:
        r1 = os.path.join(config["output_dir"],"metqc/cutadapt","{sample}_r1_trimmed.fastq"),
        r2 = os.path.join(config["output_dir"],"metqc/cutadapt","{sample}_r2_trimmed.fastq")
    conda: "QC"
    shell:
            "cutadapt -m {config[minlength]} --max-n {config[maxn]} -a {config[fwd_adapter]} -A {config[rev_adapter]} "
            "-j {config[num_cpus]} -o {output.r1} -p {output.r2} "
            "{input.r1} {input.r2}"

rule prinseq:
    input:
        r1 = os.path.join(config["output_dir"],"metqc/cutadapt","{sample}_r1_trimmed.fastq") if config["run_cutadapt"] else os.path.join(config["input_dir"],"{sample}"+config["forward_read_suffix"]),
        r2 = os.path.join(config["output_dir"],"metqc/cutadapt","{sample}_r2_trimmed.fastq") if config["run_cutadapt"] else os.path.join(config["input_dir"],"{sample}"+config["reverse_read_suffix"])
    params:
        prefix = os.path.join(config["output_dir"],"metqc/prinseq","{sample}_filtered")
    output:
        r1 = os.path.join(config["output_dir"],"metqc/prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"metqc/prinseq","{sample}_filtered_2.fastq")
    conda: "../envs/prinseq_env.yaml"
    shell:
            "perl utils/scripts/prinseq-lite.pl -fastq {input.r1} -fastq2 {input.r2} "
            "-trim_left {config[trimleft]} -trim_right {config[trimright]} "
            "-out_good {params.prefix} -out_bad null -lc_method {config[lc_method]} -lc_threshold {config[lc_threshold]} "
            "-derep 1 -trim_qual_type {config[trim_qual_type]} -trim_qual_window "
            "{config[trim_qual_window]} -trim_qual_step {config[trim_qual_step]} "
            "-trim_qual_rule {config[trim_qual_rule]} -trim_qual_left {config[trim_qual_left]} "
            "-trim_qual_right {config[trim_qual_right]} -min_len {config[minlength]} "
			"-ns_max_n {config[maxn]}"

rule bmtagger:
    input:
        r1 = os.path.join(config["output_dir"],"metqc/prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"metqc/prinseq","{sample}_filtered_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_1.fastq"),
        r2 = os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_2.fastq"),
    params:
        n = os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged"),
        r3 = os.path.join(config["output_dir"],"metqc/bmtagger","bmtagger_complete.txt")
    conda: "../envs/bmtagger_env.yaml"
    shell:
        "bmtagger.sh -b {config[bmfilter_ref]} -x {config[srprism_ref]} -q 1 -1 {input.r1} -2 {input.r2} -o {params.n} -X;"
        " touch  {params.r3}"


rule fastqc_prinseq_filt:
    input:
        r1 = os.path.join(config["output_dir"],"metqc/prinseq","{sample}_filtered_1.fastq"),
        r2 = os.path.join(config["output_dir"],"metqc/prinseq","{sample}_filtered_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"metqc/prinseq","fastqc","{sample}_filtered_1_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"metqc/prinseq","fastqc","{sample}_filtered_2_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"metqc/prinseq","fastqc/")
    conda: "QC"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"

rule multiqc_prinseq_filt:
    input:
        r1 = expand(os.path.join(config["output_dir"],"metqc/prinseq","fastqc","{sample}_filtered_1_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"metqc/prinseq","fastqc","{sample}_filtered_2_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"metqc/multiqc","multiqc_report_prinseq_filtered.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"metqc/prinseq","fastqc/"),
        multiqc_dir = os.path.join(config["output_dir"],"metqc/multiqc/")
    conda: "QC"
    shell: "multiqc -c utils/envs/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_prinseq_filtered.html"



rule fastqc_bmtagger_filt:
    input:
        r1 = os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_1.fastq"),
        r2 = os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_2.fastq")
    output:
        r1 = os.path.join(config["output_dir"],"metqc/bmtagger","fastqc","{sample}_bmtagged_1_fastqc.html"),
        r2 = os.path.join(config["output_dir"],"metqc/bmtagger","fastqc","{sample}_bmtagged_2_fastqc.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"metqc/bmtagger","fastqc/")
    conda: "QC"
    shell: "fastqc -o {params.fastqc_dir} {input.r1} {input.r2}"


rule multiqc_bmtagger_filt:
    input:
        r1 = expand(os.path.join(config["output_dir"],"metqc/bmtagger","fastqc","{sample}_bmtagged_1_fastqc.html"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"metqc/bmtagger","fastqc","{sample}_bmtagged_2_fastqc.html"), sample=SAMPLES)
    output: os.path.join(config["output_dir"],"metqc/multiqc","multiqc_report_bmtagger_filtered.html")
    params:
        fastqc_dir = os.path.join(config["output_dir"],"metqc/bmtagger","fastqc/"),
        multiqc_dir = os.path.join(config["output_dir"],"metqc/multiqc/")
    conda: "QC"
    shell: "multiqc -c utils/envs/multiqc_config.yaml -f {params.fastqc_dir} -o {params.multiqc_dir} -n multiqc_report_bmtagger_filtered.html"
	

rule seqkit:
     input:
        r1 = expand(os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_1.fastq"), sample=SAMPLES),
        r2 = expand(os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_2.fastq"), sample=SAMPLES),
     output:
        raw=config["output_dir"]+"/metqc/seqkit/seq_kit_raw.csv",
        prinseq=config["output_dir"] +"/metqc/seqkit/seq_kit_prinseq.csv",
        bmtagger=config["output_dir"]+"/metqc/seqkit/seq_kit_bmtagger.csv",
        complete=config["output_dir"]+"/metqc/seqkit/seq_kit_complete.csv"
     params:
        raw=config["input_dir"],
        prinseq=directory(os.path.join(config["output_dir"],"metqc/prinseq")),
        bmtagger=directory(os.path.join(config["output_dir"],"metqc/bmtagger"))
     conda: "../envs/seqkit.yaml"
     shell:
         "seqkit stats -j {config[num_cpus]} {params.prinseq}/*_[0-9].fastq -o {output.prinseq};"
         "seqkit stats -j {config[num_cpus]} {params.bmtagger}/*.fastq -o {output.bmtagger};"
         "seqkit stats -j {config[num_cpus]} {params.raw}/*.fastq.gz -o {output.raw};"
         "touch {output.complete};"


rule host_contamination:
     input:
        raw=config["output_dir"]+"/metqc/seqkit/seq_kit_raw.csv",
        prinseq=config["output_dir"] +"/metqc/seqkit/seq_kit_prinseq.csv",
        bmtagger=config["output_dir"]+"/metqc/seqkit/seq_kit_bmtagger.csv",
        complete=config["output_dir"]+"/metqc/seqkit/seq_kit_complete.csv"
     params:
        r1=forward_read_num, #config["reverse_read_suffix"],
        r2=reverse_read_num #config["forward_read_suffix"]
     conda: "../envs/python3_8.yaml"
     output:
         hc=config["output_dir"]+"/metqc/seqkit/qc_seqkit.csv"
     script:"../scripts/host_contamination.py"


rule merge_reads:
    input:
        r1 = os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_1.fastq"),
        r2 = os.path.join(config["output_dir"],"metqc/bmtagger","{sample}_bmtagged_2.fastq")
    output:
        config["output_dir"] + "/metqc/merged_data/{sample}.fastq"
    shell:
        "cat {input.r1} {input.r2} > {output}"
