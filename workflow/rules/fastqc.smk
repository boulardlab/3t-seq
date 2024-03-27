

rule fastqc_raw:
    input:
        get_fastq,
    output:
        fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.zip"),
        report(
            fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.html"),
            category="FastQC",
            subcategory="Raw reads",
            labels={"serie": "{serie}", "sample": "{sample}"},
        ),
    params:
        fastqc_folder=lambda wildcards: os.path.join(fastqc_raw_folder, wildcards.serie),
    threads: 4
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc/{serie}/{sample}.log"),
    resources:
        runtime=60,
        mem_mb=4000,
    shell:
        """
        set -e 
        fastqc -t {threads} -noextract -o {params.fastqc_folder} {input}
        """


rule multiqc_raw:
    input:
        unpack(get_fastqc),
    output:
        report(
            multiqc_raw_folder.joinpath("{serie}", "multiqc_report.html"),
            category="MultiQC",
            subcategory="Raw reads",
            labels={"serie": "{serie}"},
        ),
    params:
        fastqc_folder=fastqc_raw_folder,
        multiqc_folder=multiqc_raw_folder,
    log:
        log_folder.joinpath("multiqc-raw", "multiqc-{serie}.log"),
    conda:
        "../env/qc.yml"
    resources:
        runtime=20,
        mem_mb=2048,
    shell:
        """
        set -e 
        multiqc --fullnames --dirs --export -f \
        -o {params.multiqc_folder}/{wildcards.serie} \
        {params.fastqc_folder}/{wildcards.serie} |& tee {log}
        """
