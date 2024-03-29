rule filter_bam:
    input:
        alignment=starTE_folder.joinpath("{serie}/{method}/{sample}.Aligned.out.bam"),
        annotation=rmsk_folder.joinpath(
            "{0}.{1}".format(config["genome"]["label"], "bed")
        ),
    output:
        starTE_folder.joinpath("{serie}/filter/{method}/{sample}.TEonly.bam"),
    log:
        log_folder.joinpath("samtools_view/{serie}/{method}/{sample}.log"),
    threads: 4
    resources:
        runtime=30,
        mem_mb=16000,
    conda:
        "../env/samtools.yml"
    shell:
        """
        set -e 
        samtools view -L {input.annotation} -@ {threads} -o {output} {input.alignment} |& tee {log}
        """
