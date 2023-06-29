rule filter_bam:
    input:
        alignment=starTE_folder.joinpath("{serie}/{method}/{sample}.Aligned.out.bam"),
        annotation=rmsk_bed,
    output:
        starTE_folder.joinpath("{serie}/filter/{method}/{sample}.TEonly.bam"),
    log:
        log_folder.joinpath("samtools_view/{serie}/{method}/{sample}.log"),
    threads: 2
    conda:
        "../env/samtools.yml"
    shell:
        """
        set -x
        samtools view -L {input.annotation} -@ {threads} -o {output} {input.alignment} |& tee {log}
        """
