rule coverage:
    input:
        bam=star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
        bai=star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam.bai"),
    output:
        star_folder.joinpath("{serie}", "{sample}.bw"),
    log:
        log_folder.joinpath("bamCoverage_se-{serie}-{sample}.log"),
    conda:
        "../env/deeptools.yml"
    threads: 2
    resources:
        runtime=120,
        mem_mb=16000,
    params:
        others=get_bamcoverage_params,
    shell:
        """

        bamCoverage -b {input.bam} \
        -o {output} -of bigwig \
        -p {threads} \
        {params.others} |& \
        tee {log}
        """
