rule coverage:
    input:
        bam=star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
        bai=star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam.bai"),
    output:
        star_folder.joinpath("{serie}", "{sample}.bw"),
    conda:
        "../env/deeptools.yml"
    params:
        others=lambda wildcards: get_params(wildcards, "bamCoverage"),
    threads: 2
    resources: 
        runtime=120,
        mem_mb=16000
    log:
        log_folder.joinpath("bamCoverage_se-{serie}-{sample}.log"),
    shell:
        """

        bamCoverage -b {input.bam} \
        -o {output} -of bigwig \
        -p {threads} \
        --effectiveGenomeSize 2652783500 \
        {params.others} |& \
        tee {log}
        """
