
localrules:
    mk_genome_tsv,


rule mk_genome_tsv:
    input:
        star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        star_folder.joinpath("{serie}", "{sample}.genome"),
    conda:
        "../env/samtools.yml"
    log:
        log_folder.joinpath("mk_genome_tsv/{serie}/{sample}.log"),
    threads: 1
    resources:
        runtime=10,
        mem_mb=2048,
    shell:
        """
        samtools view -H {input} | grep SQ | cut -f 2,3 | sed -r 's/(SN|LN)://g' | tr " " "\t" > {output}
        """


rule coverage_trna:
    input:
        bam=star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
        bai=star_folder.joinpath(
            "{serie}", "{sample}.Aligned.sortedByCoord.out.bam.bai"
        ),
        genome=star_folder.joinpath("{serie}", "{sample}.genome"),
        annotation=get_tRNA_annotation_file,
    output:
        trna_coverage_folder.joinpath("{serie}", "{sample}.bed"),
    conda:
        "../env/bedtools.yml"
    log:
        log_folder.joinpath("bedtools-trna/{serie}/{sample}.log"),
    threads: 1
    resources:
        runtime=120,
        mem_mb=16000,
    script:
        "../scripts/coverage_trna.sh"


rule build_trna_coverage_matrix:
    input:
        get_trna_coverage,
    output:
        trna_coverage_folder.joinpath("{serie}", "tRNA_matrix.txt"),
    conda:
        "../env/R.yml"
    log:
        log_folder.joinpath("bedtools-trna/build_trna_coverage_matrix-{serie}.log"),
    threads: 1
    resources:
        runtime=20,
        mem_mb=32000,
    script:
        "../scripts/build_tRNA_coverage_matrix_v1.R"
