
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


rule deseq2_tRNA:
    input:
        counts=trna_coverage_folder.joinpath("{serie}", "tRNA_matrix.txt"),
        sample_sheet=get_sample_sheet,
    output:
        dds=trna_coverage_folder.joinpath("{serie}", "tRNA_dds.rds"),
        deg_table=trna_coverage_folder.joinpath("{serie}", "tRNA_lfc.txt"),
    params:
        variable=lambda wildcards: get_deseq2_variable(wildcards),
        reference_level=lambda wildcards: get_deseq2_reference_level(wildcards),
    conda:
        "../env/R.yml"
    threads: 4
    resources:
        runtime=40,
        mem_mb=20000,
    log:
        log_folder.joinpath("R/{serie}/deseq2-trna.log"),
    script:
        "../scripts/deseq2_trna_v1.R"


localrules:
    yte_trna,
    datavzrd_trna,


rule yte_trna:
    input:
        template=workflow.source_path("../datavzrd/deg-plots-template.yaml"),
        datasets=[trna_coverage_folder.joinpath("{serie}", "tRNA_lfc.txt")],
    output:
        trna_coverage_folder.joinpath("{serie}", "datavzrd.yaml"),
    params:
        plot_name="tRNA expression",
        view_specs=[workflow.source_path("../datavzrd/volcano-ma-plot.json")],
    conda:
        "../env/yte.yml"
    log:
        log_folder.joinpath("bedtools-trna/yte-{serie}.log"),
    threads: 1
    script:
        "../scripts/yte.py"


rule datavzrd_trna:
    input:
        config=trna_coverage_folder.joinpath("{serie}", "datavzrd.yaml"),
        dataset=trna_coverage_folder.joinpath("{serie}", "tRNA_lfc.txt"),
    output:
        report(
            directory(trna_coverage_folder.joinpath("{serie}", "datavzrd")),
            category="tRNA expression",
            subcategory="Differential expression",
            labels={"serie": "{serie}", "figure": "DESeq2 analysis"},
            htmlindex="index.html",
        ),
    log:
        log_folder.joinpath("bedtools-trna/datavzrd-{serie}.log"),
    wrapper:
        "v2.6.0/utils/datavzrd"
