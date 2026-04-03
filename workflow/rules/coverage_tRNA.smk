
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
        samtools view -H {input} | grep '@SQ' | cut -f 2,3 | sed -r 's/(SN|LN)://g' | tr " " "\t" > {output}
        """


rule coverage_trna:
    input:
        bam=star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out.bam"),
        bai=star_folder.joinpath(
            "{serie}", "{sample}.Aligned.sortedByCoord.out.bam.bai"
        ),
        genome=star_folder.joinpath("{serie}", "{sample}.genome"),
        annotation=tRNA_annotation_dir.joinpath(
            "{0}-tRNAs.bed".format(config.get("genome", {}).get("label", "custom"))
        ),
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
        unpack(get_trna_coverage),
    output:
        trna_coverage_folder.joinpath("{serie}", "tRNA_matrix_standard.txt"),
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
        counts=trna_coverage_folder.joinpath("{serie}", "tRNA_matrix_{method}.txt"),
        sample_sheet=get_sample_sheet,
    output:
        dds=trna_coverage_folder.joinpath("{serie}", "tRNA_dds_{method}.rds"),
        deg_table=trna_coverage_folder.joinpath("{serie}", "tRNA_lfc_{method}.txt"),
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
        log_folder.joinpath("R/{serie}/deseq2-trna-{method}.log"),
    script:
        "../scripts/deseq2_trna_v1.R"


localrules:
    yte_trna,
    datavzrd_trna,


rule yte_trna:
    input:
        datasets=[trna_coverage_folder.joinpath("{serie}", "tRNA_lfc_{method}.txt")],
    output:
        trna_coverage_folder.joinpath("{serie}", "datavzrd_{method}.yaml"),
    params:
        template=Path(workflow.basedir) / "datavzrd/deg-plots-template.yaml",
        plot_name="tRNA expression",
        view_specs=[str(Path(workflow.basedir) / "datavzrd/volcano-ma-plot.json")],
    conda:
        "../env/yte.yml"
    log:
        log_folder.joinpath("bedtools-trna/yte-{serie}-{method}.log"),
    threads: 1
    script:
        "../scripts/yte.py"


rule datavzrd_trna:
    input:
        config=trna_coverage_folder.joinpath("{serie}", "datavzrd_{method}.yaml"),
        dataset=trna_coverage_folder.joinpath("{serie}", "tRNA_lfc_{method}.txt"),
    output:
        report(
            directory(trna_coverage_folder.joinpath("{serie}", "datavzrd_{method}")),
            category="tRNA expression",
            subcategory="Differential expression",
            labels={"serie": "{serie}", "figure": "DESeq2 analysis"},
            htmlindex="index.html",
        ),
    log:
        log_folder.joinpath("bedtools-trna/datavzrd-{serie}-{method}.log"),
    wrapper:
        "v9.4.1/utils/datavzrd"


# --- mim-tRNA-seq Rules --- #


rule build_mimseq_index:
    input:
        genome_fasta=fasta_path,
        annotation=tRNA_annotation_dir.joinpath(
            "{0}-tRNAs.bed".format(config.get("genome", {}).get("label", "custom"))
        ),
    output:
        idx_dir=directory(tRNA_annotation_dir.joinpath("mimseq_index")),
    conda:
        "../env/mim-trna-seq.yml"
    params:
        species="custom",
        cluster_id=lambda wildcards: config.get("tRNA_quantification", {})
        .get("mimseq_params", {})
        .get("cluster_identity", 0.95),
    threads: 4
    resources:
        runtime=120,
        mem_mb=16000,
    log:
        log_folder.joinpath("mimseq/build_index.log"),
    shell:
        """
        set -e
        mimseq --make-index --species {params.species} --genome {input.genome_fasta} --annotation {input.annotation} --out-dir {output.idx_dir} --threads {threads} --cluster-id {params.cluster_id} &> {log}
        """


rule run_mimseq:
    input:
        fastq=get_trimmed_fastq,  # Depends on whether single or paired
        idx_dir=tRNA_annotation_dir.joinpath("mimseq_index"),
    output:
        count_file=trna_coverage_folder.joinpath(
            "mimseq_{serie}", "{sample}", "{sample}_cluster_counts.txt"
        ),
    conda:
        "../env/mim-trna-seq.yml"
    params:
        species="custom",
        max_mismatches=lambda wildcards: config.get("tRNA_quantification", {})
        .get("mimseq_params", {})
        .get("max_mismatches", 0.1),
        min_cov=lambda wildcards: config.get("tRNA_quantification", {})
        .get("mimseq_params", {})
        .get("min_cov", 10),
        out_dir=lambda wildcards, output: Path(output.count_file).parent,
    threads: 8
    resources:
        runtime=240,
        mem_mb=32000,
    log:
        log_folder.joinpath("mimseq/{serie}/{sample}.log"),
    shell:
        """
        set -e
        # mimseq takes --cluster-id string to prefix outputs. 
        # If paired-end, mimseq doesn't officially support paired natively without mapping as single-end or using specialized pipeline variants, 
        # but latest mimseq 1.3+ allows feeding standard FASTQs. We pass all FASTQs via {input.fastq} 
        mimseq --species {params.species} -i {input.idx_dir} -n {wildcards.sample} --max-mismatches {params.max_mismatches} --min-cov {params.min_cov} --threads {threads} --out-dir {params.out_dir} {input.fastq} &> {log}
        """


rule format_mimseq_output:
    input:
        mimseq_counts=get_mimseq_counts,
        sample_sheet=get_sample_sheet,
    output:
        trna_coverage_folder.joinpath("{serie}", "tRNA_matrix_mimseq.txt"),
    conda:
        "../env/R.yml"
    threads: 1
    resources:
        runtime=20,
        mem_mb=8000,
    log:
        log_folder.joinpath("R/{serie}/format_mimseq.log"),
    script:
        "../scripts/format_mimseq_matrix.R"
