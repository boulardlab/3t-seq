rule starTE_shared_random:
    """Deduplicated starTE alignment (random mode)."""
    input:
        fastq=lambda wildcards: get_star_input(Struct(**HASH_TO_PARAMS["starTE"][wildcards.starte_hash])),
        star_index_folder=references_folder.joinpath("STAR"),
    output:
        bam=protected(starTE_folder.joinpath("_shared", "{starte_hash}", "random", "{sample}.Aligned.out.bam")),
        log=protected(starTE_folder.joinpath("_shared", "{starte_hash}", "random", "{sample}.Log.final.out")),
    threads: 8
    resources:
        runtime=lambda wildcards, input, attempt: get_star_runtime(wildcards, input.fastq, attempt),
        mem_mb=lambda wildcards, input: get_star_mem_mb(wildcards, input.fastq),
    params:
        libtype=lambda wildcards: (
            "SINGLE" if HASH_TO_PARAMS["starTE"][wildcards.starte_hash]["paired"] == False else "PAIRED"
        ),
        out_prefix=lambda wildcards: str(starTE_folder.joinpath("_shared", wildcards.starte_hash, "random", wildcards.sample)) + ".",
        outFilterMultimapNmax=lambda wildcards: get_resolved_param(wildcards, "starTE_random", "outFilterMultimapNmax", default=5000),
        winAnchorMultimapNmax=lambda wildcards: get_resolved_param(wildcards, "starTE_random", "winAnchorMultimapNmax", default=5000),
        alignTranscriptsPerWindowNmax=lambda wildcards: get_resolved_param(wildcards, "starTE_random", "alignTranscriptsPerWindowNmax", default=300),
    conda:
        "../env/alignment.yml"
    shell:
        """
         set -e 
         TMP_FOLDER=$(mktemp -u -p {resources.tmpdir})
         
         STAR \
            --outSAMtype BAM Unsorted \
            --runMode alignReads \
            --outFilterMultimapNmax {params.outFilterMultimapNmax} \
            --outSAMmultNmax 1 \
            --outFilterMismatchNmax 3 \
            --outMultimapperOrder Random \
            --winAnchorMultimapNmax {params.winAnchorMultimapNmax} \
            --alignEndsType EndToEnd \
            --alignIntronMax 1 \
            --alignMatesGapMax 350 \
            --seedSearchStartLmax 30 \
            --alignTranscriptsPerReadNmax 30000 \
            --alignWindowsPerReadNmax 30000 \
            --alignTranscriptsPerWindowNmax {params.alignTranscriptsPerWindowNmax} \
            --seedPerReadNmax 3000 \
            --seedPerWindowNmax 300 \
            --seedNoneLociPerWindow 1000 \
            --readFilesCommand zcat \
            --genomeLoad NoSharedMemory \
            --outTmpDir $TMP_FOLDER \
            --runThreadN {threads} \
            --genomeDir {input.star_index_folder} \
            --outFileNamePrefix {params.out_prefix} \
            --readFilesIn {input.fastq} \
            --limitBAMsortRAM {resources.mem_mb} \
            --outBAMcompression -1

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
         """

rule symlink_starTE_random:
    input:
        bam=lambda wildcards: get_shared_starTE_path(get_sample_hash(wildcards.serie, wildcards.sample, "starTE"), wildcards.sample, "random", ".Aligned.out.bam"),
        log=lambda wildcards: get_shared_starTE_path(get_sample_hash(wildcards.serie, wildcards.sample, "starTE"), wildcards.sample, "random", ".Log.final.out"),
    output:
        bam=starTE_folder.joinpath("{serie}/random/{sample}.Aligned.out.bam"),
        log=starTE_folder.joinpath("{serie}/random/{sample}.Log.final.out"),
    shell:
        """
        ln -sfr {input.bam} {output.bam}
        ln -sfr {input.log} {output.log}
        """


rule featureCounts_random:
    input:
        bam=lambda wildcards: expand(
            starTE_folder.joinpath("{serie}/filter/random/{sample}.TEonly.bam"),
            serie=wildcards.serie,
            sample=get_samples_names(wildcards),
        ),
        annotation=rmsk_folder.joinpath(
            "{0}.{1}".format(config["genome"].get("label", "custom"), "gtf")
        ),
    output:
        starTE_folder.joinpath("{serie}/featureCount/random.txt"),
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("featureCounts/{serie}/random.log"),
    threads: 4
    resources:
        runtime=360,
        mem_mb=lambda wildcards, input: get_featurecounts_mem_mb(wildcards, input.bam),
    params:
        strandedness=lambda wildcards: get_resolved_param(wildcards, "strandedness", default=0),
    shell:
        """
         set -e 
         featureCounts -M -F GTF -T {threads} -s {params.strandedness} -a {input.annotation} -o {output} -g repName {input.bam}
        """


rule deseq2_starTE_random:
    input:
        counts=starTE_folder.joinpath("{serie}/featureCount/random.txt"),
        sample_sheet=get_sample_sheet,
    output:
        dds=starTE_folder.joinpath("{serie}", "DESeq2", "dds_random.rds"),
        deg_table=starTE_folder.joinpath("{serie}", "DESeq2", "lfc.txt"),
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
        log_folder.joinpath("R/{serie}/deseq2-starTE-random.log"),
    script:
        "../scripts/deseq2_starTE_random_v1.R"


localrules:
    yte_starTE_random,
    datavzrd_starTE_random,


rule yte_starTE_random:
    input:
        datasets=[starTE_folder.joinpath("{serie}", "DESeq2", "lfc.txt")],
    output:
        starTE_folder.joinpath("{serie}", "datavzrd.yaml"),
    params:
        template=Path(workflow.basedir) / "datavzrd/deg-plots-template.yaml",
        plot_name="starTE-random DESeq2",
        view_specs=[str(Path(workflow.basedir) / "datavzrd/volcano-ma-plot.json")],
    conda:
        "../env/yte.yml"
    log:
        log_folder.joinpath("starTE/random/{serie}/yte.log"),
    threads: 1
    script:
        "../scripts/yte.py"


rule datavzrd_starTE_random:
    input:
        config=starTE_folder.joinpath("{serie}", "datavzrd.yaml"),
        dataset=starTE_folder.joinpath("{serie}", "DESeq2", "lfc.txt"),
    output:
        report(
            directory(starTE_folder.joinpath("{serie}", "random", "datavzrd")),
            category="starTE",
            subcategory="Random",
            labels={"serie": "{serie}", "figure": "DESeq2 analysis"},
            htmlindex="index.html",
        ),
    log:
        log_folder.joinpath("starTE/random/{serie}/datavzrd.log"),
    wrapper:
        "v2.6.0/utils/datavzrd"


rule starTE_shared_multihit:
    """Deduplicated starTE alignment (multihit mode)."""
    input:
        fastq=lambda wildcards: get_star_input(Struct(**HASH_TO_PARAMS["starTE"][wildcards.starte_hash])),
        star_index_folder=references_folder.joinpath("STAR"),
    output:
        bam=protected(starTE_folder.joinpath("_shared", "{starte_hash}", "multihit", "{sample}.Aligned.out.bam")),
        log=protected(starTE_folder.joinpath("_shared", "{starte_hash}", "multihit", "{sample}.Log.final.out")),
    threads: 8
    resources:
        runtime=lambda wildcards, input, attempt: get_star_runtime(wildcards, input.fastq, attempt),
        mem_mb=lambda wildcards, input: get_star_mem_mb(wildcards, input.fastq),
    params:
        libtype=lambda wildcards: (
            "SINGLE" if HASH_TO_PARAMS["starTE"][wildcards.starte_hash]["paired"] == False else "PAIRED"
        ),
        out_prefix=lambda wildcards: str(starTE_folder.joinpath("_shared", wildcards.starte_hash, "multihit", wildcards.sample)) + ".",
        outFilterMultimapNmax=lambda wildcards: get_resolved_param(wildcards, "starTE_multihit", "outFilterMultimapNmax", default=1),
        winAnchorMultimapNmax=lambda wildcards: get_resolved_param(wildcards, "starTE_multihit", "winAnchorMultimapNmax", default=5000),
        alignTranscriptsPerWindowNmax=lambda wildcards: get_resolved_param(wildcards, "starTE_multihit", "alignTranscriptsPerWindowNmax", default=3000),
    conda:
        "../env/alignment.yml"
    shell:
        """
         set -e 
         TMP_FOLDER=$(mktemp -u -p {resources.tmpdir})
                  
         STAR \
            --outSAMtype BAM Unsorted \
            --runMode alignReads \
            --outFilterMultimapNmax {params.outFilterMultimapNmax} \
            --outFilterMismatchNmax 3 \
            --outMultimapperOrder Random \
            --winAnchorMultimapNmax {params.winAnchorMultimapNmax} \
            --alignEndsType EndToEnd \
            --alignIntronMax 1 \
            --alignMatesGapMax 350 \
            --seedSearchStartLmax 30 \
            --alignTranscriptsPerReadNmax 30000 \
            --alignWindowsPerReadNmax 30000 \
            --alignTranscriptsPerWindowNmax {params.alignTranscriptsPerWindowNmax} \
            --seedPerReadNmax 3000 \
            --seedPerWindowNmax 300 \
            --seedNoneLociPerWindow 1000 \
            --outTmpDir $TMP_FOLDER \
            --runThreadN {threads} \
            --genomeDir {input.star_index_folder} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.out_prefix} \
            --readFilesIn {input.fastq} \
            --limitBAMsortRAM {resources.mem_mb} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression -1

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
         """

rule symlink_starTE_multihit:
    input:
        bam=lambda wildcards: get_shared_starTE_path(get_sample_hash(wildcards.serie, wildcards.sample, "starTE"), wildcards.sample, "multihit", ".Aligned.out.bam"),
        log=lambda wildcards: get_shared_starTE_path(get_sample_hash(wildcards.serie, wildcards.sample, "starTE"), wildcards.sample, "multihit", ".Log.final.out"),
    output:
        bam=starTE_folder.joinpath("{serie}/multihit/{sample}.Aligned.out.bam"),
        log=starTE_folder.joinpath("{serie}/multihit/{sample}.Log.final.out"),
    shell:
        """
        ln -sfr {input.bam} {output.bam}
        ln -sfr {input.log} {output.log}
        """


rule featureCounts_multihit:
    input:
        bam=lambda wildcards: expand(
            starTE_folder.joinpath("{serie}/filter/multihit/{sample}.TEonly.bam"),
            serie=wildcards.serie,
            sample=get_samples_names(wildcards),
        ),
        annotation=rmsk_folder.joinpath(
            "{0}.{1}".format(config["genome"].get("label", "custom"), "gtf")
        ),
    output:
        starTE_folder.joinpath("{serie}/featureCount/multihit.txt"),
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("featureCounts/{serie}/multihit.log"),
    threads: 4
    resources:
        runtime=360,
        mem_mb=lambda wildcards, input: get_featurecounts_mem_mb(wildcards, input.bam),
    params:
        strandedness=lambda wildcards: get_resolved_param(wildcards, "strandedness", default=0),
    shell:
        """
         set -e 
         featureCounts -M --fraction -F GTF -T {threads} -s {params.strandedness} -a {input.annotation} -o {output} -g repName {input.bam}
         """
