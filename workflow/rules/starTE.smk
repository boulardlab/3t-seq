rule starTE_random:
    input:
        bam=get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
    output:
        starTE_folder.joinpath("{serie}/random/{sample}.Aligned.out.bam"),
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: 1440 * attempt,
        mem_mb=32000,
    params:
        libtype=lambda wildcards: (
            "SINGLE" if wildcards.serie in library_names_single else "PAIRED"
        ),
        alignments_folder=starTE_folder,
        tmp_folder=tmp_folder,
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("starTE/random/{serie}/{sample}.log"),
    shell:
        """
         set -e 
         TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
         sleep 10
         
         STAR \
            --outSAMtype BAM Unsorted \
            --runMode alignReads \
            --outFilterMultimapNmax 5000 \
            --outSAMmultNmax 1 \
            --outFilterMismatchNmax 3 \
            --outMultimapperOrder Random \
            --winAnchorMultimapNmax 5000 \
            --alignEndsType EndToEnd \
            --alignIntronMax 1 \
            --alignMatesGapMax 350 \
            --seedSearchStartLmax 30 \
            --alignTranscriptsPerReadNmax 30000 \
            --alignWindowsPerReadNmax 30000 \
            --alignTranscriptsPerWindowNmax 300 \
            --seedPerReadNmax 3000 \
            --seedPerWindowNmax 300 \
            --seedNoneLociPerWindow 1000 \
            --readFilesCommand zcat \
            --genomeLoad NoSharedMemory \
            --outTmpDir $TMP_FOLDER \
            --runThreadN {threads} \
            --genomeDir {input.star_index_folder} \
            --outFileNamePrefix {params.alignments_folder}/{wildcards.serie}/random/{wildcards.sample}. \
            --readFilesIn {input.bam} \
            --limitBAMsortRAM {resources.mem_mb} \
            --outBAMcompression -1 |& \
         tee {log}

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
         """


rule featureCounts_random:
    input:
        bam=lambda wildcards: expand(
            starTE_folder.joinpath("{serie}/filter/random/{sample}.TEonly.bam"),
            serie=wildcards.serie,
            sample=(
                samples["single"][wildcards.serie]
                if wildcards.serie in samples["single"]
                else samples["paired"][wildcards.serie]
            ),
        ),
        annotation=rmsk_path,
    output:
        starTE_folder.joinpath("{serie}/featureCount/random.txt"),
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("featureCounts/{serie}/random.log"),
    threads: 4
    resources:
        runtime=360,
        mem_mb=16000,
    shell:
        """
         set -e 
         featureCounts -M -F GTF -T {threads} -g repName -s 0 -a {input.annotation} -o {output} {input.bam}
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
        template=workflow.source_path("../datavzrd/deg-plots-template.yaml"),
        datasets=[starTE_folder.joinpath("{serie}", "DESeq2", "lfc.txt")],
    output:
        starTE_folder.joinpath("{serie}", "datavzrd.yaml"),
    params:
        plot_name="starTE-random DESeq2",
        view_specs=[workflow.source_path("../datavzrd/volcano-ma-plot.json")],
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


rule starTE_multihit:
    input:
        bam=get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: 1440 * attempt,
        mem_mb=32000,
    output:
        starTE_folder.joinpath("{serie}/multihit/{sample}.Aligned.out.bam"),
    params:
        libtype=lambda wildcards: (
            "SINGLE" if wildcards.serie in library_names_single else "PAIRED"
        ),
        alignments_folder=starTE_folder,
        tmp_folder=tmp_folder,
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("starTE/{serie}/multihit/{sample}.log"),
    shell:
        """
         set -e 
         TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
         sleep 10
                  
         STAR \
            --outSAMtype BAM Unsorted \
            --runMode alignReads \
            --outFilterMultimapNmax 1 \
            --outFilterMismatchNmax 3 \
            --outMultimapperOrder Random \
            --winAnchorMultimapNmax 5000 \
            --alignEndsType EndToEnd \
            --alignIntronMax 1 \
            --alignMatesGapMax 350 \
            --seedSearchStartLmax 30 \
            --alignTranscriptsPerReadNmax 30000 \
            --alignWindowsPerReadNmax 30000 \
            --alignTranscriptsPerWindowNmax 3000 \
            --seedPerReadNmax 3000 \
            --seedPerWindowNmax 300 \
            --seedNoneLociPerWindow 1000 \
            --outTmpDir $TMP_FOLDER \
            --runThreadN {threads} \
            --genomeDir {input.star_index_folder} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.alignments_folder}/{wildcards.serie}/multihit/{wildcards.sample}. \
            --readFilesIn {input.bam} \
            --limitBAMsortRAM {resources.mem_mb} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression -1 |& \
         tee {log}

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
         """


rule featureCounts_multihit:
    input:
        bam=lambda wildcards: expand(
            starTE_folder.joinpath("{serie}/filter/multihit/{sample}.TEonly.bam"),
            serie=wildcards.serie,
            sample=(
                samples["single"][wildcards.serie]
                if wildcards.serie in samples["single"]
                else samples["paired"][wildcards.serie]
            ),
        ),
        annotation=rmsk_path,
    output:
        starTE_folder.joinpath("{serie}/featureCount/multihit.txt"),
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("featureCounts/{serie}/multihit.log"),
    threads: 4
    resources:
        runtime=360,
        mem_mb=16000,
    shell:
        """
         set -e 
         featureCounts -M --fraction -F GTF -T {threads} -g repName -s 0 -a {input.annotation} -o {output} {input.bam}
         """
