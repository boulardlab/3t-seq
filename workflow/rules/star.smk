rule star_genome_preparation:
    input:
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path,
    output:
        directory(references_folder.joinpath("STAR")),
    log:
        log_folder.joinpath("star/genome_preparation.log"),
    cache: True
    conda:
        "../env/alignment.yml"
    threads: 8
    resources:
        runtime=120,
        mem_mb=256000,
    shell:
        """
        set -e 

        TMP_FOLDER=$(mktemp -d -u -p {resources.tmpdir})
        
        STAR --runMode genomeGenerate \
        --outTmpDir $TMP_FOLDER \
        --runThreadN {threads} \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.genome_annotation_file} \
        --sjdbOverhang 100         
        
        if [ -f {output}/Log.out ]; then
          cp {output}/Log.out {log}
        elif [ -f Log.out ]; then
          cp Log.out {log}
        fi
        """


rule star_shared:
    """Deduplicated STAR alignment."""
    input:
        fastq=lambda wildcards: get_star_input(
            Struct(**HASH_TO_PARAMS["alignment"][wildcards.align_hash])
        ),
        star_index_folder=references_folder.joinpath("STAR"),
        genome_annotation_file=gtf_path,
    output:
        protected(
            multiext(
                str(star_folder.joinpath("_shared", "{align_hash}", "{sample}")),
                ".Aligned.sortedByCoord.out.bam",
                ".Aligned.toTranscriptome.out.bam",
                ".ReadsPerGene.out.tab",
                ".SJ.out.tab",
                ".Signal.Unique.str1.out.wig",
                ".Signal.Unique.str2.out.wig",
                ".Signal.UniqueMultiple.str1.out.wig",
                ".Signal.UniqueMultiple.str2.out.wig",
                ".Log.final.out",
            )
        ),
    log:
        log_folder.joinpath("star/{align_hash}/{sample}.log"),
    conda:
        "../env/alignment.yml"
    threads: 8
    resources:
        runtime=lambda wildcards, input, attempt: get_star_runtime(
            wildcards, input.fastq, attempt
        ),
        mem_mb=lambda wildcards, input: get_star_mem_mb(wildcards, input.fastq),
    params:
        out_prefix=lambda wildcards: str(
            star_folder.joinpath("_shared", wildcards.align_hash, wildcards.sample)
        )
        + ".",
        tmp_folder=tmp_folder,
        params_others=lambda wildcards: get_params(
            Struct(**HASH_TO_PARAMS["alignment"][wildcards.align_hash]),
            "star",
            default="--seedSearchStartLmax 30 --outFilterMismatchNoverReadLmax 0.04 --winAnchorMultimapNmax 40",
        ),
    shell:
        """
         set -e 
         TMP_FOLDER=$(mktemp -d -u -p {resources.tmpdir})

         STAR --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir $TMP_FOLDER \
         --outSAMtype BAM SortedByCoordinate \
         --sjdbGTFfile {input.genome_annotation_file} \
         --runThreadN {threads} \
         --chimOutType WithinBAM \
         --twopassMode Basic \
         --genomeDir {input.star_index_folder} \
         --readFilesCommand zcat \
         --outFileNamePrefix {params.out_prefix} \
         --readFilesIn {input.fastq} \
         --limitBAMsortRAM $(({resources.mem_mb} * 1024 * 1024)) \
         --genomeLoad NoSharedMemory \
         --outSAMunmapped Within \
         --outReadsUnmapped FastX \
         --outBAMsortingThreadN {threads} \
         --bamRemoveDuplicatesType UniqueIdentical \
         --quantTranscriptomeBAMcompression -1 \
         --outBAMcompression -1 --outWigType wiggle \
         {params.params_others}

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
         """


rule symlink_star:
    """Links shared STAR results back to per-serie folders."""
    input:
        bam=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".Aligned.sortedByCoord.out.bam",
        ),
        transcriptome_bam=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".Aligned.toTranscriptome.out.bam",
        ),
        gene_counts=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".ReadsPerGene.out.tab",
        ),
        sj=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".SJ.out.tab",
        ),
        wig1=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".Signal.Unique.str1.out.wig",
        ),
        wig2=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".Signal.Unique.str2.out.wig",
        ),
        wig1_mult=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".Signal.UniqueMultiple.str1.out.wig",
        ),
        wig2_mult=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".Signal.UniqueMultiple.str2.out.wig",
        ),
        log_final=lambda wildcards: get_shared_star_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "alignment"),
            wildcards.sample,
            ".Log.final.out",
        ),
    output:
        bam=star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
        transcriptome_bam=star_folder.joinpath(
            "{serie}/{sample}.Aligned.toTranscriptome.out.bam"
        ),
        gene_counts=star_folder.joinpath("{serie}/{sample}.ReadsPerGene.out.tab"),
        sj=star_folder.joinpath("{serie}/{sample}.SJ.out.tab"),
        wig1=star_folder.joinpath("{serie}/{sample}.Signal.Unique.str1.out.wig"),
        wig2=star_folder.joinpath("{serie}/{sample}.Signal.Unique.str2.out.wig"),
        wig1_mult=star_folder.joinpath(
            "{serie}/{sample}.Signal.UniqueMultiple.str1.out.wig"
        ),
        wig2_mult=star_folder.joinpath(
            "{serie}/{sample}.Signal.UniqueMultiple.str2.out.wig"
        ),
        log_final=star_folder.joinpath("{serie}/{sample}.Log.final.out"),
    log:
        log_folder.joinpath("symlink/star/{serie}/{sample}.log"),
    shell:
        """
        ln -sfr {input.bam} {output.bam}
        ln -sfr {input.transcriptome_bam} {output.transcriptome_bam}
        ln -sfr {input.gene_counts} {output.gene_counts}
        ln -sfr {input.sj} {output.sj}
        ln -sfr {input.wig1} {output.wig1}
        ln -sfr {input.wig2} {output.wig2}
        ln -sfr {input.wig1_mult} {output.wig1_mult}
        ln -sfr {input.wig2_mult} {output.wig2_mult}
        ln -sfr {input.log_final} {output.log_final}
        """


rule fastqc_star:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        fastqc_star_folder.joinpath(
            "{serie}", "{sample}.Aligned.sortedByCoord.out_fastqc.zip"
        ),
        report(
            fastqc_star_folder.joinpath(
                "{serie}", "{sample}.Aligned.sortedByCoord.out_fastqc.html"
            ),
            category="FastQC",
            subcategory="Aligned reads",
            labels={"serie": "{serie}", "sample": "{sample}"},
        ),
    log:
        log_folder.joinpath("fastqc_star/{serie}/{sample}.log"),
    conda:
        "../env/qc.yml"
    threads: 2
    resources:
        runtime=20,
        mem_mb=4000,
    params:
        fastqc_folder=fastqc_star_folder,
    shell:
        """

        set -e 

        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}

        """


rule index_bam:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam.bai"),
    log:
        log_folder.joinpath("index-bam/{serie}/{sample}.log"),
    conda:
        "../env/samtools.yml"
    threads: 4
    resources:
        runtime=30,
        mem_mb=8000,
    shell:
        """

        samtools index -@{threads} {input}

        """
