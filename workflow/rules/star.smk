rule star_genome_preparation:
    input:
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path,
    output:
        directory(references_folder.joinpath("STAR")),
    cache: True
    conda:
        "../env/alignment.yml"
    threads: 8
    resources:
        runtime=120,
        mem_mb=256000,
    log:
        log_folder.joinpath("star/genome_preparation.log"),
    shell:
        """
        set -e 
        
        STAR --runMode genomeGenerate \
        --outTmpDir $(mktemp -d -u) \
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
        fastq=lambda wildcards: get_star_input(Struct(**HASH_TO_PARAMS["alignment"][wildcards.align_hash])),
        star_index_folder=references_folder.joinpath("STAR"),
        genome_annotation_file=gtf_path,
    output:
        protected(multiext(
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
        )),
    threads: 8
    resources:
        runtime=lambda wildcards, input, attempt: get_star_runtime(wildcards, input.fastq, attempt),
        mem_mb=lambda wildcards, input: get_star_mem_mb(wildcards, input.fastq),
    params:
        out_prefix=lambda wildcards: str(star_folder.joinpath("_shared", wildcards.align_hash, wildcards.sample)) + ".",
        tmp_folder=tmp_folder,
        params_others=lambda wildcards: get_params(
            Struct(**HASH_TO_PARAMS["alignment"][wildcards.align_hash]), 
            "star", 
            default="--seedSearchStartLmax 30 --outFilterMismatchNoverReadLmax 0.04 --winAnchorMultimapNmax 40"
        ),
    conda:
        "../env/alignment.yml"
    shell:
        """
         set -e 
         TMP_FOLDER=$(mktemp -u -p {resources.tmpdir})

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
        lambda wildcards: expand(
            str(star_folder.joinpath("_shared", get_sample_hash(wildcards.serie, wildcards.sample, "alignment"), wildcards.sample)) + "{ext}",
            ext=[
                ".Aligned.sortedByCoord.out.bam",
                ".Aligned.toTranscriptome.out.bam",
                ".ReadsPerGene.out.tab",
                ".SJ.out.tab",
                ".Signal.Unique.str1.out.wig",
                ".Signal.Unique.str2.out.wig",
                ".Signal.UniqueMultiple.str1.out.wig",
                ".Signal.UniqueMultiple.str2.out.wig",
                ".Log.final.out",
            ]
        )
    output:
        multiext(
            str(star_folder.joinpath("{serie}", "{sample}")),
            ".Aligned.sortedByCoord.out.bam",
            ".Aligned.toTranscriptome.out.bam",
            ".ReadsPerGene.out.tab",
            ".SJ.out.tab",
            ".Signal.Unique.str1.out.wig",
            ".Signal.Unique.str2.out.wig",
            ".Signal.UniqueMultiple.str1.out.wig",
            ".Signal.UniqueMultiple.str2.out.wig",
            ".Log.final.out",
        ),
    shell:
        """
        ln -sft {input[0]} {output[0]}
        ln -sft {input[1]} {output[1]}
        ln -sft {input[2]} {output[2]}
        ln -sft {input[3]} {output[3]}
        ln -sft {input[4]} {output[4]}
        ln -sft {input[5]} {output[5]}
        ln -sft {input[6]} {output[6]}
        ln -sft {input[7]} {output[7]}
        ln -sft {input[8]} {output[8]}
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
    params:
        fastqc_folder=fastqc_star_folder,
    threads: 2
    resources:
        runtime=20,
        mem_mb=4000,
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc_star/{serie}/{sample}.log"),
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
    threads: 4
    resources:
        runtime=30,
        mem_mb=8000,
    log:
        log_folder.joinpath("index-bam/{serie}/{sample}.log"),
    conda:
        "../env/samtools.yml"
    shell:
        """

        samtools index -@{threads} {input}

        """


