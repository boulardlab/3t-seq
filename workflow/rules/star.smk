rule validate_genome_and_annotation:
    input:
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path,
    output:
        touch(references_folder.joinpath("genome-and-annotation-validated.done")),
    conda:
        "../env/bash.yml"
    threads: 1
    resources:
        runtime=20,
        mem_mb=1024,
    log:
        log_folder.joinpath("star/validata_genome_and_annotation.log"),
    script:
        "../scripts/validate_genome_and_annotation.sh"


rule star_genome_preparation:
    input:
        references_folder.joinpath("genome-and-annotation-validated.done"),
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path,
    output:
        directory(references_folder.joinpath("STAR")),
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
        --sjdbOverhang 100 |& \
        tee {log}
        """


rule star:
    input:
        bam=get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
        genome_annotation_file=gtf_path,
    output:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
        star_folder.joinpath("{serie}/{sample}.Aligned.toTranscriptome.out.bam"),
        star_folder.joinpath("{serie}/{sample}.ReadsPerGene.out.tab"),
        star_folder.joinpath("{serie}/{sample}.SJ.out.tab"),
        star_folder.joinpath("{serie}/{sample}.Signal.Unique.str1.out.wig"),
        star_folder.joinpath("{serie}/{sample}.Signal.Unique.str2.out.wig"),
        star_folder.joinpath("{serie}/{sample}.Signal.UniqueMultiple.str1.out.wig"),
        star_folder.joinpath("{serie}/{sample}.Signal.UniqueMultiple.str2.out.wig"),
        star_folder.joinpath("{serie}/{sample}.Log.final.out"),
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: 1440 * attempt,
        mem_mb=32000,
    params:
        libtype=lambda wildcards: (
            "SINGLE" if wildcards.serie in library_names_single else "PAIRED"
        ),
        alignments_folder=star_folder,
        tmp_folder=tmp_folder,
        others=lambda wildcards: get_params(wildcards, "star"),
        mem_mb=giga_to_byte(32),
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("star/{serie}/{sample}.log"),
    shell:
        """
         set -e 
         TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})

         STAR --quantMode TranscriptomeSAM GeneCounts \
         --outTmpDir $TMP_FOLDER \
         --outSAMtype BAM SortedByCoordinate \
         --sjdbGTFfile {input.genome_annotation_file} \
         --runThreadN {threads} \
         --chimOutType WithinBAM \
         --twopassMode Basic \
         --genomeDir {input.star_index_folder} \
         --readFilesCommand zcat \
         --outFileNamePrefix {params.alignments_folder}/{wildcards.serie}/{wildcards.sample}. \
         --readFilesIn {input.bam} \
         --limitBAMsortRAM {params.mem_mb} \
         --genomeLoad NoSharedMemory \
         --outSAMunmapped Within \
         --outReadsUnmapped FastX \
         --outBAMsortingThreadN {threads} \
         --bamRemoveDuplicatesType UniqueIdentical \
         --quantTranscriptomeBAMcompression -1 \
         --outBAMcompression -1 \
         --outWigType wiggle \
         {params.others} |& \
         tee {log}

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
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


rule verify_star:
    input:
        lambda wildcards: expand(
            star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
            serie=wildcards.serie,
            sample=(
                samples["single"][wildcards.serie]
                if wildcards.serie in samples["single"]
                else samples["paired"][wildcards.serie]
            ),
        ),
    output:
        touch(star_folder.joinpath("{serie}.done")),


rule index_bam:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam.bai"),
    threads: 1
    resources:
        runtime=30,
        mem_mb=8000,
    log:
        log_folder.joinpath("index-bam/{serie}/{sample}.log"),
    conda:
        "../env/samtools.yml"
    shell:
        """

        samtools index {input}

        """


rule multiqc_star:
    input:
        get_star_stats,
        get_star_fastqc,
    output:
        report(
            multiqc_star_folder.joinpath("{serie}", "multiqc_report.html"),
            category="MultiQC",
            subcategory="Alignment",
            labels={"serie": "{serie}"},
        ),
    params:
        fastqc_folder=fastqc_star_folder,
        star_folder=star_folder,
        multiqc_folder=multiqc_star_folder,
    log:
        log_folder.joinpath("multiqc-star", "multiqc-{serie}.log"),
    threads: 1
    resources:
        runtime=10,
        mem_mb=2048,
    conda:
        # paths to singularity images cannot be PosixPath
        "../env/qc.yml"
    shell:
        """
        set -e 
        multiqc --fullnames --dirs --export -f -o {params.multiqc_folder}/{wildcards.serie} {params.fastqc_folder}/{wildcards.serie} {params.star_folder}/{wildcards.serie} |& tee {log}
        """
