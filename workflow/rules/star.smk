rule validate_genome_and_annotation:
    input:
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path,
    output:
        touch(references_folder.joinpath("genome-and-annotation-validated.done")),
    cache: True
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


rule star:
    input:
        references_folder.joinpath("genome-and-annotation-validated.done"),
        bam=get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
        genome_annotation_file=gtf_path,
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
        ),
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
    shadow:
        "minimal"
    conda:
        "../env/alignment.yml"
    log:
        star_folder.joinpath("{serie}", "{sample}.Log.final.out"),
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
         {params.others}

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


rule multiqc_star:
    input:
        unpack(get_multiqc_star_inputs),
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
