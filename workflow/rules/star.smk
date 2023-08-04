rule validate_genome_and_annotation:
    input:
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path,
    output:
        touch(references_folder.joinpath("genome-and-annotation-validated.done"))
    conda:
        "../env/bash.yml"
    threads: 1
    resources:
        runtime=20,
        mem_mb=1024
    log:
        log_folder.joinpath("star/genome_preparation.log"),
    shell:
        """
        set -e

        head -n1000 {input.genome_fasta_file} | grep -q '>chr'
        FASTA_HAS_CHR=$?

        head -n1000 {input.genome_annotation_file} | grep -v '#' | head -n1 | grep -q '^chr'
        GTF_HAS_CHR=$?

        if [ $FASTA_HAS_CHR -eq 0 ] & [ $GTF_HAS_CHR -eq 1 ]; then
            sed -i -r 's/^([^#]+)$/chr\1/g' {input.genome_annotation_file}
        fi

        if [ $FASTA_HAS_CHR -eq 1 ] & [ $GTF_HAS_CHR -eq 0 ]; then
            sed -i -r 's/^>(.+)$/>chr\1/g' {input.genome_fasta_file}
        fi

        """

rule star_genome_preparation:
    input:
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path,
    output:
        directory(references_folder.joinpath("STAR")),
    params:
        tmp_folder=tmp_folder.joinpath("STAR_genome_prep"),
    conda:
        "../env/alignment.yml"
    threads: 16
    resources:
        runtime=120,
        mem_mb=32000,
    log:
        log_folder.joinpath("star/genome_preparation.log"),
    shell:
        """
        set -e 
        
        [[ -d {params.tmp_folder} ]] && rm -rf {params.tmp_folder}
        
        STAR --runMode genomeGenerate \
        --outTmpDir {params.tmp_folder} \
        --runThreadN {threads} \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.genome_annotation_file} \
        --sjdbOverhang 100 |& \
        tee {log}
        
        [[ -d {params.tmp_folder} ]] && rm -rf {params.tmp_folder} || exit 0
        """


rule star:
    input:
        bam=get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
        genome_annotation_file=gtf_path,
    output:
        temp(star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam")),
        star_folder.joinpath("{serie}/{sample}.Aligned.toTranscriptome.out.bam"),
        star_folder.joinpath("{serie}/{sample}.ReadsPerGene.out.tab"),
        star_folder.joinpath("{serie}/{sample}.Log.final.out"),
    threads: 8
    resources:
        runtime=360,
        mem_mb=32000,
    params:
        libtype=lambda wildcards: "SINGLE"
        if wildcards.serie in library_names_single
        else "PAIRED",
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
        fastqc_star_folder.joinpath(
            "{serie}", "{sample}.Aligned.sortedByCoord.out_fastqc.html"
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
            sample=samples["single"][wildcards.serie]
            if wildcards.serie in samples["single"]
            else samples["paired"][wildcards.serie],
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
        multiqc_star_folder.joinpath("{serie}", "multiqc_report.html"),
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
