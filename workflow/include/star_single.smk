from utilities.general import giga_to_byte
from utilities.rnaseq import get_samples

rule verify_star:
    input: lambda wildcards: expand(star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"), serie=wildcards.serie, sample=samples['single'][wildcards.serie] if wildcards.serie in samples['single'] else samples['paired'][wildcards.serie])
    output: touch(star_folder.joinpath("{serie}.done"))

rule index_bam:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam")
    output:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam.bai")
    threads:
        2
    log:
        log_folder.joinpath("index-bam/{serie}/{sample}.log")
    singularity:
        str(container_folder.joinpath("samtools.sif"))
    shell:
        """
        samtools index {input}
        """


rule star:
    input:
        get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
        genome_annotation_file=gtf_path
    output:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
        star_folder.joinpath("{serie}/{sample}.Aligned.toTranscriptome.out.bam"),
        star_folder.joinpath("{serie}/{sample}.ReadsPerGene.out.tab"),
        star_folder.joinpath("{serie}/{sample}.Log.final.out")
    # resources:
    threads: 8
    params:
        libtype = lambda wildcards: "SINGLE" if wildcards.serie in library_names_single else "PAIRED",
        alignments_folder = star_folder,
        tmp_folder = tmp_folder,
        others = lambda wildcards: get_params(wildcards, "star"),
        mem_mb=giga_to_byte(32)
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    log:
        log_folder.joinpath("star/{serie}/{sample}.log")
    shell:
         """
         set -x
         TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
         
         if [ {params.libtype} == "SINGLE" ]; then
            INPUTARG="{input[0]}"
         else
            INPUTARG="{input[0]} {input[1]}"
         fi
                 
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
         --readFilesIn $INPUTARG \
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


# rule star_se:
#     input:
#         trim_reads_folder.joinpath("{se_serie}", "{sample}.fastq.gz"),
#         star_index_folder=references_folder.joinpath("STAR"),
#         genome_annotation_file=gtf_path
#     output:
#         star_folder.joinpath("{se_serie}/{sample}.Aligned.sortedByCoord.out.bam"),
#         star_folder.joinpath("{se_serie}/{sample}.Aligned.toTranscriptome.out.bam"),
#         star_folder.joinpath("{se_serie}/{sample}.ReadsPerGene.out.tab"),
#         star_folder.joinpath("{se_serie}/{sample}.Log.final.out")
#     resources:
#         mem_mb=giga_to_byte(32)
#     threads: 8
#     params:
#         alignments_folder=star_folder,
#         tmp_folder=tmp_folder,
#         others=lambda wildcards: config["sequencing_libraries"][wildcards.se_serie]["star"]
#     singularity:
#         str(container_folder.joinpath("alignment.sif"))
#     log:
#         log_folder.joinpath("star/{se_serie}/{sample}.log")
#     shell:
#          """
#          set -x
#          TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
#          STAR --quantMode TranscriptomeSAM GeneCounts \
#          --outTmpDir $TMP_FOLDER \
#          --outSAMtype BAM SortedByCoordinate \
#          --sjdbGTFfile {input.genome_annotation_file} \
#          --runThreadN {threads} \
#          --chimOutType WithinBAM \
#          --twopassMode Basic \
#          --genomeDir {input.star_index_folder} \
#          --readFilesCommand zcat \
#          --outFileNamePrefix {params.alignments_folder}/{wildcards.se_serie}/{wildcards.sample}. \
#          --readFilesIn {input[0]} \
#          --limitBAMsortRAM {resources.mem_mb} \
#          --genomeLoad NoSharedMemory \
#          --outSAMunmapped Within \
#          --outReadsUnmapped FastX \
#          --outBAMsortingThreadN {threads} \
#          --bamRemoveDuplicatesType UniqueIdentical \
#          --quantTranscriptomeBAMcompression -1 \
#          --outBAMcompression -1 \
#          --outWigType wiggle \
#          {params.others} |& \
#          tee {log}
#
#          [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
#          """
#
#
# rule star_pe:
#     input:
#         trim_reads_folder.joinpath("{pe_serie}", "{sample}_1.fastq.gz"),
#         trim_reads_folder.joinpath("{pe_serie}", "{sample}_2.fastq.gz"),
#         star_index_folder=references_folder.joinpath("STAR"),
#         genome_annotation_file=gtf_path
#     output:
#         star_folder.joinpath("{pe_serie}/{sample}.Aligned.sortedByCoord.out.bam"),
#         star_folder.joinpath("{pe_serie}/{sample}.Aligned.toTranscriptome.out.bam"),
#         star_folder.joinpath("{pe_serie}/{sample}.ReadsPerGene.out.tab"),
#         star_folder.joinpath("{pe_serie}/{sample}.Log.final.out")
#     resources:
#         mem_mb=giga_to_byte(32)
#     threads: 8
#     params:
#         alignments_folder=star_folder,
#         tmp_folder=tmp_folder,
#         others = lambda wildcards: config["sequencing_libraries"][wildcards.pe_serie]["star"]
#     singularity:
#         str(container_folder.joinpath("alignment.sif"))
#     log:
#         log_folder.joinpath("star/{pe_serie}/{sample}.log")
#     shell:
#          """
#          set -x
#          TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
#          STAR --quantMode TranscriptomeSAM GeneCounts \
#          --outTmpDir $TMP_FOLDER \
#          --outSAMtype BAM SortedByCoordinate \
#          --sjdbGTFfile {input.genome_annotation_file} \
#          --runThreadN {threads} \
#          --chimOutType WithinBAM \
#          --twopassMode Basic \
#          --genomeDir {input.star_index_folder} \
#          --readFilesCommand zcat \
#          --outFileNamePrefix {params.alignments_folder}/{wildcards.pe_serie}/{wildcards.sample}. \
#          --readFilesIn {input[0]} {input[1]} \
#          --limitBAMsortRAM {resources.mem_mb} \
#          --genomeLoad NoSharedMemory \
#          --outSAMunmapped Within \
#          --outReadsUnmapped FastX \
#          --outBAMsortingThreadN {threads} \
#          --bamRemoveDuplicatesType UniqueIdentical \
#          --quantTranscriptomeBAMcompression -1 \
#          --outBAMcompression -1 \
#          --outWigType wiggle \
#          {params.others} |& \
#          tee {log}
#
#          [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
#          """


rule fastqc_star:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam")
    output:
        fastqc_star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out_fastqc.zip"),
        fastqc_star_folder.joinpath("{serie}", "{sample}.Aligned.sortedByCoord.out_fastqc.html")
    params:
        fastqc_folder = fastqc_star_folder
    threads: 2
    singularity:
        # paths to singularity images cannot be PosixPath
        str(container_folder.joinpath("qc.sif"))
    log:
        log_folder.joinpath("fastqc_star/{serie}/{sample}.log")
    shell:
        """
        set -x
        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}
        """


def get_star_stats(wildcards):
    return expand(star_folder.joinpath("{{serie}}/{sample}.Log.final.out"), sample = get_samples(wildcards, samples))


def get_star_fastqc(wildcards):
    return expand(fastqc_star_folder.joinpath("{{serie}}", "{sample}.Aligned.sortedByCoord.out_fastqc.html"), sample = get_samples(wildcards, samples))


rule multiqc_star:
    input:
        get_star_stats,
        get_star_fastqc
    output:
        multiqc_star_folder.joinpath("{serie}","multiqc_report.html")
    params:
        fastqc_folder = fastqc_star_folder,
        star_folder = star_folder,
        multiqc_folder = multiqc_star_folder
    log:
        log_folder.joinpath("multiqc-star","multiqc-{serie}.log")
    singularity:
        # paths to singularity images cannot be PosixPath
        str(container_folder.joinpath("qc.sif"))
    shell:
        """
        set -x
        multiqc --fullnames --dirs --export -f \
        -o {params.multiqc_folder}/{wildcards.serie} \
        {params.fastqc_folder}/{wildcards.serie} \
        {params.star_folder}/{wildcards.serie} |& tee {log}
        """
