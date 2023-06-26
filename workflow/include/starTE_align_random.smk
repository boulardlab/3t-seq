from utilities.general import giga_to_byte

rule starTE_random:
    input:
        get_star_input,
        star_index_folder = references_folder.joinpath("STAR")
    output:
          starTE_folder.joinpath("{serie}/random/{sample}.Aligned.out.bam"),
    threads: 8
    params:
        libtype = lambda wildcards: "SINGLE" if wildcards.serie in library_names_single else "PAIRED",
        alignments_folder=starTE_folder,
        tmp_folder=tmp_folder,
        mem_mb=giga_to_byte(32)
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    log:
        log_folder.joinpath("starTE/random/{serie}/{sample}.log")
    shell:
         """
         set -x
         TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
         echo 'tmp: $TMP_FOLDER'
         sleep 10
         
         if [ {params.libtype} == "SINGLE" ]; then
            INPUTARG="{input[0]}"
         else
            INPUTARG="{input[0]} {input[1]}"
         fi

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
            --readFilesIn $INPUTARG \
            --limitBAMsortRAM {params.mem_mb} \
            --outBAMcompression -1 |& \
         tee {log}

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
         """

# rule starTE_random_se:
#     input:
#         trim_reads_folder.joinpath("{se_serie}", "{sample}.fastq.gz"),
#         star_index_folder = references_folder.joinpath("STAR")
#     output:
#           starTE_folder.joinpath("{se_serie}/random/{sample}.Aligned.out.bam"),
#     resources:
#         mem_mb=giga_to_byte(32)
#     threads: 8
#     params:
#         alignments_folder=starTE_folder,
#         tmp_folder=tmp_folder
#     singularity:
#         str(container_folder.joinpath("alignment.sif"))
#     log:
#         log_folder.joinpath("starTE/random/{se_serie}/{sample}.log")
#     shell:
#          """
#          set -x
#          TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
#          echo 'tmp: $TMP_FOLDER'
#          sleep 10
#          STAR \
#             --outSAMtype BAM Unsorted \
#             --runMode alignReads \
#             --outFilterMultimapNmax 5000 \
#             --outSAMmultNmax 1 \
#             --outFilterMismatchNmax 3 \
#             --outMultimapperOrder Random \
#             --winAnchorMultimapNmax 5000 \
#             --alignEndsType EndToEnd \
#             --alignIntronMax 1 \
#             --alignMatesGapMax 350 \
#             --seedSearchStartLmax 30 \
#             --alignTranscriptsPerReadNmax 30000 \
#             --alignWindowsPerReadNmax 30000 \
#             --alignTranscriptsPerWindowNmax 300 \
#             --seedPerReadNmax 3000 \
#             --seedPerWindowNmax 300 \
#             --seedNoneLociPerWindow 1000 \
#             --readFilesCommand zcat \
#             --genomeLoad NoSharedMemory \
#             --outTmpDir $TMP_FOLDER \
#             --runThreadN {threads} \
#             --genomeDir {input.star_index_folder} \
#             --outFileNamePrefix {params.alignments_folder}/{wildcards.se_serie}/random/{wildcards.sample}. \
#             --readFilesIn {input[0]} \
#             --limitBAMsortRAM {resources.mem_mb} \
#             --outBAMcompression -1 |& \
#          tee {log}
#
#          [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
#          """
#
# rule starTE_random_pe:
#     input:
#         trim_reads_folder.joinpath("{pe_serie}", "{sample}_1.fastq.gz"),
#         trim_reads_folder.joinpath("{pe_serie}", "{sample}_2.fastq.gz"),
#         star_index_folder=references_folder.joinpath("STAR")
#     output:
#           starTE_folder.joinpath("{pe_serie}/random/{sample}.Aligned.out.bam"),
#     resources:
#         mem_mb=giga_to_byte(32)
#     threads: 8
#     params:
#         alignments_folder=starTE_folder,
#         tmp_folder=tmp_folder
#     singularity:
#         str(container_folder.joinpath("alignment.sif"))
#     log:
#         log_folder.joinpath("starTE/random/{pe_serie}/{sample}.log")
#     shell:
#          """
#          set -x
#          TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
#          echo 'tmp: $TMP_FOLDER'
#          sleep 10
#          STAR \
#             --outSAMtype BAM Unsorted \
#             --runMode alignReads \
#             --outFilterMultimapNmax 5000 \
#             --outSAMmultNmax 1 \
#             --outFilterMismatchNmax 3 \
#             --outMultimapperOrder Random \
#             --winAnchorMultimapNmax 5000 \
#             --alignEndsType EndToEnd \
#             --alignIntronMax 1 \
#             --alignMatesGapMax 350 \
#             --seedSearchStartLmax 30 \
#             --alignTranscriptsPerReadNmax 30000 \
#             --alignWindowsPerReadNmax 30000 \
#             --alignTranscriptsPerWindowNmax 300 \
#             --seedPerReadNmax 3000 \
#             --seedPerWindowNmax 300 \
#             --seedNoneLociPerWindow 1000 \
#             --readFilesCommand zcat \
#             --genomeLoad NoSharedMemory \
#             --outTmpDir $TMP_FOLDER \
#             --runThreadN {threads} \
#             --genomeDir {input.star_index_folder} \
#             --outFileNamePrefix {params.alignments_folder}/{wildcards.pe_serie}/random/{wildcards.sample}. \
#             --readFilesIn {input[0]} {input[1]} \
#             --limitBAMsortRAM {resources.mem_mb} \
#             --outBAMcompression -1 |& \
#          tee {log}
#
#          [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
#          """
