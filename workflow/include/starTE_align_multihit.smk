
rule starTE_multihit:
    input:
        get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
    threads: 8
    output:
        starTE_folder.joinpath("{serie}/multihit/{sample}.Aligned.out.bam"),
    params:
        libtype=lambda wildcards: "SINGLE"
        if wildcards.serie in library_names_single
        else "PAIRED",
        alignments_folder=starTE_folder,
        tmp_folder=tmp_folder,
        mem_mb=giga_to_byte(32),
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("starTE/{serie}/multihit/{sample}.log"),
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
            --readFilesIn $INPUTARG \
            --limitBAMsortRAM {params.mem_mb} \
            --genomeLoad NoSharedMemory \
            --outBAMcompression -1 |& \
         tee {log}

         [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
         """


# rule starTE_multihit_se:
#     input:
#         trim_reads_folder.joinpath("{se_serie}", "{se_sample}.fastq.gz"),
#         star_index_folder=references_folder.joinpath("STAR")
#     output:
#           starTE_folder.joinpath("{se_serie}/multihit/{se_sample}.Aligned.out.bam")
#     resources:
#         mem_mb=giga_to_byte(32)
#     threads: 8
#     params:
#         alignments_folder=starTE_folder,
#         tmp_folder=tmp_folder
#     conda:
#         "../env/alignment.yml"
#     log:
#         log_folder.joinpath("starTE/{se_serie}/multihit/{se_sample}.log")
#     shell:
#          """
#          set -x
#          TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
#          echo 'tmp: $TMP_FOLDER'
#          sleep 10
#          STAR \
#             --outSAMtype BAM Unsorted \
#             --runMode alignReads \
#             --outFilterMultimapNmax 1 \
#             --outFilterMismatchNmax 3 \
#             --outMultimapperOrder Random \
#             --winAnchorMultimapNmax 5000 \
#             --alignEndsType EndToEnd \
#             --alignIntronMax 1 \
#             --alignMatesGapMax 350 \
#             --seedSearchStartLmax 30 \
#             --alignTranscriptsPerReadNmax 30000 \
#             --alignWindowsPerReadNmax 30000 \
#             --alignTranscriptsPerWindowNmax 3000 \
#             --seedPerReadNmax 3000 \
#             --seedPerWindowNmax 300 \
#             --seedNoneLociPerWindow 1000 \
#             --outTmpDir $TMP_FOLDER \
#             --runThreadN {threads} \
#             --genomeDir {input.star_index_folder} \
#             --readFilesCommand zcat \
#             --outFileNamePrefix {params.alignments_folder}/{wildcards.se_serie}/multihit/{wildcards.se_sample}. \
#             --readFilesIn {input[0]} \
#             --limitBAMsortRAM {resources.mem_mb} \
#             --genomeLoad NoSharedMemory \
#             --outBAMcompression -1 |& \
#          tee {log}
#
#          [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
#          """
#
# rule starTE_multihit_pe:
#     input:
#         trim_reads_folder.joinpath("{pe_serie}", "{pe_sample}_1.fastq.gz"),
#         trim_reads_folder.joinpath("{pe_serie}", "{pe_sample}_2.fastq.gz"),
#         star_index_folder=references_folder.joinpath("STAR")
#     output:
#           starTE_folder.joinpath("{pe_serie}/multihit/{pe_sample}.Aligned.out.bam")
#     resources:
#         mem_mb=giga_to_byte(32)
#     threads: 8
#     params:
#         alignments_folder=starTE_folder,
#         tmp_folder=tmp_folder
#     conda:
#         "../env/alignment.yml"
#     log:
#         log_folder.joinpath("starTE/{pe_serie}/multihit/{pe_sample}.log")
#     shell:
#          """
#          set -x
#          TMP_FOLDER=$(mktemp -u -p {params.tmp_folder})
#          echo 'tmp: $TMP_FOLDER'
#          sleep 10
#          STAR \
#             --outSAMtype BAM Unsorted \
#             --runMode alignReads \
#             --outFilterMultimapNmax 1 \
#             --outFilterMismatchNmax 3 \
#             --outMultimapperOrder Random \
#             --winAnchorMultimapNmax 5000 \
#             --alignEndsType EndToEnd \
#             --alignIntronMax 1 \
#             --alignMatesGapMax 350 \
#             --seedSearchStartLmax 30 \
#             --alignTranscriptsPerReadNmax 30000 \
#             --alignWindowsPerReadNmax 30000 \
#             --alignTranscriptsPerWindowNmax 3000 \
#             --seedPerReadNmax 3000 \
#             --seedPerWindowNmax 300 \
#             --seedNoneLociPerWindow 1000 \
#             --outTmpDir $TMP_FOLDER \
#             --runThreadN {threads} \
#             --genomeDir {input.star_index_folder} \
#             --readFilesCommand zcat \
#             --outFileNamePrefix {params.alignments_folder}/{wildcards.pe_serie}/multihit/{wildcards.pe_sample}. \
#             --readFilesIn {input[0]} {input[1]} \
#             --limitBAMsortRAM {resources.mem_mb} \
#             --genomeLoad NoSharedMemory \
#             --outBAMcompression -1 |& \
#          tee {log}
#
#          [[ -d $TMP_FOLDER ]] && rm -r $TMP_FOLDER || exit 0
#          """
