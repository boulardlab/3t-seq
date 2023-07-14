rule starTE_random:
    input:
        get_star_input,
        star_index_folder=references_folder.joinpath("STAR"),
    output:
        starTE_folder.joinpath("{serie}/random/{sample}.Aligned.out.bam"),
    threads: 8
    params:
        libtype=lambda wildcards: "SINGLE"
        if wildcards.serie in library_names_single
        else "PAIRED",
        alignments_folder=starTE_folder,
        tmp_folder=tmp_folder,
        mem_mb=giga_to_byte(32),
    conda:
        "../../env/alignment.yml"
    log:
        log_folder.joinpath("starTE/random/{serie}/{sample}.log"),
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

rule featureCounts_random:
    input:
        bam=lambda wildcards: expand(
            starTE_folder.joinpath("{serie}/filter/random/{sample}.TEonly.bam"),
            serie=wildcards.serie,
            sample=samples["single"][wildcards.serie]
            if wildcards.serie in samples["single"]
            else samples["paired"][wildcards.serie],
        ),
        annotation=rmsk_path,
    output:
        starTE_folder.joinpath("{serie}/featureCount/random.txt"),
    conda:
        "../../env/alignment.yml"
    log:
        log_folder.joinpath("featureCounts/{serie}/random.log"),
    threads: 4
    shell:
        """
         set -x
         featureCounts -M -F GTF -T {threads} -s 0 -a {input.annotation} -o {output} {input.bam}
         """

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
        "../../env/alignment.yml"
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

rule featureCounts_multihit:
    input:
        bam=lambda wildcards: expand(
            starTE_folder.joinpath("{serie}/filter/multihit/{sample}.TEonly.bam"),
            serie=wildcards.serie,
            sample=samples["single"][wildcards.serie]
            if wildcards.serie in samples["single"]
            else samples["paired"][wildcards.serie],
        ),
        annotation=rmsk_path,
    output:
        starTE_folder.joinpath("{serie}/featureCount/multihit.txt"),
    conda:
        "../../env/alignment.yml"
    log:
        log_folder.joinpath("featureCounts/{serie}/multihit.log"),
    threads: 4
    shell:
        """
         set -x
         featureCounts -M --fraction -F GTF -T {threads} -s 0 -a {input.annotation} -o {output} {input.bam}
         """
