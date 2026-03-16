

rule picard_markdup:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        bam=markdup_folder.joinpath("{serie}/{sample}.markdup.bam"),
        stats=markdup_folder.joinpath("{serie}/{sample}.markdup.stats.txt"),
    log:
        log_folder.joinpath("picard/{serie}/{sample}.log"),
    threads: 2
    resources:
        runtime=360,
        mem_mb=16000,
    conda:
        "../env/picard.yml"
    shell:
        """
        set -e 

        picard MarkDuplicates \
        I={input} \
        O={output.bam} \
        M={output.stats} |& \
        tee {log}
        """


rule fastqc_markdup:
    input:
        unpack(get_markdup_bam),
    output:
        fastqc_markdup_folder.joinpath("{serie}", "{sample}.markdup_fastqc.zip"),
        report(
            fastqc_markdup_folder.joinpath("{serie}", "{sample}.markdup_fastqc.html"),
            category="FastQC",
            subcategory="Deduplicated alignments",
            labels={"serie": "{serie}", "sample": "{sample}"},
        ),
    params:
        fastqc_folder=fastqc_markdup_folder,
    threads: 4
    resources:
        runtime=20,
        mem_mb=4000,
    conda:
        # paths to singularity images cannot be PosixPath
        "../env/qc.yml"
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc_markdup/{serie}/{sample}.log"),
    shell:
        """

        set -e 

        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input.bam}

        """


