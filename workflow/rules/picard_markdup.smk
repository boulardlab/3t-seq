

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


rule multiqc_markdup:
    input:
        unpack(get_markdup_fastqc),
    output:
        report(
            multiqc_markdup_folder.joinpath("{serie}", "multiqc_report.html"),
            category="MultiQC",
            subcategory="Deduplicated alignments",
            labels={"serie": "{serie}"},
        ),
    params:
        fastqc_folder=fastqc_markdup_folder,
        markdup_folder=markdup_folder,
        multiqc_folder=multiqc_markdup_folder,
    log:
        log_folder.joinpath("multiqc-markdup", "multiqc-{serie}.log"),
    threads: 1
    resources:
        runtime=20,
        mem_mb=2048,
    conda:
        "../env/qc.yml"
    shell:
        """

        set -e 

        multiqc --fullnames --dirs --export -f -o {params.multiqc_folder}/{wildcards.serie} {params.fastqc_folder}/{wildcards.serie} {params.markdup_folder}/{wildcards.serie} |& tee {log}

        """
