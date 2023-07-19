

rule picard_markdup:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam"),
    output:
        bam=markdup_folder.joinpath("{serie}/{sample}.markdup.bam"),
        stats=markdup_folder.joinpath("{serie}/{sample}.markdup.stats.txt"),
    log:
        log_folder.joinpath("picard/{serie}/{sample}.log"),
    threads: 2
    conda:
        "../env/picard.yml"
    shell:
        """
        set -x

        picard MarkDuplicates \
        I={input} \
        O={output.bam} \
        M={output.stats} |& \
        tee {log}
        """


rule fastqc_markdup:
    input:
        markdup_folder.joinpath("{serie}/{sample}.markdup.bam"),
    output:
        fastqc_markdup_folder.joinpath("{serie}", "{sample}.markdup_fastqc.zip"),
        fastqc_markdup_folder.joinpath("{serie}", "{sample}.markdup_fastqc.html"),
    params:
        fastqc_folder=fastqc_markdup_folder,
    threads: 2
    conda:
        # paths to singularity images cannot be PosixPath
        "../env/qc.yml"
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc_markdup/{serie}/{sample}.log"),
    shell:
        """

        set -x

        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}

        """


rule multiqc_markdup:
    input:
        get_markdup_bam,
        get_markdup_fastqc,
    output:
        multiqc_markdup_folder.joinpath("{serie}", "multiqc_report.html"),
    params:
        fastqc_folder=fastqc_markdup_folder,
        markdup_folder=markdup_folder,
        multiqc_folder=multiqc_markdup_folder,
    log:
        log_folder.joinpath("multiqc-markdup", "multiqc-{serie}.log"),
    conda:
        "../env/qc.yml"
    shell:
        """

        set -x

        multiqc --fullnames --dirs --export -f -o {params.multiqc_folder}/{wildcards.serie} {params.fastqc_folder}/{wildcards.serie} {params.markdup_folder}/{wildcards.serie} |& tee {log}

        """
