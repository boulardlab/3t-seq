from utilities.rnaseq import get_samples

rule picard_markdup:
    input:
        star_folder.joinpath("{serie}/{sample}.Aligned.sortedByCoord.out.bam")
    output:
        markdup_folder.joinpath("{serie}/{sample}.markdup.bam"),
        markdup_folder.joinpath("{serie}/{sample}.markdup.stats.txt")
    log:
        log_folder.joinpath("picard/{serie}/{sample}.log")
    threads:
        2
    singularity:
        str(container_folder.joinpath("picard.sif"))
    shell:
        # absolute path to jar file inside the container! If not using Singularity this might break.
        """
        PICARD=/usr/picard/picard.jar # path inside the container
        set -x
        java -XX:ParallelGCThreads={threads} \
        -jar $PICARD \
        MarkDuplicates I={input} \
        O={output[0]} \
        M={output[1]} |& \
        tee {log}
        """


rule fastqc_markdup:
    input:
        markdup_folder.joinpath("{serie}/{sample}.markdup.bam")
    output:
        fastqc_markdup_folder.joinpath("{serie}", "{sample}.markdup_fastqc.zip"),
        fastqc_markdup_folder.joinpath("{serie}", "{sample}.markdup_fastqc.html")
    params:
        fastqc_folder = fastqc_markdup_folder
    threads: 2
    singularity:
        # paths to singularity images cannot be PosixPath
        str(container_folder.joinpath("qc.sif"))
    log:
        log_folder.joinpath("fastqc_markdup/{serie}/{sample}.log")
    shell:
        """
        set -x
        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}
        """


def get_markdup_bam(wildcards):
    return expand(markdup_folder.joinpath("{{serie}}/{sample}.markdup.bam"),sample=get_samples(wildcards,samples))


def get_markdup_fastqc(wildcards):
    return expand(fastqc_markdup_folder.joinpath("{{serie}}","{sample}.markdup_fastqc.html"), sample = get_samples(wildcards, samples))


rule multiqc_markdup:
    input:
        get_markdup_bam,
        get_markdup_fastqc
    output:
        multiqc_markdup_folder.joinpath("{serie}","multiqc_report.html")
    params:
        fastqc_folder  = fastqc_markdup_folder,
        markdup_folder   = markdup_folder,
        multiqc_folder = multiqc_markdup_folder
    log:
        log_folder.joinpath("multiqc-markdup","multiqc-{serie}.log")
    singularity:
        # paths to singularity images cannot be PosixPath
        str(container_folder.joinpath("qc.sif"))
    shell:
        """
        set -x
        multiqc --fullnames --dirs --export -f \
        -o {params.multiqc_folder}/{wildcards.serie} \
        {params.fastqc_folder}/{wildcards.serie} \
        {params.markdup_folder}/{wildcards.serie} |& tee {log}
        """
