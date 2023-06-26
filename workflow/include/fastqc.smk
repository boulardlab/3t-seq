from utilities.rnaseq import get_samples


rule fastqc_raw_se:
    input:
        raw_reads_folder.joinpath("{serie}/{sample}.fq.gz")
    output:
        fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.zip"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.html")
    params:
        fastqc_folder=fastqc_raw_folder
    threads: 2
    singularity:
        # paths to singularity images cannot be PosixPath
        str(container_folder.joinpath("qc.sif"))
    log:
        log_folder.joinpath("fastqc/{serie}/{sample}.log")
    shell:
        """
        set -x
        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}
        """

rule fastqc_raw_pe:
    input:
        raw_reads_folder.joinpath("{serie}/{sample}_1_sequence.fq.gz"),
        raw_reads_folder.joinpath("{serie}/{sample}_2_sequence.fq.gz")
    output:
        fastqc_raw_folder.joinpath("{serie}", "{sample}_1_sequence.fq.gz_fastqc.zip"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_1_sequence.fq.gz_fastqc.html"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_2_sequence.fq.gz_fastqc.zip"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_2_sequence.fq.gz_fastqc.html")
    params:
        fastqc_folder=fastqc_raw_folder
    threads: 2
    singularity:
        # paths to singularity images cannot be PosixPath
        str(container_folder.joinpath("qc.sif"))
    log:
        log_folder.joinpath("fastqc/{serie}/{sample}.log")
    shell:
        """
        set -x
        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}
        """


def get_fastqc(wildcards):
    if wildcards.serie in library_names_single:
        ret = expand(fastqc_raw_folder.joinpath("{{serie}}/{sample}_fastqc.html"), sample=get_samples(wildcards,samples))
    else:
        ret = expand(fastqc_raw_folder.joinpath("{{serie}}/{sample}_1_sequence_fastqc.html"),sample=get_samples(wildcards,samples)) + \
              expand(fastqc_raw_folder.joinpath("{{serie}}/{sample}_2_sequence_fastqc.html"),sample=get_samples(wildcards,samples))
    return ret


rule multiqc_raw:
    input:
        get_fastqc
    output:
        multiqc_raw_folder.joinpath("{serie}", "multiqc_report.html")
    params:
        fastqc_folder = fastqc_raw_folder,
        multiqc_folder = multiqc_raw_folder
    log:
        log_folder.joinpath("multiqc-raw","multiqc-{serie}.log")
    singularity:
        # paths to singularity images cannot be PosixPath
        str(container_folder.joinpath("qc.sif"))
    shell:
        """
        set -x
        multiqc --fullnames --dirs --export -f \
        -o {params.multiqc_folder}/{wildcards.serie} \
        {params.fastqc_folder}/{wildcards.serie} |& tee {log}
        """
