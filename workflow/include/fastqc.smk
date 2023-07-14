

rule fastqc_raw:
    input:
        get_fastq
    output:
        fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.zip"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.html"),
    params:
        fastqc_folder=lambda wildcards: os.path.join(fastqc_raw_folder, wildcards.serie),
    threads: 2
    conda:
        "../../env/qc.yml"
    log:
        log_folder.joinpath("fastqc/{serie}/{sample}.log"),
    shell:
        """
        set -x
        fastqc -t {threads} -noextract -o {params.fastqc_folder} {input}
        """

use rule fastqc_raw as fastqc_raw_pe with:
    input:
        unpack(get_fastq_paired)
    output:
        fastqc_raw_folder.joinpath("{serie}", "{sample}_1_fastqc.zip"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_1_fastqc.html"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_2_fastqc.zip"),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_2_fastqc.html"),
    params:
        fastqc_folder=lambda wildcards: os.path.join(fastqc_raw_folder, wildcards.serie),
    log:
        log_folder.joinpath("fastqc/{serie}/{sample}.log"),


# rule fastqc_raw_pe:
#     input:
#         unpack(get_fastq_paired)
#     output:
#         fastqc_raw_folder.joinpath("{serie}", "{sample}_1_sequence.fq.gz_fastqc.zip"),
#         fastqc_raw_folder.joinpath("{serie}", "{sample}_1_sequence.fq.gz_fastqc.html"),
#         fastqc_raw_folder.joinpath("{serie}", "{sample}_2_sequence.fq.gz_fastqc.zip"),
#         fastqc_raw_folder.joinpath("{serie}", "{sample}_2_sequence.fq.gz_fastqc.html"),
#     params:
#         fastqc_folder=fastqc_raw_folder,
#     threads: 2
#     conda:
#         "../../env/qc.yml"
#     log:
#         log_folder.joinpath("fastqc/{serie}/{sample}.log"),
#     shell:
#         """
#         set -x
#         fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input.m1} {input.m2}
#         """


def get_fastqc(wildcards):
    s = get_samples(wildcards, samples)
    if wildcards.serie in library_names_single:
        ret = expand(
            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}_fastqc.html"),
            sample=s,
        )
    else:
        for m in supported_suffixes:
            for ext in supported_extensions:
                m1 = expand(os.path.join(raw_reads_folder, wildcards.serie, "{sample}" + f"{m[0]}.{ext}"), sample = s)
                m2 = expand(os.path.join(raw_reads_folder, wildcards.serie, "{sample}" + f"{m[1]}.{ext}"), sample = s)
                
                if all([os.path.exists(p) for p in m1 + m2]):
                    ret = expand(
                        fastqc_raw_folder.joinpath(wildcards.serie, "{sample}" + f"{m[0]}_fastqc.html"), sample=s
                        ) + expand(
                            fastqc_raw_folder.joinpath(wildcards.serie, "{sample}" + f"{m[1]}_fastqc.html"), sample=s
                        )
                                    
    return ret


rule multiqc_raw:
    input:
        get_fastqc,
    output:
        multiqc_raw_folder.joinpath("{serie}", "multiqc_report.html"),
    params:
        fastqc_folder=fastqc_raw_folder,
        multiqc_folder=multiqc_raw_folder,
    log:
        log_folder.joinpath("multiqc-raw", "multiqc-{serie}.log"),
    conda:
        "../../env/qc.yml"
    shell:
        """
        set -x
        multiqc --fullnames --dirs --export -f \
        -o {params.multiqc_folder}/{wildcards.serie} \
        {params.fastqc_folder}/{wildcards.serie} |& tee {log}
        """
