localrules: symlink_markdup

rule picard_markdup_shared:
    """Deduplicated Picard MarkDuplicates."""
    input:
        lambda wildcards: get_shared_star_path(
            HASH_TO_PARAMS["markdup"][wildcards.markdup_hash]["align_hash"],
            wildcards.sample,
            ".Aligned.sortedByCoord.out.bam",
        ),
    output:
        bam=protected(
            markdup_folder.joinpath(
                "_shared", "{markdup_hash}", "{sample}.markdup.bam"
            )
        ),
        stats=protected(
            markdup_folder.joinpath(
                "_shared", "{markdup_hash}", "{sample}.markdup.stats.txt"
            )
        ),
    log:
        protected(markdup_folder.joinpath("_shared", "{markdup_hash}", "{sample}.log")),
    conda:
        "../env/picard.yml"
    threads: 2
    resources:
        runtime=360,
        mem_mb=16000,
    shell:
        """
        set -e 

        picard MarkDuplicates \
        I={input} \
        O={output.bam} \
        M={output.stats} |& \
        tee {log}
        """


rule symlink_markdup:
    """Links shared markdup results back to per-serie folders."""
    input:
        bam=lambda wildcards: get_shared_markdup_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "markdup"),
            wildcards.sample,
            ".markdup.bam",
        ),
        stats=lambda wildcards: get_shared_markdup_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "markdup"),
            wildcards.sample,
            ".markdup.stats.txt",
        ),
        log=lambda wildcards: get_shared_markdup_path(
            get_sample_hash(wildcards.serie, wildcards.sample, "markdup"),
            wildcards.sample,
            ".log",
        ),
    output:
        bam=markdup_folder.joinpath("{serie}/{sample}.markdup.bam"),
        stats=markdup_folder.joinpath("{serie}/{sample}.markdup.stats.txt"),
        log=log_folder.joinpath("picard/{serie}/{sample}.log"),
    log:
        log_folder.joinpath("symlink/markdup/{serie}/{sample}.log"),
    shell:
        """
        ln -sfr {input.bam} {output.bam}
        ln -sfr {input.stats} {output.stats}
        ln -sfr {input.log} {output.log}
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
    log:
        log_folder.joinpath("fastqc_markdup/{serie}/{sample}.log"),
    conda:
        "../env/qc.yml"
    threads: 4
    resources:
        runtime=20,
        mem_mb=4000,
    params:
        fastqc_folder=fastqc_markdup_folder,
    shell:
        """

        set -e 

        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input.bam}

        """
