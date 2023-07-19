ruleorder: trimmomatic_pe > trimmomatic_se


rule trimmomatic_pe:
    wildcard_constraints:
        serie="|".join(
            library_names_paired if len(library_names_paired) > 0 else ["none"]
        ),
    input:
        unpack(get_fastq_paired),
    output:
        paired1=trim_reads_folder.joinpath("{serie}", "{sample}_1.fastq.gz"),
        paired2=trim_reads_folder.joinpath("{serie}", "{sample}_2.fastq.gz"),
        unpaired1=trim_reads_folder.joinpath("{serie}", "{sample}_1.unpaired.fastq.gz"),
        unpaired2=trim_reads_folder.joinpath("{serie}", "{sample}_2.unpaired.fastq.gz"),
        summary=trim_reads_folder.joinpath("{serie}", "{sample}.summary.txt"),
        stats=trim_reads_folder.joinpath("{serie}", "{sample}.stats.txt"),
    params:
        lambda wildcards: get_params(wildcards, "trimmomatic"),
    threads: 2
    log:
        log_folder.joinpath("trimmomatic_pe", "{serie}", "{sample}.log"),
    conda:
        "../env/trimmomatic.yml"
    shell:
        """
        trimmomatic PE \
        -threads {threads} -trimlog {log} \
        -summary {output.summary} \
        {input.m1} {input.m2} \
        {output.paired1} {output.unpaired1} \
        {output.paired2} {output.unpaired2} \
        {params} |& tee {output.stats}
        """


rule trimmomatic_se:
    wildcard_constraints:
        serie="|".join(
            library_names_single if len(library_names_single) > 0 else ["none"]
        ),
    input:
        get_fastq,
    output:
        fastq=trim_reads_folder.joinpath("{serie}", "{sample}.fastq.gz"),
        summary=trim_reads_folder.joinpath("{serie}", "{sample}.summary.txt"),
        stats=trim_reads_folder.joinpath("{serie}", "{sample}.stats.txt"),
    params:
        lambda wildcards: get_params(wildcards, "trimmomatic"),
    threads: 2
    log:
        log_folder.joinpath("trimmomatic_se-{serie}-{sample}.log"),
    conda:
        "../env/trimmomatic.yml"
    shell:
        """
        trimmomatic SE \
        -threads {threads} -trimlog {log} \
        -summary {output.summary} \
        {input} \
        {output.fastq} \
        {params} |& tee {output.stats}
        """


rule fastqc_trim:
    wildcard_constraints:
        serie="|".join(library_names_single),
    input:
        get_trimmed_fastq,
    output:
        fastqc_trim_folder.joinpath("{serie}", "{sample}_fastqc.zip"),
        fastqc_trim_folder.joinpath("{serie}", "{sample}_fastqc.html"),
    params:
        fastqc_folder=fastqc_trim_folder,
    threads: 2
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc_trim/{serie}/{sample}.log"),
    shell:
        """
        set -x
        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}
        """


rule fastqc_trim_pe:
    wildcard_constraints:
        serie="|".join(library_names_paired),
    input:
        [
            trim_reads_folder.joinpath("{serie}/{sample}_1.fastq.gz"),
            trim_reads_folder.joinpath("{serie}/{sample}_2.fastq.gz"),
        ],
    output:
        fastqc_trim_folder.joinpath("{serie}", "{sample}_1_fastqc.zip"),
        fastqc_trim_folder.joinpath("{serie}", "{sample}_1_fastqc.html"),
        fastqc_trim_folder.joinpath("{serie}", "{sample}_2_fastqc.zip"),
        fastqc_trim_folder.joinpath("{serie}", "{sample}_2_fastqc.html"),
    params:
        fastqc_folder=fastqc_trim_folder,
    threads: 2
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc_trim/{serie}/{sample}.log"),
    shell:
        """
        set -x
        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}
        """


rule multiqc_trim:
    input:
        get_trimmomatic_stats,
        get_trimmed_fastqc,
    output:
        multiqc_trim_folder.joinpath("{serie}", "multiqc_report.html"),
    params:
        fastqc_folder=fastqc_trim_folder,
        reads_folder=trim_reads_folder,
        multiqc_folder=multiqc_trim_folder,
    log:
        log_folder.joinpath("multiqc-trim", "multiqc-{serie}.log"),
    conda:
        "../env/qc.yml"
    shell:
        """
        set -x
        multiqc --fullnames --dirs --export -f \
        -o {params.multiqc_folder}/{wildcards.serie} \
        {params.reads_folder}/{wildcards.serie} \
        {params.fastqc_folder}/{wildcards.serie} |& tee {log}
        """