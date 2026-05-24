

rule fastqc_raw_se:
    input:
        get_fastq,
    output:
        fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.zip"),
        report(
            fastqc_raw_folder.joinpath("{serie}", "{sample}_fastqc.html"),
            category="FastQC",
            subcategory="Raw reads",
            labels={"serie": "{serie}", "sample": "{sample}"},
        ),
    log:
        log_folder.joinpath("fastqc/{serie}/{sample}.log"),
    wildcard_constraints:
        serie="|".join(
            library_names_single if len(library_names_single) > 0 else ["none"]
        ),
    conda:
        "../env/qc.yml"
    threads: 8
    resources:
        runtime=10,
        mem_mb=4000,
    params:
        fastqc_folder=lambda wildcards: os.path.join(fastqc_raw_folder, wildcards.serie),
    shell:
        """
        set -e
        fastqc -t {threads} -noextract -o {params.fastqc_folder} {input}

        # FastQC names outputs from the input stem (strips compression then format extension).
        # Strip .gz/.bz2 then the remaining extension to match what FastQC produces.
        f=$(basename {input}); f="${{f%.gz}}"; f="${{f%.bz2}}"; f="${{f%.*}}"
        mv {params.fastqc_folder}/${{f}}_fastqc.zip {output[0]}
        mv {params.fastqc_folder}/${{f}}_fastqc.html {output[1]}
        """


rule fastqc_raw_pe:
    input:
        unpack(get_fastq_paired),
    output:
        fastqc_raw_folder.joinpath("{serie}", "{sample}_1_fastqc.zip"),
        report(
            fastqc_raw_folder.joinpath("{serie}", "{sample}_1_fastqc.html"),
            category="FastQC",
            subcategory="Raw reads",
            labels={"serie": "{serie}", "sample": "{sample}", "mate": "1"},
        ),
        fastqc_raw_folder.joinpath("{serie}", "{sample}_2_fastqc.zip"),
        report(
            fastqc_raw_folder.joinpath("{serie}", "{sample}_2_fastqc.html"),
            category="FastQC",
            subcategory="Raw reads",
            labels={"serie": "{serie}", "sample": "{sample}", "mate": "2"},
        ),
    log:
        log_folder.joinpath("fastqc/{serie}/{sample}_pe.log"),
    wildcard_constraints:
        serie="|".join(
            library_names_paired if len(library_names_paired) > 0 else ["none"]
        ),
    conda:
        "../env/qc.yml"
    threads: 8
    resources:
        runtime=10,
        mem_mb=4000,
    params:
        fastqc_folder=lambda wildcards: os.path.join(fastqc_raw_folder, wildcards.serie),
    shell:
        """
        set -e
        fastqc -t {threads} -noextract -o {params.fastqc_folder} {input.m1} {input.m2}

        # FastQC names outputs from the input stem (strips compression then format extension).
        # Strip .gz/.bz2 then the remaining extension to match what FastQC produces.
        m1=$(basename {input.m1}); m1="${{m1%.gz}}"; m1="${{m1%.bz2}}"; m1="${{m1%.*}}"
        m2=$(basename {input.m2}); m2="${{m2%.gz}}"; m2="${{m2%.bz2}}"; m2="${{m2%.*}}"
        mv {params.fastqc_folder}/${{m1}}_fastqc.zip {output[0]}
        mv {params.fastqc_folder}/${{m1}}_fastqc.html {output[1]}
        mv {params.fastqc_folder}/${{m2}}_fastqc.zip {output[2]}
        mv {params.fastqc_folder}/${{m2}}_fastqc.html {output[3]}
        """
