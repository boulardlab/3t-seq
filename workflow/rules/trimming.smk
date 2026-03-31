ruleorder: trimmomatic_shared_pe > trimmomatic_shared_se
ruleorder: symlink_trim_pe > symlink_trim_se

def get_trimmomatic_adaptive_input(wildcards):
    """Returns the adaptive params JSON if adaptive trimming is enabled."""
    if wildcards.trim_hash not in HASH_TO_PARAMS["trim"]:
        return []
    params = HASH_TO_PARAMS["trim"][wildcards.trim_hash]
    t_params = params.get("params", {})
    if isinstance(t_params, dict) and t_params.get("adaptive"):
        # We use a path that includes the sample name to avoid collisions if multiple samples share a hash
        # but have slightly different adaptive needs (though they shouldn't if the hash is correct).
        # Actually, samples with the same hash should have the same adaptive params if they are the same data.
        # But FastQC is sample-specific. So we use {sample} in the path.
        return trim_reads_folder.joinpath("_shared", wildcards.trim_hash, f"{wildcards.sample}.adaptive_params.json")
    return []

def get_trimmomatic_command(wildcards):
    """Dynamic derivation of the trimmomatic command string."""
    params = HASH_TO_PARAMS["trim"][wildcards.trim_hash]
    t_params = params.get("params", {})
    
    # 1. Adaptive mode
    if isinstance(t_params, dict) and t_params.get("adaptive"):
        adaptive_json = trim_reads_folder.joinpath("_shared", wildcards.trim_hash, f"{wildcards.sample}.adaptive_params.json")
        with open(adaptive_json, 'r') as f:
            data = json.load(f)
            return data["trimmomatic_params"]
    
    # 2. Fixed string mode (explicit)
    if isinstance(t_params, str) and t_params != "default":
        return t_params
        
    # 3. Default fallback
    is_paired = params.get("paired", False)
    suffix = "PE" if is_paired else "SE"
    return f"ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-{suffix}.fa:2:30:10 SLIDINGWINDOW:20:22 MAXINFO:4:20 LEADING:3 TRAILING:3 MINLEN:36"

rule derive_trim_params:
    """Derives adaptive trimming parameters from FastQC results."""
    input:
        fastqc=lambda wildcards: get_fastqc_raw(Struct(**HASH_TO_PARAMS["trim"][wildcards.trim_hash]))["fastqc"]
    output:
        json=trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}.adaptive_params.json")
    params:
        is_paired_flag=lambda wildcards: "--is-paired" if HASH_TO_PARAMS["trim"][wildcards.trim_hash].get("paired", False) else "",
        extra_params=lambda wildcards: HASH_TO_PARAMS["trim"][wildcards.trim_hash].get("params", {}).get("extra_params", "") if isinstance(HASH_TO_PARAMS["trim"][wildcards.trim_hash].get("params"), dict) else ""
    conda:
        "../env/qc.yml"
    shell:
        """
        python workflow/scripts/derive_trim_params.py \
            --fastqc {input.fastqc} \
            {params.is_paired_flag} \
            --extra-params "{params.extra_params}" \
            --output {output.json}
        """

rule trimmomatic_shared_pe:
    """Deduplicated trimmomatic for paired-end reads."""
    input:
        reads=unpack(lambda wildcards: get_fastq_paired(
            Struct(**HASH_TO_PARAMS["trim"][wildcards.trim_hash])
        )),
        adaptive=get_trimmomatic_adaptive_input
    output:
        paired1=protected(trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}_1.fastq.gz")),
        paired2=protected(trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}_2.fastq.gz")),
        unpaired1=protected(trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}_1.unpaired.fastq.gz")),
        unpaired2=protected(trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}_2.unpaired.fastq.gz")),
        summary=protected(trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}.summary.txt")),
    params:
        trimmomatic_cmd=lambda wildcards: get_trimmomatic_command(wildcards)
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: 240 * attempt,
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    log:
        trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}.stats.txt"),
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
        {params.trimmomatic_cmd} |& tee {log}
        """

rule symlink_trim_pe:
    input:
        paired1=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample, "_1"),
        paired2=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample, "_2"),
        unpaired1=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample, "_1.unpaired"),
        unpaired2=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample, "_2.unpaired"),
        summary=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample, "").with_suffix(".summary.txt"),
        stats=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample, "").with_suffix(".stats.txt"),
    output:
        paired1=trim_reads_folder.joinpath("{serie}", "{sample}_1.fastq.gz"),
        paired2=trim_reads_folder.joinpath("{serie}", "{sample}_2.fastq.gz"),
        unpaired1=trim_reads_folder.joinpath("{serie}", "{sample}_1.unpaired.fastq.gz"),
        unpaired2=trim_reads_folder.joinpath("{serie}", "{sample}_2.unpaired.fastq.gz"),
        summary=trim_reads_folder.joinpath("{serie}", "{sample}.summary.txt"),
        stats=trim_reads_folder.joinpath("{serie}", "{sample}.stats.txt"),
    shell:
        """
        ln -sfr {input.paired1} {output.paired1}
        ln -sfr {input.paired2} {output.paired2}
        ln -sfr {input.unpaired1} {output.unpaired1}
        ln -sfr {input.unpaired2} {output.unpaired2}
        ln -sfr {input.summary} {output.summary}
        ln -sfr {input.stats} {output.stats}
        """


rule trimmomatic_shared_se:
    """Deduplicated trimmomatic for single-end reads."""
    input:
        reads=lambda wildcards: get_fastq(
            Struct(**HASH_TO_PARAMS["trim"][wildcards.trim_hash])
        ),
        adaptive=get_trimmomatic_adaptive_input
    output:
        fastq=protected(trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}.fastq.gz")),
        summary=protected(trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}.summary.txt")),
    params:
        trimmomatic_cmd=lambda wildcards: get_trimmomatic_command(wildcards)
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: 180 * attempt,
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    log:
        trim_reads_folder.joinpath("_shared", "{trim_hash}", "{sample}.stats.txt"),
    conda:
        "../env/trimmomatic.yml"
    shell:
        """
        trimmomatic SE \
        -threads {threads} -trimlog {log} \
        -summary {output.summary} \
        {input.reads} \
        {output.fastq} \
        {params.trimmomatic_cmd} |& tee {log}
        """

rule symlink_trim_se:
    input:
        fastq=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample),
        summary=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample).with_suffix(".summary.txt"),
        stats=lambda wildcards: get_shared_trim_path(get_sample_hash(wildcards.serie, wildcards.sample, "trim"), wildcards.sample).with_suffix(".stats.txt"),
    output:
        fastq=trim_reads_folder.joinpath("{serie}", "{sample}.fastq.gz"),
        summary=trim_reads_folder.joinpath("{serie}", "{sample}.summary.txt"),
        stats=trim_reads_folder.joinpath("{serie}", "{sample}.stats.txt"),
    shell:
        """
        ln -sfr {input.fastq} {output.fastq}
        ln -sfr {input.summary} {output.summary}
        ln -sfr {input.stats} {output.stats}
        """


rule fastqc_trim:
    wildcard_constraints:
        serie="|".join(library_names_single),
    input:
        get_trimmed_fastq,
    output:
        fastqc_trim_folder.joinpath("{serie}", "{sample}_fastqc.zip"),
        report(
            fastqc_trim_folder.joinpath("{serie}", "{sample}_fastqc.html"),
            category="FastQC",
            subcategory="Trimmed reads",
            labels={"serie": "{serie}", "sample": "{sample}"},
        ),
    params:
        fastqc_folder=fastqc_trim_folder,
    threads: 4
    resources:
        runtime=90,
        mem_mb=4000,
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc_trim/{serie}/{sample}.log"),
    shell:
        """
        set -e 
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
        report(
            fastqc_trim_folder.joinpath("{serie}", "{sample}_1_fastqc.html"),
            category="FastQC",
            subcategory="Trimmed reads",
            labels={"serie": "{serie}", "sample": "{sample}", "mate": "1"},
        ),
        fastqc_trim_folder.joinpath("{serie}", "{sample}_2_fastqc.zip"),
        report(
            fastqc_trim_folder.joinpath("{serie}", "{sample}_2_fastqc.html"),
            category="FastQC",
            subcategory="Trimmed reads",
            labels={"serie": "{serie}", "sample": "{sample}", "mate": "2"},
        ),
    params:
        fastqc_folder=fastqc_trim_folder,
    threads: 4
    resources:
        runtime=90,
        mem_mb=4000,
    conda:
        "../env/qc.yml"
    log:
        log_folder.joinpath("fastqc_trim/{serie}/{sample}.log"),
    shell:
        """
        set -e 
        fastqc -t {threads} -noextract -o {params.fastqc_folder}/{wildcards.serie} {input}
        """


