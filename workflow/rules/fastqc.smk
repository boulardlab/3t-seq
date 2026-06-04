

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
        runtime=get_fastqc_runtime,
        mem_mb=4000,
    params:
        fastqc_folder=lambda wildcards: os.path.join(fastqc_raw_folder, wildcards.serie),
    shell:
        """
        set -e
        # Symlink input under the sample name so FastQC embeds the correct sample
        # name inside the zip regardless of whether the original path is absolute
        # or relative and regardless of the original filename.
        src=$(realpath {input})
        ext="${{src##*.}}"
        [ "$ext" = "gz" ] && link={wildcards.sample}.fastq.gz || link={wildcards.sample}.fastq

        if [ ! -z $TMPDIR ]; then
            echo "Using TMPDIR: $TMPDIR" | tee -a {log}
            ln -sf "$src" $TMPDIR/$link
            fastqc -t {threads} -noextract -o $TMPDIR $TMPDIR/$link | tee -a {log}

            sleep $(( 5 + RANDOM % 5 ))
            mv $TMPDIR/{wildcards.sample}_fastqc.zip {output[0]}
            mv $TMPDIR/{wildcards.sample}_fastqc.html {output[1]}
        else
            echo "Using default output directory: {params.fastqc_folder}" | tee -a {log}
            ln -sf "$src" {params.fastqc_folder}/$link
            fastqc -t {threads} -noextract -o {params.fastqc_folder} {params.fastqc_folder}/$link | tee -a {log}

            sleep $(( 5 + RANDOM % 5 ))
            mv {params.fastqc_folder}/{wildcards.sample}_fastqc.zip {output[0]}
            mv {params.fastqc_folder}/{wildcards.sample}_fastqc.html {output[1]}
            rm -f {params.fastqc_folder}/$link
        fi
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
        runtime=get_fastqc_runtime,
        mem_mb=4000,
    params:
        fastqc_folder=lambda wildcards: os.path.join(fastqc_raw_folder, wildcards.serie),
    shell:
        """
        set -e

        # Symlink inputs under the sample name so FastQC embeds the correct sample
        # name inside the zip regardless of whether the original path is absolute
        # or relative and regardless of the original filename.
        src1=$(realpath {input.m1})
        src2=$(realpath {input.m2})
        ext1="${{src1##*.}}"
        ext2="${{src2##*.}}"
        [ "$ext1" = "gz" ] && link1={wildcards.sample}_1.fastq.gz || link1={wildcards.sample}_1.fastq
        [ "$ext2" = "gz" ] && link2={wildcards.sample}_2.fastq.gz || link2={wildcards.sample}_2.fastq

        if [ ! -z $TMPDIR ]; then
            echo "Using TMPDIR: $TMPDIR" | tee -a {log}
            ln -sf "$src1" $TMPDIR/$link1
            ln -sf "$src2" $TMPDIR/$link2
            fastqc -t {threads} -noextract -o $TMPDIR $TMPDIR/$link1 $TMPDIR/$link2 | tee -a {log}

            sleep $(( 5 + RANDOM % 5 ))
            mv $TMPDIR/{wildcards.sample}_1_fastqc.zip {output[0]}
            mv $TMPDIR/{wildcards.sample}_1_fastqc.html {output[1]}
            mv $TMPDIR/{wildcards.sample}_2_fastqc.zip {output[2]}
            mv $TMPDIR/{wildcards.sample}_2_fastqc.html {output[3]}

        else
            echo "Using default output directory: {params.fastqc_folder}" | tee -a {log}
            ln -sf "$src1" {params.fastqc_folder}/$link1
            ln -sf "$src2" {params.fastqc_folder}/$link2
            fastqc -t {threads} -noextract -o {params.fastqc_folder} {params.fastqc_folder}/$link1 {params.fastqc_folder}/$link2 | tee -a {log}

            sleep $(( 5 + RANDOM % 5 ))
            mv {params.fastqc_folder}/{wildcards.sample}_1_fastqc.zip {output[0]}
            mv {params.fastqc_folder}/{wildcards.sample}_1_fastqc.html {output[1]}
            mv {params.fastqc_folder}/{wildcards.sample}_2_fastqc.zip {output[2]}
            mv {params.fastqc_folder}/{wildcards.sample}_2_fastqc.html {output[3]}
            rm -f {params.fastqc_folder}/$link1 {params.fastqc_folder}/$link2
        fi
        """
