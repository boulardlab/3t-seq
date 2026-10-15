

rule fastqc_raw_se:
    input:
        fastq=get_fastq,
    output:
        zip=fastqc_raw_folder.joinpath("{serie}", "{sample}_raw_fastqc.zip"),
        html=report(
            fastqc_raw_folder.joinpath("{serie}", "{sample}_raw_fastqc.html"),
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
        src=$(realpath {input.fastq})
        ext="${{src##*.}}"
        [ "$ext" = "gz" ] && link={wildcards.sample}.fastq.gz || link={wildcards.sample}.fastq

        if [ ! -z $TMPDIR ]; then
            echo "Using TMPDIR: $TMPDIR" | tee -a {log}
            ln -sf "$src" $TMPDIR/$link
            fastqc -t {threads} -noextract -o $TMPDIR $TMPDIR/$link | tee -a {log}

            sleep $(( 5 + RANDOM % 5 ))
            mv $TMPDIR/{wildcards.sample}_fastqc.zip {output.zip}
            mv $TMPDIR/{wildcards.sample}_fastqc.html {output.html}
        else
            echo "Using default output directory: {params.fastqc_folder}" | tee -a {log}
            ln -sf "$src" {params.fastqc_folder}/$link
            fastqc -t {threads} -noextract -o {params.fastqc_folder} {params.fastqc_folder}/$link | tee -a {log}

            sleep $(( 5 + RANDOM % 5 ))
            mv {params.fastqc_folder}/{wildcards.sample}_fastqc.zip {output.zip}
            mv {params.fastqc_folder}/{wildcards.sample}_fastqc.html {output.html}
            rm -f {params.fastqc_folder}/$link
        fi
        """


rule fastqc_raw_pe:
    input:
        unpack(get_fastq_paired),
    output:
        zip1=fastqc_raw_folder.joinpath("{serie}", "{sample}_raw_1_fastqc.zip"),
        html1=report(
            fastqc_raw_folder.joinpath("{serie}", "{sample}_raw_1_fastqc.html"),
            category="FastQC",
            subcategory="Raw reads",
            labels={"serie": "{serie}", "sample": "{sample}", "mate": "1"},
        ),
        zip2=fastqc_raw_folder.joinpath("{serie}", "{sample}_raw_2_fastqc.zip"),
        html2=report(
            fastqc_raw_folder.joinpath("{serie}", "{sample}_raw_2_fastqc.html"),
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
            mv $TMPDIR/{wildcards.sample}_1_fastqc.zip {output.zip1}
            mv $TMPDIR/{wildcards.sample}_1_fastqc.html {output.html1}
            mv $TMPDIR/{wildcards.sample}_2_fastqc.zip {output.zip2}
            mv $TMPDIR/{wildcards.sample}_2_fastqc.html {output.html2}

        else
            echo "Using default output directory: {params.fastqc_folder}" | tee -a {log}
            ln -sf "$src1" {params.fastqc_folder}/$link1
            ln -sf "$src2" {params.fastqc_folder}/$link2
            fastqc -t {threads} -noextract -o {params.fastqc_folder} {params.fastqc_folder}/$link1 {params.fastqc_folder}/$link2 | tee -a {log}

            sleep $(( 5 + RANDOM % 5 ))
            mv {params.fastqc_folder}/{wildcards.sample}_1_fastqc.zip {output.zip1}
            mv {params.fastqc_folder}/{wildcards.sample}_1_fastqc.html {output.html1}
            mv {params.fastqc_folder}/{wildcards.sample}_2_fastqc.zip {output.zip2}
            mv {params.fastqc_folder}/{wildcards.sample}_2_fastqc.html {output.html2}
            rm -f {params.fastqc_folder}/$link1 {params.fastqc_folder}/$link2
        fi
        """
