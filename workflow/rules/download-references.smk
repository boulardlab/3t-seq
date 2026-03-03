rule refgenie_init:
    output:
        refgenie_cfg=str(references_folder.joinpath("refgenie", "genome_config.yaml")),
    log:
        log_folder.joinpath("download/genome/refgenie_init.log"),
    conda:
        "../env/refgenie.yml"
    threads: 1
    resources:
        runtime=10,
        mem_mb=1000,
    shell:
        """
        export REFGENIE={output.refgenie_cfg}
        if [ ! -f "$REFGENIE" ]; then
            refgenie init -c $REFGENIE &> {log}
        fi
        """

rule refgenie_fetch_fasta:
    input:
        refgenie_cfg=str(references_folder.joinpath("refgenie", "genome_config.yaml")),
    output:
        fasta=protected(str(references_folder.joinpath(config["genome"]["label"], "fasta.fa"))),
    params:
        genome_id=config["genome"]["label"],
        selected_chrs=" ".join(config["genome"]["selected_chromosomes"]) if config["genome"]["selected_chromosomes"] else "",
    conda:
        "../env/refgenie.yml"
    log:
        log_folder.joinpath("download/genome/{}_fasta.log".format(config["genome"]["label"])),
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000,
    shell:
        """
        export REFGENIE={input.refgenie_cfg}
        refgenie pull {params.genome_id}/fasta &> {log}
        FASTA_REAL_PATH=$(refgenie seek {params.genome_id}/fasta)
        
        bash ../scripts/format_fasta.sh "$FASTA_REAL_PATH" "{output.fasta}" "{params.selected_chrs}"
        """

rule refgenie_fetch_gtf:
    input:
        refgenie_cfg=str(references_folder.joinpath("refgenie", "genome_config.yaml")),
    output:
        gtf=protected(str(references_folder.joinpath(config["genome"]["label"], "annotation.gtf"))),
    cache: True
    conda:
        "../env/refgenie.yml"
    params:
        genome_id=config["genome"]["label"],
        selected_chrs=",".join(config["genome"]["selected_chromosomes"]) if config["genome"]["selected_chromosomes"] else "",
    log:
        log_folder.joinpath("download/genome/{}_gtf.log".format(config["genome"]["label"])),
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000,
    shell:
        """
        export REFGENIE={input.refgenie_cfg}
        refgenie pull {params.genome_id}/ensembl_gtf &> {log} || refgenie pull {params.genome_id}/gencode_gtf &>> {log}
        GTF_REAL_PATH=$(refgenie seek {params.genome_id}/ensembl_gtf) || GTF_REAL_PATH=$(refgenie seek {params.genome_id}/gencode_gtf)
        
        bash ../scripts/format_gtf.sh "$GTF_REAL_PATH" "{output.gtf}" "{params.selected_chrs}"
        """

if "fasta_path" in config["genome"]:

    rule format_custom_fasta:
        input:
            fasta=config["genome"]["fasta_path"],
        output:
            fasta=protected(str(fasta_path)),
        params:
            selected_chrs=" ".join(config["genome"]["selected_chromosomes"]) if config["genome"]["selected_chromosomes"] else "",
        conda:
            "../env/bash.yml"
        log:
            log_folder.joinpath("download/genome/custom_fasta.log"),
        threads: 1
        resources:
            runtime=60,
            mem_mb=4000,
        shell:
            """
            bash ../scripts/format_fasta.sh "{input.fasta}" "{output.fasta}" "{params.selected_chrs}"
            """

    rule format_custom_gtf:
        input:
            gtf=config["genome"]["gtf_path"],
        output:
            gtf=protected(str(gtf_path)),
        conda:
            "../env/bash.yml"
        params:
            selected_chrs=",".join(config["genome"]["selected_chromosomes"]) if config["genome"]["selected_chromosomes"] else "",
        log:
            log_folder.joinpath("download/genome/custom_gtf.log"),
        threads: 1
        resources:
            runtime=60,
            mem_mb=4000,
        shell:
            """
            bash ../scripts/format_gtf.sh "{input.gtf}" "{output.gtf}" "{params.selected_chrs}"
            """




rule download_repeatmasker_annotation_file:
    output:
        protected(
            multiext(
                str(rmsk_folder.joinpath(config["genome"]["label"])), ".gtf", ".bed"
            )
        ),
    cache: True
    params:
        genome_id=config["genome"]["label"],
        selected_chromosome=config["genome"]["selected_chromosomes"],
    conda:
        "../env/pandas.yml"  # use a Python env, the script does not really use Pandas
    log:
        log_folder.joinpath("download/genome/rmsk.log"),
    threads: 1
    resources:
        runtime=20,
        mem_mb=4000,
    script:
        "../scripts/get_rmsk.py"


GTRNADB_EXTS = [
    "-filtered-tRNAs.fa",
    "-mature-tRNAs.fa",
    "-tRNAs_name_map.txt",
    "-tRNAs-confidence-set.out",
    "-tRNAs-confidence-set.ss",
    "-tRNAs-detailed.out",
    "-tRNAs-detailed.ss",
    "-tRNAs.bed",
    "-tRNAs.fa",
]

if config["genome"]["selected_chromosomes"]:
    gtrnadb_raw_dir = tRNA_annotation_dir.joinpath("raw")

    rule download_gtRNAdb:
        output:
            temp(
                multiext(
                    str(gtrnadb_raw_dir.joinpath(config["genome"]["label"])),
                    *GTRNADB_EXTS
                )
            ),
        cache: True
        params:
            url=config["genome"]["gtrnadb_url"],
            output_dir=str(gtrnadb_raw_dir),
        log:
            log_folder.joinpath("download/genome/gtrnadb.log"),
        conda:
            "../env/wget.yml"
        script:
            "../scripts/download-gtrnadb.sh"

    rule filter_gtRNAdb:
        input:
            multiext(
                str(gtrnadb_raw_dir.joinpath(config["genome"]["label"])),
                *GTRNADB_EXTS
            )
        output:
            protected(
                multiext(
                    str(tRNA_annotation_dir.joinpath(config["genome"]["label"])),
                    *GTRNADB_EXTS
                )
            )
        params:
            selected_chromosomes=config["genome"]["selected_chromosomes"]
        log:
            log_folder.joinpath("download/genome/filter_gtrnadb.log"),
        script:
            "../scripts/filter_gtrnadb.py"

else:
    rule download_gtRNAdb:
        output:
            protected(
                multiext(
                    str(tRNA_annotation_dir.joinpath(config["genome"]["label"])),
                    *GTRNADB_EXTS
                )
            ),
        cache: True
        params:
            url=config["genome"]["gtrnadb_url"],
            output_dir=str(tRNA_annotation_dir),
        log:
            log_folder.joinpath("download/genome/gtrnadb.log"),
        conda:
            "../env/wget.yml"
        script:
            "../scripts/download-gtrnadb.sh"
