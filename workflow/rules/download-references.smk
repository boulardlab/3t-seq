rule refgenie_fetch_fasta:
    output:
        fasta=protected(str(references_folder.joinpath("{genome}", "fasta.fa"))),
    params:
        refgenie_cfg=str(references_folder.joinpath("refgenie", "genome_config.yaml")),
    conda:
        "../env/refgenie.yml"
    log:
        log_folder.joinpath("download/genome/{genome}_fasta.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000,
    shell:
        """
        export REFGENIE={params.refgenie_cfg}
        if [ ! -f "$REFGENIE" ]; then
            refgenie init -c $REFGENIE
        fi
        refgenie pull {wildcards.genome}/fasta &> {log}
        FASTA_REAL_PATH=$(refgenie seek {wildcards.genome}/fasta)
        ln -s $FASTA_REAL_PATH {output.fasta}
        """

rule refgenie_fetch_gtf:
    output:
        gtf=protected(str(references_folder.joinpath("{genome}", "annotation.gtf"))),
    cache: True
    conda:
        "../env/refgenie.yml"
    params:
        refgenie_cfg=str(references_folder.joinpath("refgenie", "genome_config.yaml")),
    log:
        log_folder.joinpath("download/genome/{genome}_gtf.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000,
    shell:
        """
        export REFGENIE={params.refgenie_cfg}
        if [ ! -f "$REFGENIE" ]; then
            refgenie init -c $REFGENIE
        fi
        refgenie pull {wildcards.genome}/ensembl_gtf &> {log} || refgenie pull {wildcards.genome}/gencode_gtf &>> {log}
        GTF_REAL_PATH=$(refgenie seek {wildcards.genome}/ensembl_gtf) || GTF_REAL_PATH=$(refgenie seek {wildcards.genome}/gencode_gtf)
        ln -s $GTF_REAL_PATH {output.gtf}
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


rule download_gtRNAdb:
    output:
        protected(
            multiext(
                str(tRNA_annotation_dir.joinpath(config["genome"]["label"])),
                "-filtered-tRNAs.fa",
                "-mature-tRNAs.fa",
                "-tRNAs_name_map.txt",
                "-tRNAs-confidence-set.out",
                "-tRNAs-confidence-set.ss",
                "-tRNAs-detailed.out",
                "-tRNAs-detailed.ss",
                "-tRNAs.bed",
                "-tRNAs.fa",
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
