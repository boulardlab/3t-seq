rule download_genome_fasta_file:
    output:
        protected(fasta_path),
    params:
        url=config["genome"]["fasta_url"],
    cache: True
    conda:
        "../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/fasta.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000,
    script:
        "../scripts/download-fasta.sh"


rule download_genome_annotation_file:
    output:
        protected(gtf_path),
    cache: True
    params:
        url=config["genome"]["gtf_url"],
    conda:
        "../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/gtf.log"),
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000,
    script:
        "../scripts/download-gtf.sh"


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
        output_dir=tRNA_annotation_dir,
    log:
        log_folder.joinpath("download/genome/gtrnadb.log"),
    conda:
        "../env/wget.yml"
    script:
        "../scripts/download-gtrnadb.sh"


# rule download_gaf_file:
#     output:
#         gaf_path,
#     params:
#         url=config["genome"]["gaf_url"],
#     conda:
#         "../env/wget.yml"
#     log:
#         log_folder.joinpath("download/genome/gaf.log"),
#     threads: 4
#     script:
#         "../scripts/download-gaf.sh"
