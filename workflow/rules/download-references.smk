rule download_genome_fasta_file:
    output:
        protected(fasta_path),
    params:
        url=config["genome"]["fasta_url"],
    conda:
        "../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/fasta.log"),
    threads: 1
    resources: 
        runtime=60,
        mem_mb=4000    
    script:
        "../scripts/download-fasta.sh"


rule download_genome_annotation_file:
    output:
        protected(gtf_path),
    params:
        url=config["genome"]["gtf_url"],
        tmp=config["globals"]["tmp_folder"],
    conda:
        "../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/gtf.log"),
    threads: 1
    resources: 
        runtime=60,
        mem_mb=4000
    script:
        "../scripts/download-gtf.sh"


rule download_repeatmasker_annotation_file:
    output:
        protected(rmsk_path),
    params:
        url=config["genome"]["rmsk_url"],
    conda:
        "../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/rmsk.log"),
    threads: 1
    resources: 
        runtime=20,
        mem_mb=4000
    script:
        "../scripts/download-rmsk.sh"


checkpoint download_gtRNAdb:
    output:
        protected(directory(tRNA_annotation_dir)),
    params:
        url=config["genome"]["gtrnadb_url"],
    log:
        log_folder.joinpath("download/genome/gtrnadb.log"),
    conda:
        "../env/wget.yml"
    script:
        "../scripts/download-gtrnadb.sh"


rule download_gaf_file:
    output:
        gaf_path,
    params:
        url=config["genome"]["gaf_url"],
    conda:
        "../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/gaf.log"),
    threads: 4
    script:
        "../scripts/download-gaf.sh"