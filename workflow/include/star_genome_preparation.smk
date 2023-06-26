rule star_genome_preparation:
    input:
        genome_fasta_file=fasta_path,
        genome_annotation_file=gtf_path
    output:
        directory(references_folder.joinpath("STAR"))
    params:
        tmp_folder=tmp_folder.joinpath("STAR_genome_prep")
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    threads:
        8
    log:
        log_folder.joinpath("star/genome_preparation.log")
    shell:
        """
        set -x
        
        [[ -d {params.tmp_folder} ]] && rm -rf {params.tmp_folder}
        
        STAR --runMode genomeGenerate \
        --outTmpDir {params.tmp_folder} \
        --runThreadN {threads} \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome_fasta_file} \
        --sjdbGTFfile {input.genome_annotation_file} \
        --sjdbOverhang 100 |& \
        tee {log}
        
        [[ -d {params.tmp_folder} ]] && rm -rf {params.tmp_folder} || exit 0
        """
