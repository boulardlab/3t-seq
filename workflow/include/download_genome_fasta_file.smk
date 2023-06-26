rule download_genome_fasta_file:
    output:
        fasta_path
    params:
        url=config["genome"]["fasta_link"]
    singularity: str(container_folder.joinpath("alignment.sif"))
    log: log_folder.joinpath("download/genome/fasta.log")
    threads: 8
    shell:
        """
        set -x
        URL={params.url}
        OUTPUT={output}
        [[ ${{URL: -3}} == ".gz" && ! ${{OUTPUT: -3}} == ".gz" ]] &&  \
        ( wget -O - $URL | pigz -d -p {threads} > $OUTPUT |& tee {log} ) || \
        ( wget -O $OUTPUT $URL |& tee {log} )       
        """
