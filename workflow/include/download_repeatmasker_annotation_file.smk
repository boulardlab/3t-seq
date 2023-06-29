rule download_repeatmasker_annotation_file:
    output:
        rmsk_path,
    params:
        url=config["genome"]["rmsk_link"],
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("download/genome/rmsk.log"),
    threads: 4
    shell:
        """
        set -x
        URL={params.url}
        OUTPUT={output}
        [[ ${{URL: -3}} == ".gz" && ! ${{OUTPUT: -3}} == ".gz" ]] &&  \
        ( wget -O - $URL | pigz -d -p {threads} > $OUTPUT |& tee {log} ) || \
        ( wget -O $OUTPUT $URL |& tee {log} )         
        """
