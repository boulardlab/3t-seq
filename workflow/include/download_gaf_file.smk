rule download_gaf_file:
    output:
        gaf_path
    params:
        url=config["genome"]["gaf_link"]
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    log:
        log_folder.joinpath("download/genome/gaf.log")
    threads:
        4
    shell:
        """
        set -x
        URL={params.url}
        OUTPUT={output}
        
        if [ ${{URL: -3}} == ".gz" ]; then
            if [ ${{OUTPUT: -3}} == ".gz" ]; then
                wget -O - $URL | pigz -d -p {threads} | grep -v ! | pigz -c -p {threads} > $OUTPUT |& tee {log}
            else
                wget -O - $URL | pigz -d -p {threads} | grep -v ! > $OUTPUT |& tee {log}
            fi
        else
            if [ ${{OUTPUT: -3}} == ".gz" ]; then
                wget -O - $URL | grep -v ! | pigz -c -p {threads} > $OUTPUT |& tee {log}
            else
                wget -O - $URL | grep -v ! > $OUTPUT |& tee {log}
            fi
        fi 
        """
