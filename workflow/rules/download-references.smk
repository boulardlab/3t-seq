rule download_genome_fasta_file:
    output:
        fasta_path,
    params:
        url=config["genome"]["fasta_url"],
    conda:
        "../../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/fasta.log"),
    threads: 1
    shell:
        """
        set -x
        URL="{params.url}"
        OUTPUT={output}
        [[ ${{URL: -3}} == ".gz" && ! ${{OUTPUT: -3}} == ".gz" ]] &&  \
        ( wget -q -O - "$URL" | gunzip -c > $OUTPUT |& tee {log} ) || \
        ( wget -q -O $OUTPUT "$URL" |& tee {log} )       
        """


rule download_genome_annotation_file:
    output:
        gtf_path,
    params:
        url=config["genome"]["gtf_url"],
        tmp=config["globals"]["tmp_folder"],
    conda:
        "../../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/gtf.log"),
    threads: 1
    script:
        "../scripts/download-gtf.sh"


rule download_repeatmasker_annotation_file:
    output:
        rmsk_path,
    params:
        url=config["genome"]["rmsk_link"],
    conda:
        "../../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/rmsk.log"),
    threads: 1
    shell:
        """
        set -x
        URL={params.url}
        OUTPUT={output}
        [[ ${{URL: -3}} == ".gz" && ! ${{OUTPUT: -3}} == ".gz" ]] &&  \
        ( wget -q -O - $URL | gunzip -c > $OUTPUT |& tee {log} ) || \
        ( wget -q -O $OUTPUT $URL |& tee {log} )         
        """


checkpoint download_gtRNAdb:
    output:
        directory(tRNA_annotation_dir),
    params:
        url=config["genome"]["gtrnadb_url"],
    log:
        log_folder.joinpath("download/gtrnadb.log"),    
    conda:
        "../../env/wget.yml"
    shell:
        """
        mkdir -p {output}
        cd {output}
        F=$(basename {params.url})
        wget -q {params.url} |& tee {log}
        tar xvf $F |& tee -a {log}
        """


rule download_gaf_file:
    output:
        gaf_path,
    params:
        url=config["genome"]["gaf_url"],
    conda:
        "../../env/wget.yml"
    log:
        log_folder.joinpath("download/genome/gaf.log"),
    threads: 4
    shell:
        """
        set -x
        URL={params.url}
        OUTPUT={output}
        
        if [ ${{URL: -3}} == ".gz" ]; then
            if [ ${{OUTPUT: -3}} == ".gz" ]; then
                wget -q -O - $URL | zgrep -v ! | pigz -c -p {threads} > $OUTPUT |& tee {log}
            else
                wget -q -O - $URL | zgrep -v ! > $OUTPUT |& tee {log}
            fi
        else
            if [ ${{OUTPUT: -3}} == ".gz" ]; then
                wget -q -O - $URL | grep -v ! | pigz -c -p {threads} > $OUTPUT |& tee {log}
            else
                wget -q -O - $URL | grep -v ! > $OUTPUT |& tee {log}
            fi
        fi 
        """
