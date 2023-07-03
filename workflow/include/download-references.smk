rule download_genome_fasta_file:
    output:
        fasta_path,
    params:
        url=config["genome"]["fasta_url"],
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("download/genome/fasta.log"),
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


rule download_genome_annotation_file:
    output:
        gtf_path,
    params:
        url=config["genome"]["gtf_url"],
        tmp=config["globals"]["tmp_folder"],
    conda:
        "../env/alignment.yml"
    log:
        log_folder.joinpath("download/genome/gtf.log"),
    threads: 8
    shell:
        """
        set -x 
        URL={params.url}
        TMP=$(mktemp -p {params.tmp})
        
        OUTPUT={output}
        if [ ${{URL: -3}} == ".gz" ] && [ ! ${{OUTPUT: -3}} == ".gz" ]; then
            wget -q -O - $URL | pigz -d -p {threads} > $TMP |& tee {log}
        else
            wget -q -O $TMP $URL |& tee {log}
        fi
        
        N=$(grep -v '#' $TMP | head -n 10 | grep -c 'chr')
        
        if [ ! $N == 10 ]; then
            echo "Adding \"chr\" to first column, then mv'ing to $OUTPUT"
            awk '{{if("#"~\$0){{print $0}}else{{print "chr"$0}}}}' $TMP > $OUTPUT
        else
            echo "Mv'ing to $OUTPUT"
            mv $TMP $OUTPUT
        fi
        """


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


checkpoint download_gtRNAdb:
    output:
        directory(tRNA_annotation_dir),
    params:
        url=config["genome"]["gtrnadb_url"]
    conda:
        "../env/alignment.yml"
    shell:
        """
        mkdir -p {output}
        cd {output}
        F=$(basename {params.url})
        wget {params.url}
        tar xf $F
        """

rule download_gaf_file:
    output:
        gaf_path,
    params:
        url=config["genome"]["gaf_url"],
    conda:
        "../env/alignment.yml"
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
