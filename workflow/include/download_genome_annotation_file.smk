rule download_genome_annotation_file:
    output:
        gtf_path
    params:
        url = config["genome"]["gtf_link"],
        tmp = config["globals"]["tmp_folder"]
    singularity:
        str(container_folder.joinpath("alignment.sif"))
    log:
       log_folder.joinpath("download/genome/gtf.log")
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