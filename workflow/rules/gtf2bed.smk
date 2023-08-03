rule gtf2bed:
    input:
        gtf=rmsk_path,
    output:
        bed=rmsk_bed,
    threads: 1
    resources: 
        runtime=60,
        mem_mb=4000
    log:
        log_folder.joinpath("samtools_view/bed2gtf.log"),
    conda:
        "../env/samtools.yml"
    shell:
        """
        IN={input}
        if [[ "$IN" = *.gz ]]; then
            T=$(mktemp -u)
            gunzip -c $IN > $T
        else
            T="$IN"
        fi
        
        gtfToGenePred $T temp.genePred
        genePredToBed temp.genePred tmp.bed 
        
        OUT={output}
        if [[ $OUT = *.gz ]]; then
            gzip -c tmp.bed > $OUT
            rm tmp.bed            
        else
            mv tmp.bed $OUT
        fi
        
        rm temp.genePred
        """
