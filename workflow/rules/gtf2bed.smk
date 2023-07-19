rule gtf2bed:
    input:
        gtf=rmsk_path,
    output:
        bed=rmsk_bed,
    threads: 1
    log:
        log_folder.joinpath("samtools_view/bed2gtf.log"),
    conda:
        "../env/samtools.yml"
    shell:
        """
        IN={input}
        if [[ ${{IN: -3}} == ".gz" ]]; then
            gunzip -c $IN > temp.gtf
        else
            ln -s $IN temp.gtf
        fi
        
        gtfToGenePred temp.gtf temp.genePred
        genePredToBed temp.genePred tmp.bed 
        
        OUT={output}
        if [[ ${{OUT: -3}} == ".gz" ]]; then
            gzip -c tmp.bed > $OUT
            rm tmp.bed            
        else
            mv tmp.bed $OUT
        fi        
        
        rm temp.gtf temp.genePred
        """
