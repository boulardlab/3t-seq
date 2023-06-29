checkpoint get_runinfo:
    output:
        data_folder.joinpath("{serie}_sra.csv"),
    conda:
        "../env/ncbi.yml"
    log:
        log_folder.joinpath("download/{serie}/runinfo.log"),
    shell:
        """ 
        set -x
        esearch -db gds -query {wildcards.serie}[ACCN] 2>{log} | \
        elink -name gds_sra -db gds -target sra 2>>{log} | \
        efetch -db sra -format runinfo 2>>{log} >{output} 
        """
