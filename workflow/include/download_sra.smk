checkpoint download_sra:
    input:
        # data_folder.joinpath("{serie}/{sample}.txt")
        data_folder.joinpath("{serie}_sra.csv")
    output:
        # dynamic(raw_reads_folder.joinpath("{serie}/{sample}_{mate}.fastq"))
        directory(raw_reads_folder.joinpath("{serie}"))
    singularity:
        # paths to singularity images cannot be PosixPaths.
        str(container_folder.joinpath("ncbi.sif"))
    threads: 8
    log:
        # log_folder.joinpath("{serie}/download_{sample}_{mate}.log")
        log_folder.joinpath("download/{serie}/sra.log")
    shell:
        """
        set -x
        cat {input} | cut -d "," -f 1 | grep SRR | head | \
        parallel -j {threads} fastq-dump -I -X 1000 -v -L info --split-files --skip-technical -O {output} {{}} |& \
        tee {log}
        """
