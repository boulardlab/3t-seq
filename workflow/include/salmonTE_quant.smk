rule salmonTE_quant:
    input:
        raw_reads_folder.joinpath("{serie}")
    output:
        directory(data_folder.joinpath("salmonTE/quant/{serie}")),
        salmonTE_folder.joinpath("quant/{serie}/EXPR.csv"),
        salmonTE_folder.joinpath("quant/{serie}/MAPPING_INFO.csv"),
        salmonTE_folder.joinpath("quant/{serie}/clades.csv"),
        salmonTE_folder.joinpath("quant/{serie}/condition.csv")
    params:
          # 13/10/2020 available references: hs mm dr dm
          # https://github.com/LiuzLab/SalmonTE#running-the-quant-mode-to-collect-te-expressions
          # TODO: validate SalmonTE reference genome string
        reference_genome="mm"
    log:
        log_folder.joinpath("salmonTE/{serie}/quant.log")
    singularity:
        str(container_folder.joinpath("SalmonTE.sif"))
    threads:
        8
    shell:
        """
        python /opt/SalmonTE/SalmonTE.py quant \
        --reference={params.reference_genome} \
        --outpath={output[0]} \
        --num_threads={threads} {input}
        """
