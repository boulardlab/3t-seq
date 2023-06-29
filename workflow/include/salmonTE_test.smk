rule salmonTE_test:
    input:
        salmonTE_folder.joinpath("quant/{serie}"),
        salmonTE_folder.joinpath("quant/{serie}/edit_condition.done"),
    output:
        directory(salmonTE_folder.joinpath("de_analysis/{serie}")),
    log:
        log_folder.joinpath("salmonTE/{serie}/test_de.log"),
    singularity:
        "docker://registry.embl.de/tabaro/snakemake-rnaseq/salmontTE:latest"
    shell:
        """
        python /opt/SalmonTE/SalmonTE.py test \
        --inpath={input[0]} \
        --outpath={output} \
        --tabletype=csv \
        --figtype=png \
        --analysis_type=DE \
        --conditions=control,treatment
        """
