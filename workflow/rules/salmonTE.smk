localrules:
    edit_condition_file,


rule edit_condition_file:
    input:
        data_folder.joinpath("salmonTE/quant/{serie}"),
        get_sample_sheet,
    output:
        touch(data_folder.joinpath("salmonTE/quant/{serie}/edit_condition.done")),
    params:
        variable=lambda wildcards: get_deseq2_variable(wildcards),
    log:
        log_folder.joinpath("salmonTE/{serie}/edit_condition.log"),
    conda:
        "../../env/pandas.yml"
    script:
        "../scripts/edit_condition_file.py"


checkpoint salmonTE_quant:
    input:
        get_salmonTE_quant_input,
    output:
        outfolder=directory(salmonTE_folder.joinpath("quant/{serie}")),
        expr=salmonTE_folder.joinpath("quant/{serie}/EXPR.csv"),
        mapping_info=salmonTE_folder.joinpath("quant/{serie}/MAPPING_INFO.csv"),
        clades=salmonTE_folder.joinpath("quant/{serie}/clades.csv"),
        condition=salmonTE_folder.joinpath("quant/{serie}/condition.csv"),
    params:
        # 13/10/2020 available references: hs mm dr dm
        # https://github.com/LiuzLab/SalmonTE#running-the-quant-mode-to-collect-te-expressions
        reference_genome=lambda w: set_salmonTE_genome(),
    log:
        log_folder.joinpath("salmonTE/{serie}/quant.log"),
    container:
        "docker://ftabaro/salmonte:latest"
    threads: 8
    shell:
        """
        set -x 
        I=""
        T=$(mktemp -d)
        for F in {input}; do
            BN=$(basename $F)
            if [[ $BN == *.gz ]]; then
                gunzip -c $F > $T/${{BN%.gz}}
                I="$I $T"
            else
                I="$I $F"
            fi
        done

        python /opt/SalmonTE/SalmonTE.py quant \
        --reference={params.reference_genome} \
        --outpath={output.outfolder} \
        --num_threads={threads} $I |& \
        tee {log}
        """


rule salmonTE_test:
    input:
        infolder=salmonTE_folder.joinpath("quant/{serie}"),
        condition_file=salmonTE_folder.joinpath("quant/{serie}/edit_condition.done"),
    output:
        directory(salmonTE_folder.joinpath("de_analysis/{serie}")),
    log:
        log_folder.joinpath("salmonTE/{serie}/test_de.log"),
    container:
        "docker://ftabaro/salmonte:latest"
    shell:
        """
        python /opt/SalmonTE/SalmonTE.py test \
        --inpath={input.infolder} \
        --outpath={output} \
        --tabletype=csv \
        --figtype=png \
        --analysis_type=DE \
        --conditions=control,treatment |& \
        tee {log}
        """
