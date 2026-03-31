localrules:
    edit_condition_file,


rule edit_condition_file:
    input:
        condition_file=salmonTE_folder.joinpath("quant/{serie}/condition.csv"),
        sample_sheet=get_sample_sheet,
    output:
        touch(data_folder.joinpath("salmonTE/quant/{serie}/edit_condition.done")),
    params:
        variable=lambda wildcards: get_deseq2_variable(wildcards),
        reference_level=lambda wildcards: get_deseq2_reference_level(wildcards),
    log:
        log_folder.joinpath("salmonTE/{serie}/edit_condition.log"),
    conda:
        "../env/pandas.yml"
    threads: 1
    resources:
        runtime=10,
        mem_mb=2048,
    script:
        "../scripts/edit_condition_file.py"


rule salmonTE_quant:
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
    resources:
        runtime=720,
        mem_mb=16000,
    script:
        "../scripts/salmonTE-quant.sh"


rule salmonTE_test:
    input:
        infolder=salmonTE_folder.joinpath("quant/{serie}"),
        condition_file=data_folder.joinpath(
            "salmonTE/quant/{serie}/edit_condition.done"
        ),
    output:
        report(
            directory(salmonTE_folder.joinpath("de_analysis/{serie}")),
            category="SalmonTE",
            subcategory="{serie}",
            patterns=["{name}.png"],
            labels={"serie": "{serie}", "file": "{name}"},
        ),
    log:
        log_folder.joinpath("salmonTE/{serie}/test_de.log"),
    container:
        "docker://ftabaro/salmonte:latest"
    threads: 4
    resources:
        runtime=360,
        mem_mb=16000,
    script:
        "../scripts/salmonTE-test.sh"
