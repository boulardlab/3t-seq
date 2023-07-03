localrules:
    edit_condition_file,


rule edit_condition_file:
    input:
        data_folder.joinpath("salmonTE/quant/{serie}"),
        get_sample_sheet,
    output:
        touch(data_folder.joinpath("salmonTE/quant/{serie}/edit_condition.done")),
    log:
        log_folder.joinpath("salmonTE/{serie}/edit_condition.log"),
    run:
        import pandas as pd
        from pathlib import Path

        serie_folder = Path(input[0])
        sample_sheet_path = Path(input[1])

        salmon_condition_sheet = pd.read_csv(serie_folder.joinpath("condition.csv"))
        salmon_condition_sheet = salmon_condition_sheet.set_index("SampleID")

        sample_sheet = pd.read_csv(sample_sheet_path)

        if "filename" in sample_sheet.columns:
            sample_sheet = sample_sheet.set_index("filename")
        else:
            if all(sample_sheet.filename_1.apply(lambda x: x.endswith("_sequence"))):
                sample_sheet["tmp"] = sample_sheet.filename_1.apply(
                    lambda x: x.split("_sequence")[0]
                )
                sample_sheet = sample_sheet.set_index("tmp")
            else:
                sample_sheet = sample_sheet.set_index("filename_1")

        sample_sheet = sample_sheet.reindex(index=salmon_condition_sheet.index)

        joined = sample_sheet.join(salmon_condition_sheet)
        joined["condition"] = joined.apply(
            lambda row: "control" if row["genotype"] == "WT" else "treatment", axis=1
        )
        joined = joined.reset_index()

        out = joined[["SampleID", "condition"]]
        out.to_csv(serie_folder.joinpath("condition.csv"), index=False)


checkpoint salmonTE_quant:
    input:
        raw_reads_folder.joinpath("{serie}"),
    output:
        directory(salmonTE_folder.joinpath("quant/{serie}")),
        salmonTE_folder.joinpath("quant/{serie}/EXPR.csv"),
        salmonTE_folder.joinpath("quant/{serie}/MAPPING_INFO.csv"),
        salmonTE_folder.joinpath("quant/{serie}/clades.csv"),
        salmonTE_folder.joinpath("quant/{serie}/condition.csv"),

    params:
        # 13/10/2020 available references: hs mm dr dm
        # https://github.com/LiuzLab/SalmonTE#running-the-quant-mode-to-collect-te-expressions
        # TODO: validate SalmonTE reference genome string
        reference_genome="mm",
    log:
        log_folder.joinpath("salmonTE/{serie}/quant.log"),
    container:
        "docker://registry.git.embl.de/tabaro/snakemake-rna-seq/salmonte:latest"
    threads: 8
    shell:
        """
        python /opt/SalmonTE/SalmonTE.py quant \
        --reference={params.reference_genome} \
        --outpath={output[0]} \
        --num_threads={threads} {input} |& \
        tee {log}
        """


rule salmonTE_test:
    input:
        salmonTE_folder.joinpath("quant/{serie}"),
        salmonTE_folder.joinpath("quant/{serie}/edit_condition.done"),
    output:
        directory(salmonTE_folder.joinpath("de_analysis/{serie}")),
    log:
        log_folder.joinpath("salmonTE/{serie}/test_de.log"),
    container:
        "docker://registry.git.embl.de/tabaro/snakemake-rna-seq/salmonte:latest"
    shell:
        """
        python /opt/SalmonTE/SalmonTE.py test \
        --inpath={input[0]} \
        --outpath={output} \
        --tabletype=csv \
        --figtype=png \
        --analysis_type=DE \
        --conditions=control,treatment |& \
        tee {log}
        """
