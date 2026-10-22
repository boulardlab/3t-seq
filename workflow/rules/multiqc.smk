
rule multiqc:
    input:
        unpack(get_multiqc_inputs),
    output:
        report(
            multiqc_folder.joinpath("{serie}", "multiqc_report.html"),
            category="MultiQC",
            labels={"serie": "{serie}"},
        ),
    log:
        log_folder.joinpath("multiqc", "multiqc-{serie}.log"),
    conda:
        "../env/qc.yml"
    resources:
        runtime=get_multiqc_runtime,
        mem_mb=2048,
    params:
        config=workflow.source_path("../config/multiqc_config.yaml"),
        outdir=lambda wc: multiqc_folder / wc.serie,
        title=lambda wc: wc.serie,
    shell:
        """
        set -e
        multiqc --force --dirs --export \
          --config {params.config} \
          --cl-config "title: '{params.title}'" \
          -o {params.outdir} \
          {input} |& tee {log}
        """
