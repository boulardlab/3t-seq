from utilities.general import get_filename

rule deseq2_report:
    input:
        analysis_folder.joinpath("deseq2-{serie}.done"),
        gaf_path,
        ann_path = str(rdata_folder.joinpath("deseq2/ann.rds"))
    output:
        touch(analysis_folder.joinpath("deseq2-report-{serie}.done"))
        # notebooks_folder.joinpath("{serie}/%s.html"%get_filename(deseq2_notebook_input_path, stem=True))
    params:
        dds_path = str(rdata_folder.joinpath("deseq2/{serie}/dds.rds")),
        # ann_path = str(rdata_folder.joinpath("deseq2/ann.rds")),
        notebook_path = deseq2_notebook_input_path,
        # root_path = deseq2_working_directory,
        output_dir = notebooks_folder,
        annotation_type = config["genome"]["annotation_type"]
    singularity:
        str(container_folder.joinpath("R.sif"))
    log:
        log_folder.joinpath("R/{serie}/deseq2_report.log")
    shell:
         """
         set -x
         if [ -f {params.dds_path} ]; then
            OUTPUT_DIR="{params.output_dir}/{wildcards.serie}"
            R --no-echo --vanilla \
            -e "rmarkdown::render('{params.notebook_path}', output_dir = '$OUTPUT_DIR', output_format = 'html_document', params = list(output_dir = '$OUTPUT_DIR', ann_path = '{input.ann_path}', dds_path = '{params.dds_path}', annotation_type = '{params.annotation_type}', gaf_path = '{input[1]}'))" |& \
            tee {log}
        else
            exit 0
        fi
         """

# knit_root_dir = '{params.root_path}'

