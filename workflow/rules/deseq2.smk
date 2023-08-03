rule subset_gtf:
    input:
        gtf_path,
    output:
        rdata_folder.joinpath("deseq2/ann.rds"),
    conda:
        "../env/R.yml"
    log:
        log_folder.joinpath("R/subset_gtf.log"),
    script:
        "../scripts/subset_gtf_v1.R"


rule deseq2:
    input:
        star_flag=star_folder.joinpath("{serie}.done"),
        #counts_folder=star_folder.joinpath("{serie}"),
        # runinfo_file=data_folder.joinpath("{serie}_sra.csv"),
        annotation_file=gtf_path,
        sample_sheet=get_sample_sheet,
    output:
        # samples_clustering=pictures_folder.joinpath("deseq2/{serie}/samples_distance.pdf"),
        # pca_plot=pictures_folder.joinpath("deseq2/{serie}/samples_pca.pdf"),
        # deg_heatmap=pictures_folder.joinpath("deseq2/{serie}/deg_heatmap.pdf"),
        touch(analysis_folder.joinpath("deseq2-{serie}.done")),
    params:
        deg_table=str(tables_folder.joinpath("deseq2/{serie}/results.csv")),
        deg_table_shrink=str(
            tables_folder.joinpath("deseq2/{serie}/results.shrink.csv")
        ),
        dds=str(rdata_folder.joinpath("deseq2/{serie}/dds.rds")),
        # quantile_threshold=0.25,
        annotation_type=config["genome"]["annotation_type"],
        test=lambda wildcards: get_deseq2_test(wildcards),
        variable=lambda wildcards: get_deseq2_variable(wildcards),
    conda:
        "../env/R.yml"
    log:
        log_folder.joinpath("R/{serie}/deseq2.log"),
    script:
        "../scripts/deseq2_v1.R"


# rule deseq2_report:
#     input:
#         analysis_folder.joinpath("deseq2-{serie}.done"),
#         gaf_path,
#         ann_path=str(rdata_folder.joinpath("deseq2/ann.rds")),
#     output:
#         touch(analysis_folder.joinpath("deseq2-report-{serie}.done")),
#         # notebooks_folder.joinpath("{serie}/%s.html"%get_filename(deseq2_notebook_input_path, stem=True))
#     params:
#         dds_path=str(rdata_folder.joinpath("deseq2/{serie}/dds.rds")),
#         # ann_path = str(rdata_folder.joinpath("deseq2/ann.rds")),
#         notebook_path=deseq2_notebook_input_path,
#         # root_path = deseq2_working_directory,
#         output_dir=notebooks_folder,
#         annotation_type=config["genome"]["annotation_type"],
#     conda:
#         "../env/R.yml"
#     log:
#         log_folder.joinpath("R/{serie}/deseq2_report.log"),
#     shell:
#         """
#          set -e 
#          if [ -f {params.dds_path} ]; then
#             OUTPUT_DIR="{params.output_dir}/{wildcards.serie}"
#             R --no-echo --vanilla \
#             -e "rmarkdown::render('{params.notebook_path}', output_dir = '$OUTPUT_DIR', output_format = 'html_document', params = list(output_dir = '$OUTPUT_DIR', ann_path = '{input.ann_path}', dds_path = '{params.dds_path}', annotation_type = '{params.annotation_type}', gaf_path = '{input[1]}'))" |& \
#             tee {log}
#         else
#             exit 0
#         fi
#          """
