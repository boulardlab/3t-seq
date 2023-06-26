def get_deseq2_test(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    return deseq2_params["test"]

def get_deseq2_variable(wildcards):
    deseq2_params = get_params(wildcards, "deseq2")
    return deseq2_params["variable"]


rule deseq2:
    input:
         star_flag=star_folder.joinpath("{serie}.done"),
         #counts_folder=star_folder.joinpath("{serie}"),
         # runinfo_file=data_folder.joinpath("{serie}_sra.csv"),
         annotation_file=gtf_path,
         sample_sheet=get_sample_sheet
    output:
          # samples_clustering=pictures_folder.joinpath("deseq2/{serie}/samples_distance.pdf"),
          # pca_plot=pictures_folder.joinpath("deseq2/{serie}/samples_pca.pdf"),
          # deg_heatmap=pictures_folder.joinpath("deseq2/{serie}/deg_heatmap.pdf"),
        touch(analysis_folder.joinpath("deseq2-{serie}.done"))
    params:
        deg_table = str(tables_folder.joinpath("deseq2/{serie}/results.csv")),
        deg_table_shrink = str(tables_folder.joinpath("deseq2/{serie}/results.shrink.csv")),
        dds = str(rdata_folder.joinpath("deseq2/{serie}/dds.rds")),
        # quantile_threshold=0.25,
        annotation_type="mgi",
        test = lambda wildcards: get_deseq2_test(wildcards),
        variable = lambda wildcards: get_deseq2_variable(wildcards)
    singularity:
        str(container_folder.joinpath("R.sif"))
    log:
        log_folder.joinpath("R/{serie}/deseq2.log")
    script:
          "../scripts/deseq2_v1.R"
    # shell:
    #       "R CMD BATCH --vanilla '--args counts_folder={input.counts_folder} annotation_file_path={input.annotation_file} sample_sheet={input.sample_sheet} dds_rds_path={output.dds} deg_table_path={output.deg_table}' ../../src/R/deseq2_v1.R"
