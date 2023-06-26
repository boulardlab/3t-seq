rule subset_gtf:
    input:
        gtf_path,
    output:
        rdata_folder.joinpath("deseq2/ann.rds")
    singularity:
        str(container_folder.joinpath("R.sif"))
    log:
        log_folder.joinpath("R/subset_gtf.log")
    script:
          "../scripts/subset_gtf_v1.R"
