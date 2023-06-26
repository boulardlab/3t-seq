# Created by: Francesco Tabaro
# Created on: 10/22/20

input_path <- snakemake@input[[1]]
output_path <- snakemake@output[[1]]

message("input file:", input_path)
message("output file:", output_path)

library(rtracklayer)
ann <- import(input_path)
ann <- subset(ann, type %in% c("gene", "pseudogene"))
saveRDS(ann, output_path)
