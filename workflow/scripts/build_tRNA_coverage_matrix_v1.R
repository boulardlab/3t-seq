library(tidyverse)
coverage_files <- snakemake@input
coverage_files %>%
  map(read_tsv, col_names = FALSE) %>%
  map(select, c(X4, X14)) %>%
  purrr::reduce(inner_join, by = "X4") %>%
  rename_with(function(x) c("Name", basename(as.character(coverage_files)))) %>%
  write_tsv(file = as.character(snakemake@output))
