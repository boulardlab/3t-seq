sink(file = stderr())

library(tidyverse)
coverage_files <- as.character(snakemake@input[["bed"]])

sample_sheet <- read_csv(as.character(snakemake@input[["sample_sheet"]]))
print(sample_sheet)

# Extract the name component from paths
path_names <- gsub(".*/(\\w+_\\w+).bed", "\\1", coverage_files)

# Create a mapping between paths and names
coverage_files <- setNames(coverage_files, path_names)

# Sort paths based on the name column of the data frame
coverage_files <- coverage_files[match(sample_sheet$name, names(coverage_files))]
print(coverage_files)


out <- coverage_files %>%
  map(read_tsv, col_names = FALSE) %>%
  map(select, c(X4, X14)) %>%
  purrr::reduce(inner_join, by = "X4") %>%
  rename_with(function(x) c("Name", basename(as.character(coverage_files)))) 
  
print(head(out), file = stderr())
  
write_tsv(x = out, file = as.character(snakemake@output))
