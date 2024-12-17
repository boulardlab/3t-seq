sink(file = stderr())

library(tidyverse)
coverage_files <- as.character(snakemake@input[["bed"]])

sample_sheet <- read_csv(as.character(snakemake@input[["sample_sheet"]]))
print(sample_sheet)

# Extract the name component from paths
path_names <- gsub(".bed", "", basename(coverage_files))
print(path_names)
# Create a mapping between paths and names
coverage_files <- setNames(coverage_files, path_names)
print(coverage_files)
# Sort paths based on the name column of the data frame
coverage_files <- coverage_files[match(sample_sheet$name, names(coverage_files))]
print(coverage_files)

out <- coverage_files %>%
  map(read_tsv, col_names = FALSE) %>%
  map(dplyr::select, c(X4, X14))

ncol <- length(out)
nrow <- nrow(out[[1]]) # assume that all elements in out have the same number of rows

m <- matrix(0, nrow, ncol, dimnames = list(out[[1]]$X4, basename(as.character(coverage_files))))
for (i in seq(length(out))) {
    df <- out[[i]]
    if (nrow(df) != nrow(m)) stop(sprintf("Incorrect dimensions! Element %d has %d rows instead of %d.",
        i, nrow(df), nrow(m)))
    #df <- df[match(rownames(m), df$X4),]
    m[,i] <- df$X14
}
m <- as_tibble(m, rownames="Name")

write_tsv(x = m, file = as.character(snakemake@output))
