sink(file = stderr())

library(tidyverse)

# Inputs from snakemake
count_files <- as.character(snakemake@input[["mimseq_counts"]])
sample_sheet <- read_csv(as.character(snakemake@input[["sample_sheet"]]), show_col_types = FALSE)

# Extract the sample name from the file path. Assuming files end in '_cluster_counts.txt'
path_names <- gsub("_cluster_counts.txt|\\.tsv|\\.txt", "", basename(count_files))

# Create a mapping between paths and names
count_files <- setNames(count_files, path_names)

# Sort paths based on the name column of the sample sheet to ensure order matches experimental design
count_files <- count_files[match(sample_sheet$name, names(count_files))]

# Read each mim-tRNA-seq output file and extract the expected columns
# mimseq output typically contains columns like 'Name', 'Length', and 'Count'
out <- count_files %>%
  map(read_tsv, show_col_types = FALSE) %>%
  map(dplyr::select, c("Name", "Count"))

ncol <- length(out)
nrow <- nrow(out[[1]]) # assume that all elements in out have the same number of rows due to using the same index

# Initialize matrix
m <- matrix(0, nrow, ncol, dimnames = list(out[[1]]$Name, names(count_files)))

# Populate matrix
for (i in seq_along(out)) {
    df <- out[[i]]
    if (nrow(df) != nrow(m)) {
        stop(sprintf("Incorrect dimensions! Element %d has %d rows instead of %d.", i, nrow(df), nrow(m)))
    }
    m[,i] <- df$Count
}

# Convert back to tibble and set the index column name to "Name" to perfectly mimic legacy output
m <- as_tibble(m, rownames="Name")

write_tsv(x = m, file = as.character(snakemake@output))
