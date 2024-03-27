library(DESeq2)

dds_rds_path <- snakemake@output[["dds"]]
deg_table_path <- snakemake@output[["deg_table"]]

count_matrix <- read.table(
    snakemake@input[["counts"]], 
    sep="\t", 
    header=TRUE,
    comment = "#",
    row.names = "Geneid",
    check.names = FALSE
)
count_matrix <- count_matrix[,-which(colnames(count_matrix) %in% c("Chr", "Start", "End", "Strand", "Length"))]
colnames(count_matrix) <- basename(colnames(count_matrix))
colnames(count_matrix) <- sub(".TEonly.bam", "", colnames(count_matrix))

print(head(count_matrix))

count_matrix <- as.matrix(count_matrix)

print(head(count_matrix))

sample_sheet <- read.csv(
    snakemake@input[["sample_sheet"]], 
    sep=",", 
    header=TRUE,
    row.names = "name"    
)


if ("filename_1" %in% colnames(sample_sheet)) {
    colnames_order <- sapply(colnames(count_matrix), grep, x = sample_sheet$filename_1 )
}else{
    colnames_order <- sapply(colnames(count_matrix), grep, x = sample_sheet$filename)
}

colnames(count_matrix)[colnames_order] <- rownames(sample_sheet)


sample_sheet <- sample_sheet[match(colnames(count_matrix), rownames(sample_sheet)),]

design_variable <- snakemake@params[["variable"]]
reference_level <- snakemake@params[["reference_level"]]

if (!design_variable %in% colnames(sample_sheet)) {
  message <- sprintf("Could not find design variable in columns of sample_sheet.\nvariable: %s\nsample sheet: %s", design_variable, snakemake@input[["sample_sheet"]])
  stop(message)
}


design_formula <- as.formula(sprintf("~ %s", design_variable))

sample_sheet[,design_variable] <- as.factor(sample_sheet[,design_variable])
sample_sheet[,design_variable] <- relevel(sample_sheet[,design_variable], ref = reference_level)

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_sheet,
    design = design_formula
)

dds <- dds[rowSums(counts(dds)) > ncol(dds), ]
dds <- DESeq(dds, test = "Wald")

saveRDS(dds, file = dds_rds_path)

## Extract results
results <- results(dds, name = resultsNames(dds)[2])

results$gene_name <- rownames(results)

columns <- colnames(results)
columns <- c("gene_name", columns[-which(columns == "gene_name")])
results <- results[,columns]


# saveRDS(results, file = results_rds_path)
write.csv(results, file = deg_table_path)

