library(DESeq2)

dds_rds_path <- snakemake@output[["dds"]]
deg_table_path <- snakemake@output[["deg_table"]]

count_matrix <- read.table(
    snakemake@input[["counts"]], 
    sep="\t", 
    header=TRUE
)

# Remove duplicated tRNA names.
# Sum multiple instances with the same name
# Assign unique row names to the count matrix
grouping_variable <- as.factor(count_matrix$Name)
cm <- count_matrix[,-which(colnames(count_matrix) == "Name")]
count_matrix <- aggregate(cm, list(Name = grouping_variable), sum)
rownames(count_matrix) <- count_matrix$Name
count_matrix <- count_matrix[,-which(colnames(count_matrix) == "Name")]

count_matrix <- as.matrix(count_matrix)
colnames(count_matrix) <- sub(".bed", "", colnames(count_matrix))

sample_sheet <- read.csv(
    snakemake@input[["sample_sheet"]], 
    sep=",", 
    header=TRUE,
    row.names = "sample"    
)
sample_sheet <- sample_sheet[match(colnames(count_matrix), rownames(sample_sheet)),]

design_variable <- snakemake@params[["variable"]]
reference_level <- snakemake@params[["reference_level"]]

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
results$enrichment <- -log10(results$padj)
results$Name <- rownames(results)

# saveRDS(results, file = results_rds_path)
write.csv(results, file = deg_table_path)

