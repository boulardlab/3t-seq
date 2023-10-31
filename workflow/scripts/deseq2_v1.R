options(error = quote({
  fn <- "deseq2_emergency_dump"
  dump.frames(dumpto = fn, to.file = TRUE, include.GlobalEnv = TRUE)
  cat(sprintf("Emergency dump saved in %s.rds", fn), file = stderr())
  quit(status = 1)
}))

library(data.table)
library(rtracklayer)
library(DESeq2)

## Get file paths from Snakemake
## INPUTS
# counts_folder <- snakemake@input[["counts_folder"]]
star_flag <- snakemake@input[["star_flag"]]
counts_folder <- sub(".done", "", star_flag)
annotation_file_path <- snakemake@input[["annotation_file"]]
sample_sheet <- snakemake@input[["sample_sheet"]]


# OUTPUTS
dds_rds_path <- snakemake@output[["dds"]]
deg_table_path <- snakemake@output[["deg_table"]]
deg_table_shrink_path <- snakemake@output[["deg_table_shrink"]]

if (!dir.exists(dirname(dds_rds_path))) {
  dir.create(dirname(dds_rds_path), recursive = TRUE)
}

if (!dir.exists(dirname(deg_table_path))) {
  dir.create(dirname(deg_table_path), recursive = TRUE)
}

annotation_type <- snakemake@params[["annotation_type"]]
test_name <- snakemake@params[["test"]]
design_variable <- snakemake@params[["variable"]]
reference_level <- snakemake@params[["reference_level"]]

if (length(design_variable) == 1) {
  design_formula <- as.formula(sprintf("~ %s", design_variable))
} else if (length(design_variable) == 2) {
  design_formula <- as.formula(sprintf("~ %s", paste(design_variable, collapse = " * "))) # test for a + b + a:b
} else {
  stop("Detected more than 2 design variables. This feature is currently not supported.")
}

if (test_name == "LRT") {
  if (length(design_variable) == 1) {
    reduced_model <- as.formula("~ 1")
  } else if (length(design_variable) == 2) {
    reduced_model <- as.formula(sprintf("~ %s", paste(design_variable, collapse = " + "))) # remove the interaction term (a:b)
  }
}

###

## Import files
### GTF
ann <- import(annotation_file_path)

### Sample sheet
colData <- fread(sample_sheet)

if ("filename" %in% colnames(colData)) {
  libtype <- "single"
  if (!all(c("sample", "filename", design_variable) %in% colnames(colData))) {
    stop(
      sprintf("Wrong columns in sample sheet. Colnames must be: \"sample\", \"filename\" and \"%s\".\nRemove extensions from filename column. Detected column names: %s", design_variable, colnames(colData))
    )
  }
} else if ("filename_1" %in% colnames(colData)) {
  libtype <- "paired"
  if (!all(c("sample", "filename_1", "filename_2", design_variable) %in% colnames(colData))) {
    stop(
      sprintf("Wrong columns in sample sheet. Colnames must be: \"sample\", \"filename_1\", \"filename_2\", and \"%s\".\nRemove extensions from filename column. Detected column names: %s", design_variable, colnames(colData)))
  }
} else {
  stop(sprintf("Could not determine sequencing library type. From the following column names: %s\nPlease, check your column names. They should be all lower case.", colnames(colData)))
}

# build colData DataFrame
colData <- DataFrame(colData)
# convert design columns to factors
for (i in seq_along(design_variable)) {
  column <- design_variable[i]
  colData[, column] <- as.factor(colData[, column])
  colData[, column] <- relevel(colData[, column], ref = reference_level)  
}

if (libtype == "single") {
  rownames(colData) <- colData[, "filename"]
  colnames_order <- sapply(colnames(count_matrix), grep, x = sample_sheet$filename)
} else {
  rownames(colData) <- gsub("_1.*$", "", colData[, "filename_1"])
  colnames_order <- sapply(colnames(count_matrix), grep, x = sample_sheet$filename_1 )
}

# Rename columns to match sample sheet sample column
colnames(count_matrix)[colnames_order] <- rownames(sample_sheet)


## Import STAR counts and generate a count matrix
ls <- list.files(counts_folder, pattern = "ReadsPerGene.out.tab", full.names = TRUE)
ls <- setNames(ls, sub("(^.+)\\.ReadsPerGene\\.out\\.tab", "\\1", basename(ls)))
mat <- sapply(ls, function(p) {
  dt <- fread(p, skip = 4)
  setNames(dt$V2, dt$V1)
})

## Extract annotations data for the quantified genes
if (annotation_type == "gencode" || annotation_type == "ensembl") {
  ann <- ann[ann$gene_id %in% rownames(mat) & ann$type == "gene"]
} else if (annotation_type == "mgi") {
  ann <- ann[ann$gene_id %in% rownames(mat) & ann$source == "MGI"]
} else {
  ann <- ann[ann$gene_id %in% rownames(mat)]
}
names(ann) <- ann$gene_id
ann <- ann[match(rownames(mat), names(ann))]

if (!all(rownames(mat) %in% ann$gene_id)) {
  stop("Missing 1:1 mapping between count matrix and gene_annotation. Please check annotation_type in config file and gtf file.")
}

## Match mat colnames to sample_sheet rownames
colData <- colData[match(colnames(mat), rownames(colData)), ]


## Build DESeqDataSet object
se <- SummarizedExperiment(
  assays = SimpleList(counts = mat),
  rowRanges = ann, colData = colData
)
dds <- DESeqDataSet(se, design = design_formula)

# apply light prefiltering (as shown in vignette, quantile filtering may be too harsh)
dds <- dds[rowSums(counts(dds)) > ncol(dds), ]

## Run DESeq2. Export modified dds for later.
if (test_name == "LRT") {
  dds <- DESeq(dds, test = test_name, reduced = reduced_model)
} else { # Wald t-test
  dds <- DESeq(dds, test = test_name)
}

saveRDS(dds, file = dds_rds_path)

## Extract results
results <- results(dds, name = resultsNames(dds)[2])
results_shrink <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")

gene_names_columns <- c("gene_name", "Name")
k <- gene_names_columns %in% colnames(mcols(rowRanges(dds)))
if (any(k)) {
  cn <- gene_names_columns[k]

  df <- subset(mcols(rowRanges(dds)), select = cn)
  df$gene_id <- rownames(dds)

  results$gene_id <- rownames(results)
  results_shrink$gene_id <- rownames(results_shrink)

  results <- merge(as.data.frame(results), as.data.frame(df), by = "gene_id")
  results_shrink <- merge(as.data.frame(results_shrink), as.data.frame(df), by = "gene_id")

  idx <- which(colnames(results) == cn)
  colnames(results)[idx] <- "gene_name"

  idx <- which(colnames(results_shrink) == cn)
  colnames(results_shrink)[idx] <- "gene_name"

  columns <- colnames(results)
  columns <- c("gene_name", columns[-which(columns == "gene_name")])
  results <- results[,columns]

  columns <- colnames(results_shrink)
  columns <- c("gene_name", columns[-which(columns == "gene_name")])
  results_shrink <- results_shrink[,columns]


} else {

  results$gene_name <- rownames(results)
  results_shrink$gene_name <- rownames(results_shrink)
}

write.csv(results, file = deg_table_path)
write.csv(results, file = deg_table_shrink_path)
