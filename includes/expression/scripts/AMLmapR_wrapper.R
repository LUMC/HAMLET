#!/usr/bin/env Rscript

# seAMLess Analysis Pipeline
# Refactored for command-line usage using only base R

suppressPackageStartupMessages({
  library(AMLmapR)
  library(caret)
})


# Helper: verbose printing
verbose_print <- function(..., verbose = FALSE) {
  if (verbose) message(paste0(...))
}


scale_data <- function(matrix,d){
  matrix <- predict(d$scaler, matrix[,d$genes])
}

#'
deseq_normalise <- function(matrix,d){
  pseudo_reference <- d[["keep"]][[2]]
  keep <- d$keep[[1]]
  ratio_to_ref  <- apply(matrix[,keep],1, function(x) x/pseudo_reference)
  sizeFactor  <- apply(ratio_to_ref,2, function(x) stats::median(x))
  matrix <- matrix/sizeFactor
  matrix <- log(matrix + 1)
  matrix
}

#' @import kernlab
#'
pred_clusters <- function(matrix,d){
  predictions <- lapply(d$models, function(model){
    predictions <- -1 * predict(model, newdata = matrix, type = "decision")
  })
  predictions <- do.call(cbind, predictions)
  colnames(predictions) <- names(d$models)
  predictions <- data.frame(predictions)
  predictions$prediction <- apply(predictions, 1, function(x) colnames(predictions)[which.max(x)])
  predictions$pass_cutoff <- apply(predictions[,1:length(d$models)], 1, function(x) max(x) > -0.04)
  predictions
}

# Print usage
usage <- function() {
  cat("
Usage: Rscript pipeline.R \\
  -c <counts_file> [-c <counts_file> ...] \\
  -o <outdir> \\
  [-s <sample_name>] \\
  [-v]

  -c, --counts        STAR count file(s) to process (can specify multiple)
  -o, --outdir        Output directory prefix for result files
  -s, --sample-name   Name of the input sample (defaults to basename of each file)
  -v, --verbose       Enable verbose logging

")
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
counts <- character()
outdir <- NULL
sample_name <- NULL
verbose <- FALSE

i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (arg %in% c("-c", "--counts")) {
    if (i == length(args)) stop("No value provided for ", arg)
    counts <- c(counts, args[i + 1])
    i <- i + 2
  } else if (arg %in% c("-o", "--outdir")) {
    if (i == length(args)) stop("No value provided for ", arg)
    outdir <- args[i + 1]
    i <- i + 2
  } else if (arg %in% c("-s", "--sample-name")) {
    if (i == length(args)) stop("No value provided for ", arg)
    sample_name <- args[i + 1]
    i <- i + 2
  } else if (arg %in% c("-v", "--verbose")) {
    verbose <- TRUE
    i <- i + 1
  } else {
    stop("Unknown argument: ", arg)
  }
}

# Validate required inputs
if (length(counts) == 0 || is.null(outdir)) {
  usage()
  stop("Both -c/--counts and -o/--outdir must be provided.", call. = FALSE)
}

# Create output directory if needed
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Function to read STAR count files
read_star <- function(file) {
  first_line <- readLines(file, n = 1)
  skip_lines <- if (grepl("^N_unmapped", first_line)) 4 else 0
  d <- read.table(file, header = FALSE, skip = skip_lines, stringsAsFactors = FALSE)
  sn <- if (!is.null(sample_name)) sample_name else basename(file)
  df <- setNames(d[, 1:2], c("GeneID", sn))
  return(df)
}

verbose_print("Reading STAR count files...\n", verbose = verbose)
count_list <- lapply(counts, read_star)
merged_counts <- Reduce(function(a, b) merge(a, b, by = "GeneID", all = TRUE), count_list)
merged_counts[is.na(merged_counts)] <- 0

# Prepare matrix input
genecol <- "GeneID"
mat_input <- merged_counts[, c(genecol, setdiff(names(merged_counts), genecol))]
rownames(mat_input) <- mat_input$GeneID
mat_input$GeneID <- NULL
mat <- t(as.matrix(mat_input))

# Load example matrix from AMLmapR
example_matrix <- AMLmapR::example_matrix

# Align genes by stripping version suffixes
original_genes_mat     <- colnames(mat)
original_genes_example <- colnames(example_matrix)
stripped_mat     <- sub("\\..*$", "", original_genes_mat)
stripped_example <- sub("\\..*$", "", original_genes_example)

common_genes <- intersect(stripped_mat, stripped_example)
keep_cols <- which(stripped_mat %in% common_genes)
mat <- mat[, keep_cols, drop = FALSE]
stripped_mat <- stripped_mat[keep_cols]

# Map stripped names to full example names
example_map <- setNames(original_genes_example, stripped_example)

# Build a new matrix with example_matrix columns
new_mat <- matrix(0,
                  nrow = nrow(mat),
                  ncol = length(original_genes_example),
                  dimnames = list(rownames(mat), original_genes_example))

for (i in seq_along(stripped_mat)) {
  gene <- stripped_mat[i]
  if (gene %in% names(example_map)) {
    col_target <- example_map[[gene]]
    new_mat[, col_target] <- mat[, i]
  }
}

mat <- new_mat  # final input matrix
# Ensure values are non-negative
new_mat[new_mat < 0] <- 0

# Round to nearest integer
new_mat <- round(new_mat)

# Convert to integer matrix with correct dimensions and names
mat <- matrix(as.integer(new_mat), nrow = nrow(new_mat), dimnames = dimnames(new_mat))


extra_sample_row <- example_matrix[1, , drop = FALSE]
mat_two_sample <- rbind(extra_sample_row, mat)

# Run prediction
verbose_print("Running AML cluster prediction...\n", verbose = verbose)
predictions <- predict_AML_clusters(mat_two_sample)


# Remove temporary sample before saving
predictions <- predictions[rownames(predictions) != rownames(extra_sample_row), , drop = FALSE]

# Convert list columns to character and replace NaN in numeric columns
for (col in names(predictions)) {
  if (is.list(predictions[[col]])) {
    predictions[[col]] <- sapply(predictions[[col]], function(x) {
      if (length(x) == 1) return(as.character(x))
      else return(paste(as.character(x), collapse = ";"))
    })
  } else if (is.numeric(predictions[[col]])) {
    predictions[[col]][is.nan(predictions[[col]])] <- 0
  }
}

# Fill empty 'prediction' entries with "NotClassified"
if ("prediction" %in% names(predictions)) {
  predictions$prediction[predictions$prediction == ""] <- "NotClassified"
}

# Write to CSV without row names
outfile <- file.path(outdir, "aml_cluster_predictions.csv")
write.csv(predictions, file = outfile, row.names = FALSE)
verbose_print(paste0("Predictions saved to: ", outfile, "\n"), verbose = verbose)
