#!/usr/bin/env Rscript

# seAMLess Analysis Pipeline


suppressPackageStartupMessages({
  library(optparse)
  library(seAMLess)
  library(Biobase)
  library(ggplot2)
})


# Helper: verbose printing
verbose_print <- function(..., verbose = FALSE) {
  if (verbose) message(paste0(...))
}

seAMLess_modified <- function(mat,
                     scRef         = seAMLessData::scRef,
                     scRef.sample  = "Sample",
                     scRef.label   = "label.new",
                     verbose       = TRUE) {
  # Printing function
  verbosePrint <- verboseFn(verbose)
  verbosePrint(">> Loading libraries...")

  requireNamespace("randomForest", quietly = TRUE)
  requireNamespace("Biobase",     quietly = TRUE)
  requireNamespace("MuSiC",       quietly = TRUE)


  # wrangle count matrix

  mat <- wrangleMat(mat)

  ## ─── ADD MISSING REFERENCE GENES AS ZERO ROWS ─────────────────────────────
  # get all genes in the single-cell reference
  ref_genes <- rownames(scRef@assayData$exprs)
  # find which of those are absent in our bulk matrix
  missing_genes <- setdiff(ref_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    verbosePrint(">> Adding missing genes from reference: ",
                 paste(missing_genes, collapse = ", "))
    # create a zero matrix for those genes
    zero_mat <- matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(mat),
      dimnames = list(missing_genes, colnames(mat))
    )
    # bind them onto the bottom of our bulk data
    mat <- rbind(mat, zero_mat)
  }
  ## 

  # If ensembl ids are provided
  if (grepl("ENSG", rownames(mat)[1], fixed = TRUE)) {
    verbosePrint(">> Converting human ensembl ids to symbols...")
    # ens to symbol map
    ens2gene <- seAMLess::grch38
    m <- match(rownames(mat), ens2gene$ensgene)
    mapped.genes <- ens2gene$symbol[m]
    removed.genes <- duplicated(mapped.genes) | is.na(mapped.genes)
    mat <- mat[!removed.genes, ]
    rownames(mat) <- mapped.genes[!removed.genes]
  }

  verbosePrint(">> Common genes with the reference: ",
               sum(rownames(mat) %in% rownames(scRef@assayData$exprs)))

  # make mat suitable for MuSiC
  T.eset <- Biobase::ExpressionSet(assayData = as.matrix(mat))

  verbosePrint(">> Deconvoluting samples...")
  # MuSiC deconvolution
  deconv <- MuSiC::music_prop(
    bulk.eset  = T.eset,
    sc.eset    = scRef,
    clusters   = scRef.label,
    markers    = NULL,
    normalize  = FALSE,
    samples    = scRef.sample,
    verbose    = FALSE
  )$Est.prop.weighted
  verbosePrint(">> Deconvolution completed...")

  # --- NEW: harmonize cell‐type set to default scRef ------------------------
  verbosePrint(">> Checking deconvolution cell types against default reference...")
  # 1) default (built-in) cell types
  default_ct <- c("CD14 Mono",   "GMP", "T Cells", "pre B",  "LMPP", 
  "B Cells", "Early Eryth", "EMP", "Late Eryth", "pDC", 
  "CLP", "Plasma", "HSC", "pro B", "cDC", "BaEoMa", "Prog Mk", 
  "pre-pDC", "NK Cells", "pre-mDC", "CD16 Mono", "ASDC")
  # 2) what we actually got
  recv_ct    <- colnames(deconv)

  # 3) drop extras
  extra_ct <- setdiff(recv_ct, default_ct)
  if (length(extra_ct) > 0) {
    verbosePrint(">> Removing extra cell types: ", paste(extra_ct, collapse = ", "))
    deconv <- deconv[, setdiff(recv_ct, extra_ct), drop = FALSE]
  }

  # 4) add missing
  missing_ct <- setdiff(default_ct, recv_ct)
  if (length(missing_ct) > 0) {
    verbosePrint(">> Adding missing cell types: ", paste(missing_ct, collapse = ", "))
    # build a zero‐matrix for all missing types
    zero_mat <- matrix(
      0,
      nrow = nrow(deconv),
      ncol = length(missing_ct),
      dimnames = list(rownames(deconv), missing_ct)
    )
    deconv <- cbind(deconv, zero_mat)
  }

  # 5) reorder columns to the default order
  deconv <- deconv[, default_ct, drop = FALSE]
  # -------------------------------------------------------------------------

  verbosePrint(">> Predicting Venetoclax resistance...")
  veno.res <- stats::predict(seAMLess::venoModel, newdata = deconv)

  return(list(
    Deconvolution         = deconv,
    Venetoclax.resistance = veno.res
  ))
}


# Define command-line options
option_list <- list(
  make_option(c("-e", "--exprs-ref"), type = "character", default = NULL,
              help = "Path to reference expression CSV (genes x samples). If not provided, uses built-in seAMLessData::scRef.", metavar = "file"),
  make_option(c("-m", "--meta-ref"), type = "character", default = NULL,
              help = "Path to reference metadata CSV (samples x features). Required if --exprs-ref is provided.", metavar = "file"),
  make_option(c("-c", "--counts"), type = "character", action = "append", default = NULL,
              help = "STAR count file(s) to process (can specify multiple)", metavar = "file"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory prefix for result files", metavar = "dir"),
  make_option(c("-s", "--sample-name"), type = "character", default = NULL,
              help = "Name of the input sample"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Enable verbose logging")
)

parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

# Validate inputs
if (is.null(opts$c) || is.null(opts$o)) {
  print_help(parser)
  stop("Both --counts and --outdir must be provided.", call. = FALSE)
}

# If custom reference is provided, both exprs-ref and meta-ref must be set
use_custom_ref <- !is.null(opts$e)
if (use_custom_ref) {
  ext <- tolower(tools::file_ext(opts$e))
  if (ext == "rda") {
    # .rda needs no meta-ref
    require_meta <- FALSE
  } else {
    # CSV needs meta-ref
    require_meta <- TRUE
  }
  if (require_meta && is.null(opts$m)) {
    stop("--meta-ref must be provided when --exprs-ref points to a CSV file.", call. = FALSE)
  }
}

# Create output directory if it doesn't exist
if (!dir.exists(opts$o)) {
  dir.create(opts$o, recursive = TRUE)
}



# Load or default reference data
if (use_custom_ref) {
  ext <- tolower(tools::file_ext(opts$e))
  if (ext == "rda") {
    verbose_print("Loading reference from RDA: ", opts$e)
    tmp <- new.env()
    load(opts$e, envir = tmp)
    if (!exists("scRef", envir = tmp)) {
      stop("The .rda file does not contain an object named 'scRef'.", call. = FALSE)
    }
    ref_sub <- tmp$scRef

  } else {
    verbose_print("Loading reference expression (CSV): ", opts$e)
    exprs_ref <- read.csv(opts$e, header = TRUE, row.names = 1)
    exprs_ref <- data.matrix(exprs_ref)
    colnames(exprs_ref) <- sub("^X", "", colnames(exprs_ref))

    verbose_print("Loading reference metadata (CSV): ", opts$m)
    meta_ref <- read.csv(opts$m, header = TRUE, row.names = 1)

    ref_sub <- ExpressionSet(
      assayData   = exprs_ref,
      phenoData   = AnnotatedDataFrame(meta_ref)
    )
  }

} else {
  verbose_print("Using built-in single-cell reference: seAMLessData::scRef")
  suppressPackageStartupMessages(library(seAMLessData))
  ref_sub <- seAMLessData::scRef
}

# Define STAR count reader
read_star <- function(file) {
  first_line <- readLines(file, n = 1)
  skip_lines <- if (grepl("^N_unmapped", first_line)) 4 else 0
  d <- read.table(file, header = FALSE, skip = skip_lines, stringsAsFactors = FALSE)
  if (!is.null(opts$s)){
    sample_name <- opts$s
  }else{
    sample_name <- basename(file)
  }
  df <- setNames(d[, 1:2], c("GeneID", sample_name))
  
  return(df)
}

# Merge counts
verbose_print("Reading and merging STAR count files...", verbose = opts$verbose)
count_list <- lapply(opts$c, read_star)
merged_counts <- Reduce(function(a, b) merge(a, b, by = "GeneID", all = TRUE), count_list)
merged_counts[is.na(merged_counts)] <- 0
sample_names <- colnames(merged_counts)[-1]
genecol <- "GeneID"
mat_input <- merged_counts[, c(genecol, sample_names)]
colnames(mat_input)[1] <- "X"
mat_input$temporarysample <- mat_input[[2]]


# Run deconvolution and resistance analysis
verbose_print("Running seAMLess_modified...", verbose = opts$verbose)
# run the analysis
res <- seAMLess_modified(mat = mat_input, scRef = ref_sub, verbose = opts$verbose)

# pull out the two pieces
deconv   <- res$Deconvolution
veno_res <- res$Venetoclax.resistance

# drop the dummy sample
if ("temporarysample" %in% rownames(deconv)) {
  deconv <- deconv[rownames(deconv) != "temporarysample", , drop = FALSE]
}

# if veno_res is a named vector, drop it there too:
if (! is.null(names(veno_res)) && "temporarysample" %in% names(veno_res)) {
  veno_res <- veno_res[names(veno_res) != "temporarysample"]
}
veno_res_df <- data.frame(Venetoclax.resistance = veno_res)
rownames(veno_res_df) <- names(veno_res)

out_deconv <- file.path(opts$o, "deconvolution.csv")
out_veno   <- file.path(opts$o, "venetoclax_resistance.csv")
# write them out
write.csv(deconv,  out_deconv, row.names = TRUE)
write.csv(veno_res, out_veno,   row.names = TRUE)

# Create a plot of the cell types
png(file.path(opts$o, "cell-types.png"))

tibble::enframe(deconv[1,]) |>
ggplot(aes(name, value *100, fill = name)) +
geom_bar(stat = "identity") + ylab("Cell type composition (%)") +
guides(x = guide_axis(angle = 90)) + theme_minimal() + labs(fill = "Cell Types") + theme_bw(base_size = 16) + theme(legend.position = "none") + xlab("")
dev.off()
