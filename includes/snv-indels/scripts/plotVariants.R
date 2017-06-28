#!/usr/bin/env Rscript
#
# Script for plotting variant calls from an input JSON payload.
#
# Copyright (c) 2015 Leiden University Medical Center <http://lumc.nl>
# All rights reserved.

suppressMessages(library(Gviz))

RenameEnsemblChroms <- function(annots) {
  # Renames Ensembl chromosomes in the given GenomicRanges object to use the GRCh38 naming.
  #
  # This means:
  #   - Autosomal chromosome names are prepended with `chr`
  #   - Sex chromosome names are prepended with `chr`
  #   - Mitchondrial genome is renamed from `MT` to `chrM`
  #   - All other chromosomes are left untouched
  #
  # Input:
  #   GenomicRanges object
  #
  # Returns:
  #   Character vector of the renamed chromosomes
  canonicalChroms <- c(as.character(1:22), c("X", "Y"))
  renameFunc <- function(n) {
    if (n %in% canonicalChroms) paste("chr", n, sep="")
    else if (n == "MT") "chrM"
    else n
  }
  unlist(lapply(GenomeInfoDb::seqlevels(annots), renameFunc))
}

ParseJson <- function(input.json) {
  # Parses the input JSON file containing variant data to plot.
  #
  # Args:
  #   input.json: Path to input JSON file.
  #
  # Returns:
  #   A list containing the sample IDs, gene IDs, and the parsed variants data.
  parsed <- jsonlite::fromJSON(input.json)$samples
  samples <- names(parsed)
  gene.ids <- unique(as.vector(sapply(parsed, function(o) { names(o) })))
  # for each sample
  parsed <- lapply(parsed, function(perSample) {
    # for each gene within a sample
    lapply(perSample, function(perGene) {
      varscan <- perGene$varscan
      if (is.null(varscan)) {
        NULL
      } else {
        # create list of zipped (vep, varscan) values
        mapply(function(a, b) list(varscan = a, vep = dplyr::tbl_df(b)),
               as.data.frame(as.list(t(varscan))), perGene$vep, SIMPLIFY = FALSE)
      }
    })
  })
  list(sample.ids = samples, gene.ids = gene.ids, data = parsed)
}

ParseIdMappings <- function(input.ids) {
  # Parses the id mappings file.
  #
  # ID mappings files consists of three columns: ENSG ID, gene symbol, and
  # ENST ID (in that order).
  #
  # Args:
  #   input.ids: Path to IDs file.
  #
  # Returns:
  #   The IDs stored in a vector of characters.
  idm <- read.table(input.ids, header = T, stringsAsFactors = F)
  # split toi_ids entries each into a vector
  dplyr::mutate(idm, TOI_IDS = strsplit(TOI_IDS, ","))
}

ParseGeneModel <- function(input.gtf, id.map=NULL) {
  # Parses the GTF file containing the gene models.
  #
  # Note: the gene IDs in this GTF file must be of the same type with the IDs
  # in the JSON and IDs file.
  #
  # Args:
  #   input.gtf: Path to input GTF file.
  #   id.map: data frame containing ID mappings to select
  #
  # Returns:
  #   GenomicRanges object of the exons.
  gtf <- rtracklayer::import(input.gtf)
  gtf.goi <-
    if (is.null(id.map)) gtf
    else gtf[gtf$gene_id %in% id.map[["GOI_ID"]]]
  GenomeInfoDb::seqlevels(gtf.goi) <- RenameEnsemblChroms(gtf.goi)
  # add columns required for Gviz plotting and select only exons
  gtf.goi$transcript <- gtf.goi$transcript_id
  gtf.goi$gene <- gtf.goi$gene_id
  gtf.goi$exon <- gtf.goi$exon_id
  gtf.goi$feature <- gtf.goi$type
  gtf.goi$symbol <- gtf.goi$gene_name
  # select only CDS and UTR so we can show differences between them
  gtf.goi[gtf.goi$feature %in% c("UTR", "CDS"),]
}

MakePlot <- function(parsed.data, sample.name, gene.id, gene.model, genome = "hg38",
                     amp.bed = NULL, hot.bed = NULL, out.file = NULL) {

  var.data <- parsed.data$data[[sample.name]][[gene.id]]

  if (!is.null(var.data)) {

    chrom <- var.data[[1]]$varscan$CHROM

    MakeVariantRanges <- function(chrom) {
      width <- 10
      poss <- unlist(lapply(var.data, function(n) { n$varscan$POS }))
      ranges <- GenomicRanges::GRanges(seqnames = chrom,
              ranges = IRanges::IRanges(start = poss - width, end = poss + width),
              strand = rep("*", length(poss)))
      ranges$transcript_ids <- lapply(var.data, function(n) { n$vep$Feature })
      ranges$filtered_in <- unlist(lapply(var.data, function(n) { !("SubpopAFAtLeast5Pct" %in% n$varscan$filters) & !("LowQualBases" %in% n$varscan$filters) }))
      ranges$has_high_impact <- unlist(lapply(var.data, function(n) { "HIGH" %in% n$vep$impact }))
      ranges$has_other_impact <- unlist(lapply(var.data, function(n) { !("HIGH" %in% n$vep$impact) | length(n$vep$impact[n$vep$impact != "HIGH"]) > 1 }))
      ranges$internal_variant_id <- apply(cbind(chrom, as.character(poss)), 1, function(n) { paste(n, collapse="_") })
      ranges$affects_toi <- unlist(lapply(var.data, function(n) { !is.null(n$vep$is_toi) && as.logical(any(na.omit(n$vep$is_toi))) }))
      ranges
    }

    MakePerPosAnnotation <- function(entry) {

      varscan <- entry$varscan
      vep <- entry$vep
      high.impacts <- dplyr::filter(vep, impact == "HIGH")
      high.impacts <- dplyr::filter(vep, !is.null(is_toi))
      high.impacts <- dplyr::filter(vep, is_toi)
      tids <-  apply(high.impacts, 1, function(n) { n$Feature })
      len <- nrow(high.impacts)
      fallback <- "n/a"

      Append <- function(cont, item) {
        is.fallback <- item[2:length(item)] == fallback
        if (!all(is.fallback)) {
          cont[[length(cont) + 1]] <- item
          cont
        } else { cont }
      }
      GetOrElse <- function(obj, attr.name, fun = NULL) {
        v <- obj[[attr.name]]
        if (is.null(v)) fallback
        else if (!is.null(fun)) fun(v)
        else v
      }
      CatLine <- function(n) { paste(n, collapse="\n") }
      GetHgvs <- function(obj, attr.name, id.attr.name, replacement) {
        v <- obj[[attr.name]]
        if (is.null(v)) fallback
        else {
          hgvs.ids <- GetOrElse(obj, id.attr.name)
          as.vector(unlist(mapply(function(eid, hgvs) gsub(eid, replacement, hgvs),
                                  hgvs.ids, v, SIMPLIFY = FALSE)))
        }
      }

      l <- list()
      l <- Append(l, c("Genotype", varscan$genotype, rep("", len - 1)))
      l <- Append(l, c("P-value", format(varscan$PVAL, scientific = TRUE), rep("", len - 1)))
      l <- Append(l, c("Depths (Ref | Alt | Total)", paste(c(varscan$RD[[1]], varscan$AD[[1]], varscan$DP[[1]]), collapse=" | "), rep("", len - 1)))
      l <- Append(l, c("BaseQ (Ref | Alt)", paste(c(varscan$RBQ[[1]], varscan$ABQ[[1]]), collapse=" | "), rep("", len - 1)))
      l <- Append(l, c("Position", prettyNum(varscan$POS, big.mark = ",", decimal.mark = "."), rep("", len - 1)))
      l <- Append(l, c("In Hotspot", GetOrElse(varscan, "is_in_hotspot", function(n) if (n) "yes" else "no"), rep("", len - 1)))

      l <- Append(l, c("Consequence", apply(high.impacts, 1, function(r) { CatLine(r$Consequence) })))
      l <- Append(l, c("Known As",  apply(high.impacts, 1, function(r) { GetOrElse(r, "Existing_variation", CatLine) })))
      l <- Append(l, c("Affected Protein", apply(high.impacts, 1, function(r) { CatLine(r$ENSP) })))
      l <- Append(l, c("HGVS - cDNA",  apply(high.impacts, 1, function(r) { GetHgvs(r, "HGVSc", "Feature", "{enst}") })))
      l <- Append(l, c("HGVS - Protein",  apply(high.impacts, 1, function(r) { GetHgvs(r, "HGVSp", "ENSP", "{ensp}") })))
      l <- Append(l, c("SIFT score",  apply(high.impacts, 1, function(r) { GetOrElse(r, "SIFT") })))
      l <- Append(l, c("PolyPhen",  apply(high.impacts, 1, function(r) { GetOrElse(r, "PolyPhen") })))

      if (length(varscan[["GONL_AF"]]) > 0) {
        l <- Append(l, c("GoNL_AF",  apply(high.impacts, 1, function(r) { paste(varscan[["GONL_AF"]], collapse=", ") })))
      }
      for (af1kg in c("1KG_P3_AF", "1KG_P3_AFR_AF", "1KG_P3_AMR_AF", "1KG_P3_EAS_AF", "1KG_P3_EUR_AF", "1KG_P3_SAS_AF")) {
        if (length(varscan[[af1kg]]) > 0) {
          rname <- paste(c("1KG Phase 3 - ", gsub("1KG_P3_", "", af1kg)), collapse="")
          l <- Append(l, c(rname, apply(high.impacts, 1, function(r) { paste(varscan[[af1kg]], collapse=", ") })))
        }
      }

      dd <- as.data.frame(setNames(replicate(len+1,numeric(0), simplify = F), c("param", tids)), stringsAsFactors=F)
      dd <- t(rbind(dd, as.data.frame(l, stringsAsFactors=F)))
      rownames(dd) <- dd[,1]
      annots <- as.data.frame(dd[,2:ncol(dd)])
      colnames(annots) <- tids
      list(table = annots, num.transcripts = length(tids))
    }

    SelectFunc <- function(identifier, ...) {
      toks <- unlist(strsplit(identifier, "_"))
      chrom <- toks[1]
      pos <- as.numeric(toks[2])
      entry <- Filter(function(i) { i$varscan$CHROM == chrom && i$varscan$POS == pos }, var.data)[[1]]
      if (is.null(entry$vep$is_toi)) {
        FALSE
      } else {
        ("HIGH" %in% entry$vep$impact) & as.logical(any(na.omit(entry$vep$is_toi)))
      }
    }

    my.theme <- gridExtra::ttheme_default(
      core = list(fg_params=list(cex = 1.2)),
      colhead = list(fg_params=list(cex = 1.0)),
      rowhead = list(fg_params=list(cex = 1.0)))

    MergeCols <- function(grob, row.idx, ntid) {
      to.display <- which(grob$layout$t == row.idx & grob$layout$l == 2)
      grob$layout[to.display, c("l")] <- list(2)
      grob$layout[to.display, c("r")] <- list(ntid)
      to.merge <- which(grob$layout$t == row.idx & grob$layout$l != 1 & grob$layout$l != 2)
      grob$layout[to.merge, c("l")] <- list(2)
      grob$layout[to.merge, c("r")] <- list(ntid)
      grob
    }

    PlotTableFunc <- function(identifier, ...) {
      toks <- unlist(strsplit(identifier, "_"))
      chrom <- toks[1]
      pos <- as.numeric(toks[2])
      entry <- Filter(function(i) { i$varscan$CHROM == chrom && i$varscan$POS == pos }, var.data)[[1]]
      pos.annot <- MakePerPosAnnotation(entry)
      annot.table <- pos.annot$table
      ntids <- pos.annot$num.transcripts
      to.draw <- if (ntids > 1) {
        d <- gridExtra::tableGrob(annot.table, theme = my.theme,
                                  cols = colnames(annot.table),
                                  rows = rownames(annot.table))
        rownum.to.merge <- rle(as.vector(annot.table[,2]))$lengths[1] + 1
        for (row.idx in 2:rownum.to.merge) { d <- MergeCols(d, row.idx, ntids + 1) }
        d
      } else {
        gridExtra::tableGrob(annot.table, theme = my.theme,
                             cols = colnames(annot.table),
                             rows = rownames(annot.table))
      }
      print(grid::grid.draw(to.draw), newpage = FALSE, prefix = "table")
    }

    variants <- MakeVariantRanges(chrom)

    # variants that pass the AF threshold and affect any transcripts of interest
    variants.passfilter <- variants[variants$filtered_in,]

    # variants we want to display in detail
    variants.detailed <- variants.passfilter[variants.passfilter$affects_toi,]
    variants.detailed <- variants.passfilter[variants.passfilter$has_high_impact,]

    itrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chrom, fontsize = 19)
    gtrack <- Gviz::GenomeAxisTrack(fontsize = 19)
    grtrack <- Gviz::GeneRegionTrack(gene.model[gene.model$gene_id == gene.id, ],
                                      chromosome = chrom, genome = genome,
                                      transcriptAnnotation = "transcript",
                                      thinBoxFeature = c("UTR"), name = "Transcripts",
                                      background.title = "seagreen",
                                      fill = "seagreen", stack.height = 0.5,
                                      fontsize = 18)

    valltrack <- Gviz::AnnotationTrack(range = variants, genome = genome,
                                       name = "All",
                                       id = variants$internal_variant_id,
                                       stacking = "dense",
                                       fill = "slategrey", background.title = "slategrey",
                                       fontsize = 18)

    vpasstrack <- Gviz::AnnotationTrack(range = variants.passfilter, genome = genome,
                                        name = "Pass Filter",
                                        id = variants.passfilter$internal_variant_id,
                                        stacking = "dense",
                                        fill = "dodgerblue3", background.title = "dodgerblue3",
                                        fontsize = 18)

    vhitrack <- Gviz::AnnotationTrack(range = variants.detailed, genome = genome,
                                      name = "Pass Filter + High Impact on Transcript of Interest",
                                      id = variants.detailed$internal_variant_id,
                                      fun = PlotTableFunc, selectFun = SelectFunc,
                                      details.size = 0.75, stacking = "dense",
                                      detailsConnector.col = "firebrick1", detailsBorder.col = "firebrick1",
                                      fill = "firebrick1", background.title = "firebrick1",
                                      fontsize = 18)

    tmp.plot <-

      if (length(variants) == 0) {
        list(tracks = c(itrack, gtrack, grtrack),
             sizes = c(0.5, 0.9, 3),
             height = 500, width = 800)

      } else if (length(variants.passfilter) == 0) {
        list(tracks = c(itrack, valltrack, gtrack, grtrack),
             sizes = c(0.5, 1.2, 0.9, 3),
             height = 500, width = 1000)

      } else if (length(variants.detailed) == 0) {
        list(tracks = c(itrack, vpasstrack, valltrack, gtrack, grtrack),
             sizes = c(0.5, 1.2, 1.3, 0.9, 3),
             height = 500, width = 1200)

      } else {
        list(tracks = c(itrack, vhitrack, vpasstrack, valltrack, gtrack, grtrack),
            sizes = c(0.5, 7.1, 1.5, 1.5, 0.9, 3),
            height = 1500, width = 1000)
      }

    gene.symbol <- gene.model[gene.model$gene_id == gene.id, ]$gene_name[1]
    hot.bed.intersects <-
      if (!is.null(hot.bed)) {
        GenomicRanges::intersect(hot.bed, gene.model[gene.model$gene_id == gene.id, ], ignore.strand=T)
      } else {
        list()
      }

    tmp2.plot <-

      if (length(hot.bed.intersects) > 0) {
        hottrack <- Gviz::AnnotationTrack(hot.bed.intersects,
                                          name = "Hotspots", stacking = "squish",
                                          fill = "sienna1", col = "white",
                                          background.title="sienna1")
        tracks <- c(head(tmp.plot$tracks, n=-1), hottrack, tail(tmp.plot$tracks, n=1))
        sizes <- c(head(tmp.plot$sizes, n=-1), 0.8, tail(tmp.plot$sizes, n=1))
        list(tracks = as.list(tracks), sizes = sizes,
             height = tmp.plot$height, width = tmp.plot$width)
      } else {
        list(tracks = as.list(tmp.plot$tracks), sizes = tmp.plot$sizes,
             height = tmp.plot$height, width = tmp.plot$width)
      }

    amp.bed.intersects <-
      if (!is.null(amp.bed)) {
        GenomicRanges::intersect(amp.bed, gene.model[gene.model$gene_id == gene.id, ], ignore.strand=T)
      } else {
        list()
      }

    to.plot <-

      if (length(amp.bed.intersects) > 0) {
        amptrack <- Gviz::AnnotationTrack(amp.bed.intersects,
                                          name = "Amplicon", stacking = "squish",
                                          fill = "slateblue", col = "white",
                                          background.title="slateblue")
        tracks <- c(head(tmp2.plot$tracks, n=-1), amptrack, tail(tmp2.plot$tracks, n=1))
        sizes <- c(head(tmp2.plot$sizes, n=-1), 0.8, tail(tmp2.plot$sizes, n=1))
        list(tracks = as.list(tracks), sizes = sizes,
             height = tmp2.plot$height, width = tmp2.plot$width)
      } else {
        list(tracks = as.list(tmp2.plot$tracks), sizes = tmp2.plot$sizes,
             height = tmp2.plot$height, width = tmp2.plot$width)
      }

    plot.title <- paste(c("Variant Calls - Sample '", sample.name, "'\n", gene.symbol,
                          " (", gene.id, ")"), collapse = "")

    if (!is.null(out.file) & tolower(unlist(strsplit(out.file, "\\."))[-1]) == "png") {
      png(out.file, height = to.plot$height, width = to.plot$width)
      Gviz::plotTracks(to.plot$tracks, chromosome = chrom, main = plot.title,
                col = NULL, sizes=to.plot$sizes)
      dev.off()
    } else {
      Gviz::plotTracks(to.plot$tracks, chromosome = chrom, main = plot.title,
                col = NULL, sizes=to.plot$sizes)
    }
  } else {
    NULL
  }
}

# Main script loop
if (!interactive()) {
  spec <- matrix(c(
    # Input JSON file
    "input-json", "j", 1, "character",
    # Input Ensembl IDs to visualize
    "input-gene-ids", "i", 1, "character",
    # Input GTF file containing gene models
    "input-gtf", "m", 1, "character",
    # Input BED file containing amplicon regions and their primers
    "input-amplicon-bed", "a", 1, "character",
    # Input BED file containing hotspot regions
    "input-hotspot-bed", "t", 1, "character",
    # Output directory
    "out-dir", "o", 1, "character",
    # Help
    "help", "h", 0, "logical"
  ), byrow = TRUE, ncol = 4)
  opt <- getopt::getopt(spec)

  if (!is.null(opt[["help"]])) {
      cat(getopt::getopt(spec, usage=TRUE))
    q(status=1)
  }

  if (is.null(opt[["input-json"]])) {
      message("Error: '--input-json/-j' must have a value.")
    q(status=2)
  }
  INPUT.JSON <- opt[["input-json"]]

  if (is.null(opt[["input-gene-ids"]])) {
      message("Error: '--input-gene-ids/-i' must have a value.")
    q(status=3)
  }
  INPUT.IDS <- opt[["input-gene-ids"]]

  if (is.null(opt[["input-gtf"]])) {
      message("Error: '--input-gtf/-m' must have a value.")
    q(status=4)
  }
  INPUT.GTF <- opt[["input-gtf"]]

  INPUT.AMPLICON.BED <- opt[["input-amplicon-bed"]]
  INPUT.HOTSPOT.BED <- opt[["input-hotspot-bed"]]

  msg <- function(...) { cat(sprintf(...), sep="", file=stderr()) }

  # set fallback values for optional args
  BASE.OUT.DIR <- if (!is.null(opt[["out-dir"]])) {
    out.dir <- opt[["out-dir"]]
    # create directory if it doesn't exist
    msg("Creating output directory '%s' ...\n", out.dir)
    dir.create(out.dir, showWarnings = FALSE)
    normalizePath(out.dir)
  } else {
    getwd()
  }

  msg("Reading JSON file ...\n")
  dat <- ParseJson(INPUT.JSON)
  msg("Reading IDs file ...\n")
  id.map <- ParseIdMappings(INPUT.IDS)
  msg("Reading GTF file ...\n")
  gtf <- ParseGeneModel(INPUT.GTF, id.map)
  amp.bed <-
    if (!is.null(INPUT.AMPLICON.BED)) {
      rtracklayer::import.bed(INPUT.AMPLICON.BED)
    } else { NULL }
  hot.bed <-
    if (!is.null(INPUT.HOTSPOT.BED)) {
      rtracklayer::import.bed(INPUT.HOTSPOT.BED)
    } else { NULL }

  if (length(dat$sample.ids) != 1) {
    message("JSON data contains less or more than one sample.")
    q(status=5)
  }

  sample.id <- dat$sample.ids[1]
  out.dir <- file.path(BASE.OUT.DIR, "variant_plots")
  dir.create(out.dir, showWarnings = FALSE)
  for (gene.id in names(dat$data[[sample.id]])) {
    gene.name <- unique(dat$data[[sample.id]][[gene.id]][[1]]$vep$SYMBOL)
    base.name <- paste(c(paste(c("sample", sample.id, "gene", gene.name), sep="_", collapse="_"),
                          ".png"), sep="", collapse="")
    out.file <- file.path(out.dir, base.name)
    msg("Plotting sample '%s', gene '%s' to '%s' ...\n", sample.id, gene.name, out.file)
    MakePlot(dat, sample.id, gene.id, gtf, out.file = out.file, amp.bed = amp.bed, hot.bed = hot.bed)
  }
}
