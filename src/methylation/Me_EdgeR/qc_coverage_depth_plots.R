#########################################################
## Required libraries
library(dplyr)

#########################################################
# Helper: load DGEList if input is a path
.load_dge <- function(me_dge) {
  if (is.character(me_dge) && length(me_dge) == 1L) {
    if (!file.exists(me_dge)) stop("RDS file not found at: ", me_dge)
    readRDS(me_dge)
  } else {
    me_dge
  }
}

#########################################################
# Helper: biological-sample coverage matrix (Me+Un per bio sample)
# Returns a matrix: features x n_bio_samples; colnames are base sample IDs.
.bio_cov <- function(dge) {

  # Case A: separate unmethylated slot => columns already correspond to biological samples
  if (!is.null(dge$unmethylated)) {
    return(dge$counts + dge$unmethylated)
  }

  # Case B: single matrix with -Me / -Un columns
  X <- dge$counts
  cn <- colnames(X)

  me_cols <- grep("-Me$", cn, value = TRUE)
  un_cols <- grep("-Un$", cn, value = TRUE)
  if (length(me_cols) == 0 || length(un_cols) == 0) {
    stop("Could not find both '-Me' and '-Un' columns (and no dge$unmethylated slot present).")
  }

  bases <- sort(unique(c(sub("-Me$", "", me_cols), sub("-Un$", "", un_cols))))
  me_map <- setNames(me_cols, sub("-Me$", "", me_cols))
  un_map <- setNames(un_cols, sub("-Un$", "", un_cols))

  cov_bio <- matrix(0,
                    nrow = nrow(X),
                    ncol = length(bases),
                    dimnames = list(rownames(X), bases))

  for (b in bases) {
    if (!is.na(me_map[b])) cov_bio[, b] <- cov_bio[, b] + X[, me_map[b]]
    if (!is.na(un_map[b])) cov_bio[, b] <- cov_bio[, b] + X[, un_map[b]]
  }

  cov_bio
}

#########################################################
# Helper: biological-sample group mapping aligned to .bio_cov() columns
# Returns a named character vector: names = bio sample IDs, values = Group.
.bio_group <- function(dge, bio_samples) {

  if (is.null(dge$samples) || !("Group" %in% colnames(dge$samples))) {
    stop("dge$samples must contain a 'Group' column.")
  }

  smp <- as.data.frame(dge$samples)
  smp$colname <- rownames(smp)

  # Prefer explicit Sample column (base IDs) if present
  if ("Sample" %in% colnames(smp)) {
    smp$base <- as.character(smp$Sample)
  } else {
    smp$base <- sub("-(Me|Un)$", "", smp$colname)
  }

  smp$Group <- as.character(smp$Group)

  grp_by_base <- dplyr::distinct(smp, base, Group)
  if (any(duplicated(grp_by_base$base))) {
    stop("A biological sample maps to multiple Group values in dge$samples (check metadata).")
  }

  gmap <- setNames(grp_by_base$Group, grp_by_base$base)
  gmap <- gmap[bio_samples]

  if (any(is.na(gmap))) {
    stop("Some biological samples have no Group mapping in dge$samples. Check naming consistency.")
  }

  gmap
}

#########################################################
# Helper: indices for a strand ("+" or "-")
.strand_idx <- function(dge, strand = c("+", "-", "merged")) {
  strand <- match.arg(strand)
  if (is.null(dge$genes) || !("Strand" %in% colnames(dge$genes))) {
    stop("dge$genes must contain a 'Strand' column.")
  }

  if (strand == "merged") {
    return(seq_len(nrow(dge$counts)))  # "consider all the rows"
  }

  which(as.character(dge$genes$Strand) == strand)
}

#########################################################
# Helper: panel grid for mosaic plotting
.panel_grid <- function(n) {
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)
  c(nr, nc)
}

#########################################################
# Helper: depth histogram bin edges (reads per CpG)
# Fine bins at low depth + coarser bins at higher depth.
.depth_breaks <- function(max_break = 200) {
  b1 <- 0:10
  b2 <- seq(15, 50, by = 5)
  b3 <- seq(60, max_break, by = 10)
  sort(unique(c(b1, b2, b3, max_break + 1)))
}

#########################################################
# Helper: histogram with log10 y-axis (counts + 1 to avoid log(0))
.hist_logy <- function(x, breaks, xlab, main) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    plot.new(); title(main)
    return(invisible(NULL))
  }

  # Ensure breaks span the full range of x (IMPORTANT)
  breaks <- sort(unique(breaks))
  if (max(x) >= max(breaks)) {
    breaks <- c(breaks, max(x) + 1)
  }

  h <- hist(x, breaks = breaks, plot = FALSE, right = FALSE)
  plot(h$mids, h$counts + 1,
       log = "y", type = "h", lwd = 2,
       xlab = xlab, ylab = "Number of CpGs (log10)",
       main = main)
  points(h$mids, h$counts + 1, pch = 16, cex = 0.5)
}

#########################################################
# Depth: per-sample mosaic (one panel per biological sample), strand-specific
# X-axis: log10(reads per CpG + 1), Y-axis: log10(counts)
.plot_depth_mosaic_per_sample <- function(dge, pdf_path, strand = c("+", "-", "merged"), max_break = 200) {
  strand <- match.arg(strand)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

  cov_bio <- .bio_cov(dge)                      # features x bio samples
  idx     <- .strand_idx(dge, strand)
  cov_s   <- cov_bio[idx, , drop = FALSE]

  smp      <- colnames(cov_s)
  grid     <- .panel_grid(length(smp))
  brks_raw <- .depth_breaks(max_break)          # breaks in depth space

  pdf(pdf_path, width = 3.5 * grid[2], height = 3.0 * grid[1])
  op <- par(mfrow = c(grid[1], grid[2]), mar = c(4, 4, 3, 1) + 0.1)

  for (s in smp) {
    depths   <- cov_s[, s]
    x_log    <- log10(depths + 1)
    brks_log <- log10(brks_raw + 1)

    .hist_logy(
      x      = x_log,
      breaks = brks_log,
      xlab   = "log10(reads per CpG + 1)",
      main   = paste0(s, " (", strand, ")")
    )
  }

  par(op); dev.off()
  invisible(pdf_path)
}

#########################################################
# Depth: merged distribution for selected samples (or all), strand-specific
# X-axis: log10(reads per CpG + 1), Y-axis: log10(counts)
.plot_depth_merged <- function(dge, pdf_path, strand = c("+", "-", "merged"),
                               samples_keep = NULL, label = "", max_break = 200) {
  strand <- match.arg(strand)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

  cov_bio <- .bio_cov(dge)
  idx     <- .strand_idx(dge, strand)
  cov_s   <- cov_bio[idx, , drop = FALSE]

  if (!is.null(samples_keep)) {
    cov_s <- cov_s[, samples_keep, drop = FALSE]
  }

  depths   <- as.vector(cov_s)                  # pool across features and samples
  brks_raw <- .depth_breaks(max_break)

  pdf(pdf_path, width = 7, height = 5)
  x_log    <- log10(depths + 1)
  brks_log <- log10(brks_raw + 1)

  .hist_logy(
    x      = x_log,
    breaks = brks_log,
    xlab   = "log10(reads per CpG + 1)",
    main   = paste0("Depth distribution ", label, " (", strand, ")")
  )
  dev.off()

  invisible(pdf_path)
}

#########################################################
# Depth: merged distributions by group, in one PDF (one panel per group)
.plot_depth_merged_by_group <- function(dge, pdf_path, strand = c("+", "-", "merged"),
                                        groups, gmap, max_break = 200) {
  strand <- match.arg(strand)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

  cov_bio <- .bio_cov(dge)
  idx     <- .strand_idx(dge, strand)
  cov_s   <- cov_bio[idx, , drop = FALSE]

  brks_raw <- .depth_breaks(max_break)

  pdf(pdf_path, width = 10, height = 5)
  op <- par(mfrow = c(1, length(groups)), mar = c(5, 4, 4, 1) + 0.1)

  for (g in groups) {
    keep <- names(gmap)[gmap == g]

    if (length(keep) == 0) {
      plot.new(); title(paste0(g, " (no samples)"))
      next
    }

    cov_g    <- cov_s[, keep, drop = FALSE]
    depths   <- as.vector(cov_g)
    x_log    <- log10(depths + 1)
    brks_log <- log10(brks_raw + 1)

    .hist_logy(
      x      = x_log,
      breaks = brks_log,
      xlab   = paste0("log10(reads per CpG + 1) in ", g),
      main   = paste0("Depth distribution (", g, ", ", strand, ")")
    )
  }

  par(op); dev.off()
  invisible(pdf_path)
}

#########################################################
# Shared CpG: merged distributions by group, in one PDF (one panel per group)
.plot_shared_merged_by_group <- function(dge, pdf_path, strand = c("+", "-", "merged"),
                                         groups, gmap) {
  strand <- match.arg(strand)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

  cov_bio <- .bio_cov(dge)
  idx     <- .strand_idx(dge, strand)
  cov_s   <- cov_bio[idx, , drop = FALSE]

  pdf(pdf_path, width = 10, height = 5)
  op <- par(mfrow = c(1, length(groups)), mar = c(5, 4, 4, 1) + 0.1)

  for (g in groups) {
    keep <- names(gmap)[gmap == g]

    if (length(keep) == 0) {
      plot.new(); title(paste0(g, " (no samples)"))
      next
    }

    cov_g <- cov_s[, keep, drop = FALSE]
    k     <- rowSums(cov_g > 0)                 # 0..nSamples_in_group
    brks  <- seq(-0.5, ncol(cov_g) + 0.5, by = 1)

    .hist_logy(
      x      = k,
      breaks = brks,
      xlab   = paste0("# bio samples with ≥1 read (", g, ")"),
      main   = paste0("Shared CpG distribution (", g, ", ", strand, ")")
    )
  }

  par(op); dev.off()
  invisible(pdf_path)
}

#########################################################
# Shared CpG: merged distribution across CpGs (k = #samples with >=1 read), strand-specific
# If samples_keep is provided, k is computed within that subset (group-specific sharing).
.plot_shared_merged <- function(dge, pdf_path, strand = c("+", "-", "merged"),
                                samples_keep = NULL, label = "") {
  strand <- match.arg(strand)
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

  cov_bio <- .bio_cov(dge)
  idx <- .strand_idx(dge, strand)
  cov_s <- cov_bio[idx, , drop = FALSE]

  if (!is.null(samples_keep)) {
    cov_s <- cov_s[, samples_keep, drop = FALSE]
  }

  k <- rowSums(cov_s > 0)  # 0..nSamplesKeep
  brks <- seq(-0.5, ncol(cov_s) + 0.5, by = 1)

  pdf(pdf_path, width = 7, height = 5)
  .hist_logy(
    x = k,
    breaks = brks,
    xlab = "Number of biological samples with ≥1 read at CpG",
    main = paste0("Shared CpG distribution ", label, " (", strand, ")")
  )
  dev.off()

  invisible(pdf_path)
}

#########################################################
#' Depth & Shared-CpG QC (strand-split): mosaics, per-group PDFs, and merged PDFs
#'
#' For EACH strand (+ and -), writes:
#' Depth:
#'  - per-sample mosaic (one panel per biological sample)
#'  - merged across all samples
#'  - merged per group (one PDF per group)
#'
#' Shared CpG:
#'  - per-sample mosaic (one panel per biological sample)
#'  - merged across all samples (k across all samples)
#'  - merged per group (k within group; one PDF per group)
#'
#' @param me_dge DGEList object OR path to an .rds containing the DGEList.
#' @param output_dir Output directory for PDFs.
#' @param groups Character vector of group names (default c("Control","MSUS")).
#' @param max_break Depth histogram max break (controls binning; default 200).
#' @return Invisibly returns a named list of output file paths.
#' @export
analyze_cpg_depth_qc <- function(me_dge,
                                 output_dir,
                                 groups = c("Control", "MSUS"),
                                 max_break = 200) {

  dge <- .load_dge(me_dge)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (is.null(dge$samples) || !("Group" %in% colnames(dge$samples))) stop("dge$samples must contain a 'Group' column.")
  if (is.null(dge$genes)   || !("Strand" %in% colnames(dge$genes)))   stop("dge$genes must contain a 'Strand' column.")

  # only need bio sample names + group mapping; do NOT keep cov_bio
  bio_samples <- if (!is.null(dge$unmethylated)) {
    colnames(dge$counts)
  } else {
    sub("-Me$", "", grep("-Me$", colnames(dge$counts), value = TRUE))
  }
  bio_samples <- sort(unique(bio_samples))
  gmap <- .bio_group(dge, bio_samples)

  out <- list()

  strand_char <- as.character(dge$genes$Strand)
  has_merged <- any(strand_char == "merged", na.rm = TRUE)

  strand_set <- if (has_merged) "merged" else c("+", "-")

  for (strand in strand_set) {

    strand_tag <- switch(
      strand,
      "+"      = "PStrand",
      "-"      = "NStrand",
      "merged" = "MStrand"
    )

    out[[paste0("depth_mosaic_", strand)]] <-
      file.path(output_dir, paste0("Depth_Mosaic_PerSample_", strand_tag, ".pdf"))
    .plot_depth_mosaic_per_sample(dge, out[[paste0("depth_mosaic_", strand)]],
                                  strand = strand, max_break = max_break)

    gc()

    out[[paste0("depth_all_", strand)]] <-
      file.path(output_dir, paste0("Depth_Distribution_AllSamples_", strand_tag, ".pdf"))
    .plot_depth_merged(dge, out[[paste0("depth_all_", strand)]],
                       strand = strand, samples_keep = NULL,
                       label = "(all samples)", max_break = max_break)

    gc()

    out[[paste0("depth_by_group_", strand)]] <-
      file.path(output_dir, paste0("Depth_Distribution_ByGroup_", strand_tag, ".pdf"))
    .plot_depth_merged_by_group(
      dge, pdf_path = out[[paste0("depth_by_group_", strand)]],
      strand = strand, groups = groups, gmap = gmap, max_break = max_break
    )

    gc()

    out[[paste0("shared_all_", strand)]] <-
      file.path(output_dir, paste0("SharedCpG_Distribution_AllSamples_", strand_tag, ".pdf"))
    .plot_shared_merged(dge, out[[paste0("shared_all_", strand)]],
                        strand = strand, samples_keep = NULL,
                        label = "(all samples)")

    gc()

    out[[paste0("shared_by_group_", strand)]] <-
      file.path(output_dir, paste0("SharedCpG_Distribution_ByGroup_", strand_tag, ".pdf"))
    .plot_shared_merged_by_group(
      dge, pdf_path = out[[paste0("shared_by_group_", strand)]],
      strand = strand, groups = groups, gmap = gmap
    )

    gc()
  }

  message("Depth/shared-CpG QC PDFs saved to: ", output_dir)
  invisible(out)
}
