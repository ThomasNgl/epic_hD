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
# Helper: panel grid for mosaic plotting
.panel_grid <- function(n) {
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)
  c(nr, nc)
}

#########################################################
# Helper: histogram with log10 y-axis (counts + 1 to avoid log(0))
.hist_logy <- function(x, breaks, xlab, main, ylab = "Count (log10)") {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    plot.new(); title(main)
    return(invisible(NULL))
  }

  breaks <- sort(unique(breaks))
  if (max(x) >= max(breaks)) {
    breaks <- c(breaks, max(x) + 1)
  }

  h <- hist(x, breaks = breaks, plot = FALSE, right = FALSE)
  plot(h$mids, h$counts + 1,
       log = "y", type = "h", lwd = 2,
       xlab = xlab, ylab = ylab,
       main = main)
  points(h$mids, h$counts + 1, pch = 16, cex = 0.5)
}

#########################################################
# Helper: depth histogram bin edges (reads per region)
.depth_breaks <- function(max_break = 200) {
  b1 <- 0:10
  b2 <- seq(15, 50, by = 5)
  b3 <- seq(60, max_break, by = 10)
  sort(unique(c(b1, b2, b3, max_break + 1)))
}

#########################################################
# Helper: safely obtain chromosome & width columns from dge$genes
.get_chr <- function(genes_df) {
  cn <- colnames(genes_df)
  if ("seqnames" %in% cn) return(as.character(genes_df$seqnames))
  if ("chr" %in% cn)      return(as.character(genes_df$chr))
  if ("chrom" %in% cn)    return(as.character(genes_df$chrom))
  if ("chromosome" %in% cn) return(as.character(genes_df$chromosome))
  stop("Could not find chromosome column in dge$genes. Expected one of: seqnames/chr/chrom/chromosome.")
}

.get_width <- function(genes_df) {
  cn <- colnames(genes_df)
  if ("width" %in% cn) return(as.numeric(genes_df$width))
  if (all(c("start", "end") %in% cn)) {
    # GRanges-style start/end are 1-based inclusive; width = end-start+1
    return(as.numeric(genes_df$end) - as.numeric(genes_df$start) + 1)
  }
  stop("Could not find width in dge$genes (need 'width' or both 'start' and 'end').")
}

.get_symbol <- function(genes_df) {
  cn <- colnames(genes_df)
  if ("Symbol" %in% cn) return(as.character(genes_df$Symbol))
  if ("symbol" %in% cn) return(as.character(genes_df$symbol))
  return(NULL)
}

#########################################################
# Plot: counts per chromosome (regions and, if available, unique Symbols)
.plot_counts_per_chr <- function(dge, pdf_path,
                                 primary_only = TRUE,
                                 primary_regex = "^chr([1-9]|1[0-9]|X|Y|M)$") {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)

  genes_df <- as.data.frame(dge$genes)
  chr <- .get_chr(genes_df)
  sym <- .get_symbol(genes_df)

  keep <- rep(TRUE, length(chr))
  if (primary_only) keep <- grepl(primary_regex, chr)

  chr2 <- chr[keep]

  ## ---- NEW: define desired chromosome order
  primary_order <- paste0("chr", c(1:19, "X", "Y", "M"))
  # Keep only chromosomes that are actually present, in that order
  chr_present <- intersect(primary_order, unique(chr2))

  # Fallback: if nothing matches primary_order (e.g. non-standard names),
  # just use the unique chr2 in their original order
  if (length(chr_present) == 0L) {
    chr_present <- unique(chr2)
  }

  ## Regions per chromosome, ordered by chr_present (NOT by decreasing count)
  tab_regions <- table(chr2)
  regions_per_chr <- tab_regions[chr_present]

  # Optional unique gene count if Symbol exists
  genes_per_chr <- NULL
  if (!is.null(sym)) {
    sym2 <- sym[keep]
    df <- data.frame(chr = chr2, sym = sym2, stringsAsFactors = FALSE)
    df <- df[!is.na(df$sym) & df$sym != "", , drop = FALSE]

    # Count unique symbol per chromosome first
    genes_per_chr <- df %>%
      dplyr::distinct(chr, sym) %>%
      dplyr::count(chr, name = "n")

    # Reorder rows by chr_present so bars match regions_per_chr order
    genes_per_chr <- genes_per_chr[match(chr_present, genes_per_chr$chr), , drop = FALSE]
    genes_per_chr <- genes_per_chr[!is.na(genes_per_chr$chr), , drop = FALSE]
  }

  pdf(pdf_path, width = 11, height = 6)
  op <- par(mfrow = c(1, ifelse(is.null(genes_per_chr), 1, 2)),
            mar = c(8, 4, 4, 1) + 0.1)

  ## Barplot 1: regions per chromosome, in chr1..chr19, chrX, chrY, chrM order
  barplot(regions_per_chr,
          las = 2, cex.names = 0.8,
          main = "Number of regions per chromosome",
          ylab = "Regions")

  ## Barplot 2: unique Symbols per chromosome (if available), same order
  if (!is.null(genes_per_chr) && nrow(genes_per_chr) > 0) {
    barplot(setNames(genes_per_chr$n, genes_per_chr$chr),
            las = 2, cex.names = 0.8,
            main = "Number of unique Symbols per chromosome",
            ylab = "Unique genes (Symbol)")
  }

  par(op); dev.off()
  invisible(pdf_path)
}


#########################################################
# Plot: width distributions (overall + per chromosome mosaics)
.plot_width_distributions <- function(dge, pdf_overall, pdf_by_chr,
                                      primary_only = TRUE, primary_regex = "^chr([1-9]|1[0-9]|X|Y|M)$",
                                      log10_width = TRUE) {
  dir.create(dirname(pdf_overall), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(pdf_by_chr), recursive = TRUE, showWarnings = FALSE)

  genes_df <- as.data.frame(dge$genes)
  chr <- .get_chr(genes_df)
  w <- .get_width(genes_df)

  keep <- is.finite(w) & w >= 0
  if (primary_only) keep <- keep & grepl(primary_regex, chr)

  chr <- chr[keep]
  w <- w[keep]

  # Overall
  pdf(pdf_overall, width = 7, height = 5)
  if (log10_width) {
    x <- log10(w + 1)
    brks <- seq(0, max(x, na.rm = TRUE), length.out = 80)
    .hist_logy(x, breaks = brks, xlab = "log10(region width + 1)", main = "Width distribution (all regions)")
  } else {
    brks <- pretty(w, n = 80)
    .hist_logy(w, breaks = brks, xlab = "region width (bp)", main = "Width distribution (all regions)")
  }
  dev.off()

  # Per chromosome mosaic
  chr_levels <- names(sort(table(chr), decreasing = TRUE))
  grid <- .panel_grid(length(chr_levels))

  pdf(pdf_by_chr, width = 3.5 * grid[2], height = 3.0 * grid[1])
  op <- par(mfrow = c(grid[1], grid[2]), mar = c(4, 4, 3, 1) + 0.1)

  for (cc in chr_levels) {
    ww <- w[chr == cc]
    if (log10_width) {
      x <- log10(ww + 1)
      brks <- seq(0, max(x, na.rm = TRUE), length.out = 60)
      .hist_logy(x, breaks = brks, xlab = "log10(width+1)", main = cc)
    } else {
      brks <- pretty(ww, n = 60)
      .hist_logy(ww, breaks = brks, xlab = "width (bp)", main = cc)
    }
  }

  par(op); dev.off()
  invisible(list(overall = pdf_overall, by_chr = pdf_by_chr))
}

#########################################################
# Plot: depth distribution (reads per region) and sharedness (#samples with >=1 read)
.plot_depth_and_shared <- function(dge, output_dir, groups = c("Control", "MSUS"), max_break = 200) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cov_bio <- .bio_cov(dge)               # regions x bio samples
  bio_samples <- colnames(cov_bio)
  gmap <- .bio_group(dge, bio_samples)

  # ---- Depth: pooled across region×sample ----
  depths_all <- as.vector(cov_bio)
  brks_raw <- .depth_breaks(max_break)
  x_log <- log10(depths_all + 1)
  brks_log <- log10(brks_raw + 1)

  pdf_depth_all <- file.path(output_dir, "Depth_Distribution_AllSamples_REGION.pdf")
  pdf(pdf_depth_all, width = 10, height = 5)
  op <- par(mfrow = c(1, 2), mar = c(5, 4, 4, 1) + 0.1)

  .hist_logy(
    x = x_log, breaks = brks_log,
    xlab = "log10(reads per region + 1)",
    main = "Depth distribution (pooled region×sample)",
    ylab = "Count (log10)"
  )

  # Per-region totals across samples (gives “how covered is each region overall”)
  totals_per_region <- rowSums(cov_bio)
  .hist_logy(
    x = log10(totals_per_region + 1),
    breaks = seq(0, max(log10(totals_per_region + 1), na.rm = TRUE), length.out = 80),
    xlab = "log10(total reads per region + 1)",
    main = "Total reads per region (summed across samples)",
    ylab = "Regions (log10)"
  )

  par(op); dev.off()

  # ---- Depth: by group (one panel per group) ----
  pdf_depth_by_group <- file.path(output_dir, "Depth_Distribution_ByGroup_REGION.pdf")
  pdf(pdf_depth_by_group, width = 10, height = 5)
  op <- par(mfrow = c(1, length(groups)), mar = c(5, 4, 4, 1) + 0.1)

  for (g in groups) {
    keep <- names(gmap)[gmap == g]
    if (length(keep) == 0) { plot.new(); title(paste0(g, " (no samples)")); next }

    depths_g <- as.vector(cov_bio[, keep, drop = FALSE])
    .hist_logy(
      x = log10(depths_g + 1),
      breaks = brks_log,
      xlab = paste0("log10(reads+1) in ", g),
      main = paste0("Depth distribution (", g, ")"),
      ylab = "Count (log10)"
    )
  }
  par(op); dev.off()

  # ---- Sharedness: k = # bio samples with >=1 read per region ----
  k_all <- rowSums(cov_bio > 0)
  brks_k <- seq(-0.5, ncol(cov_bio) + 0.5, by = 1)

  pdf_shared_all <- file.path(output_dir, "SharedRegions_Distribution_AllSamples_REGION.pdf")
  pdf(pdf_shared_all, width = 7, height = 5)
  .hist_logy(
    x = k_all,
    breaks = brks_k,
    xlab = "Number of biological samples with ≥1 read in region",
    main = "Shared-region distribution (all samples)",
    ylab = "Regions (log10)"
  )
  dev.off()

  # ---- Sharedness: by group (k computed within each group) ----
  pdf_shared_by_group <- file.path(output_dir, "SharedRegions_Distribution_ByGroup_REGION.pdf")
  pdf(pdf_shared_by_group, width = 10, height = 5)
  op <- par(mfrow = c(1, length(groups)), mar = c(5, 4, 4, 1) + 0.1)

  for (g in groups) {
    keep <- names(gmap)[gmap == g]
    if (length(keep) == 0) { plot.new(); title(paste0(g, " (no samples)")); next }

    k_g <- rowSums(cov_bio[, keep, drop = FALSE] > 0)
    brks <- seq(-0.5, length(keep) + 0.5, by = 1)

    .hist_logy(
      x = k_g,
      breaks = brks,
      xlab = paste0("# samples with ≥1 read (", g, ")"),
      main = paste0("Shared-region distribution (", g, ")"),
      ylab = "Regions (log10)"
    )
  }

  par(op); dev.off()

  invisible(list(
    depth_all = pdf_depth_all,
    depth_by_group = pdf_depth_by_group,
    shared_all = pdf_shared_all,
    shared_by_group = pdf_shared_by_group
  ))
}

#########################################################
#' Region-wise QC: chromosome counts, width distributions, depth and sharedness
#'
#' Writes PDFs to output_dir:
#'  - ChrCounts_REGION.pdf
#'  - Width_AllRegions_REGION.pdf
#'  - Width_ByChromosome_REGION.pdf
#'  - Depth_Distribution_AllSamples_REGION.pdf
#'  - Depth_Distribution_ByGroup_REGION.pdf
#'  - SharedRegions_Distribution_AllSamples_REGION.pdf
#'  - SharedRegions_Distribution_ByGroup_REGION.pdf
#'
#' @param me_dge DGEList object OR path to an .rds containing the DGEList.
#' @param output_dir Output directory for PDFs.
#' @param groups Character vector of group names (default c("Control","MSUS")).
#' @param max_break Depth histogram max break (default 200).
#' @param primary_only Keep only chr1-19,X,Y,M by default.
#' @param primary_regex Regex for primary chroms.
#' @param log10_width Plot width as log10(width+1).
#' @return Invisibly returns a named list of output file paths.
#' @export
analyze_region_qc <- function(me_dge,
                              output_dir,
                              groups = c("Control", "MSUS"),
                              max_break = 200,
                              primary_only = TRUE,
                              primary_regex = "^chr([1-9]|1[0-9]|X|Y|M)$",
                              log10_width = TRUE) {

  dge <- .load_dge(me_dge)
  print(dge)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (is.null(dge$genes)) stop("dge$genes is required (must include chromosome + width/start/end).")
  if (is.null(dge$samples) || !("Group" %in% colnames(dge$samples))) {
    stop("dge$samples must contain a 'Group' column.")
  }

  out <- list()

  out$chr_counts <- file.path(output_dir, "ChrCounts_REGION.pdf")
  .plot_counts_per_chr(dge, out$chr_counts, primary_only = primary_only, primary_regex = primary_regex)

  out$width_all <- file.path(output_dir, "Width_AllRegions_REGION.pdf")
  out$width_by_chr <- file.path(output_dir, "Width_ByChromosome_REGION.pdf")
  .plot_width_distributions(
    dge,
    pdf_overall = out$width_all,
    pdf_by_chr  = out$width_by_chr,
    primary_only = primary_only,
    primary_regex = primary_regex,
    log10_width = log10_width
  )

  out2 <- .plot_depth_and_shared(dge, output_dir, groups = groups, max_break = max_break)
  out <- c(out, out2)

  message("Region-wise QC PDFs saved to: ", output_dir)
  invisible(out)
}


library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr) # Required for combining plots (ggarrange)
library(scales)
library(HDF5Array)

#' Generates QC plots for region-level DGEList.
#'
#' @param me_dge_region The region-level DGEList object (unfiltered).
#' @param cpg_counts_per_region Vector of CpG counts for each region.
#' @param region_type Character. The name of the region being processed (for naming files).
#' @param output_dir Character. Directory to save the plots.
#' @return NULL (plots are saved to disk).
plotRegionQC <- function(me_dge_region, cpg_counts_per_region, region_type, output_dir) {
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- Plot 1: Distribution of Number of CpGs per Region ---
  cpg_density_path <- file.path(output_dir, paste0("me_dge_", region_type, "_cpg_density.pdf"))
  pdf(cpg_density_path, width = 7, height = 7)
  hist(cpg_counts_per_region, 
       breaks = 50, 
       main = paste("CpG Count Distribution for", region_type),
       xlab = "Number of CpGs per Region",
       ylab = "Number of Regions",
       col = "skyblue", border = "black")
  legend("topright", paste("Total Regions:", length(cpg_counts_per_region)), bty = "n")
  dev.off()
  message(paste("Saved CpG density plot to:", cpg_density_path))
  
  # --- Plot 2: Median Read Coverage per Region ---
  
  # Calculate median read count per CpG site for each region
  # Total reads for a region: rowSums(me_dge_region$counts)
  # Total CpGs for a region: cpg_counts_per_region
  region_reads_total <- rowSums(me_dge_region$counts)
  
  # Calculate Median Reads per CpG (approximation of coverage)
  median_reads_per_cpg_in_region <- region_reads_total / cpg_counts_per_region
  
  coverage_path <- file.path(output_dir, paste0("me_dge_", region_type, "_coverage.pdf"))
  pdf(coverage_path, width = 7, height = 7)
  
  # Use log scale for better visualization of coverage distribution
  hist(log10(median_reads_per_cpg_in_region), 
       breaks = 50, 
       main = paste("Median Read Coverage (log10) for", region_type),
       xlab = expression("log"[10]*" (Median Reads per CpG in Region)"),
       ylab = "Number of Regions",
       col = "lightcoral", border = "black")
  dev.off()
  message(paste("Saved coverage plot to:", coverage_path))
}