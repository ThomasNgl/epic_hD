suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(HDF5Array)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(grid)
})

qc_me_symmetry <- function(se_input,
                           assay_name       = "Bval",   # 0 / 100 / NA assay
                           output_dir,
                           boxplot_filename = "QC_pairedCpG_perChr_boxplot.pdf",
                           tables_filename  = "QC_strand_symmetry_tables_perSample.pdf") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ## ----------------------------
  ## 1. Load SummarizedExperiment
  ## ----------------------------
  se <- if (inherits(se_input, "SummarizedExperiment")) {
    se_input
  } else if (is.character(se_input) && length(se_input) == 1) {
    if (dir.exists(se_input)) {
      HDF5Array::loadHDF5SummarizedExperiment(se_input)
    } else if (grepl("\\.rds$", se_input, ignore.case = TRUE)) {
      readRDS(se_input)
    } else {
      stop("se_input must be a SummarizedExperiment, an .rds file, or an HDF5SummarizedExperiment directory.")
    }
  } else {
    stop("Invalid se_input.")
  }

  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  A  <- SummarizedExperiment::assay(se, assay_name)
  if (!is.matrix(A)) A <- as.matrix(A)

  if (!all(c("Chr", "Locus", "Strand") %in% colnames(rd))) {
    stop("rowData(se) must contain columns: Chr, Locus, Strand.")
  }

  ## ----------------------------
  ## 2. Identify + / - strand rows and CpG pairs
  ## ----------------------------
  plus_idx  <- which(rd$Strand == "+")
  minus_idx <- which(rd$Strand == "-")

  if (!length(plus_idx) || !length(minus_idx)) {
    stop("qc_me_symmetry: no '+' or '-' strand rows found in rowData(se)$Strand.")
  }

  # plus: CpG coordinate = Locus
  plus_df <- data.frame(
    Chr      = as.character(rd$Chr[plus_idx]),
    CpG      = rd$Locus[plus_idx],
    idx_plus = plus_idx,
    stringsAsFactors = FALSE
  )

  # minus: CpG coordinate = Locus - 1 (C on minus strand is second base in CpG)
  minus_df <- data.frame(
    Chr       = as.character(rd$Chr[minus_idx]),
    CpG       = rd$Locus[minus_idx] - 1L,
    idx_minus = minus_idx,
    stringsAsFactors = FALSE
  )

  # inner join: keep only loci where + and - exist for same Chr & CpG
  pairs <- dplyr::inner_join(plus_df, minus_df, by = c("Chr", "CpG"))
  if (!nrow(pairs)) {
    stop("qc_me_symmetry: no plus/minus CpG pairs found (Chr + CpG coordinate).")
  }

  ## ----------------------------
  ## 3. PAIRED CpGs: per-sample & per-chr counts
  ## ----------------------------
  Ap <- A[pairs$idx_plus,  , drop = FALSE]  # plus rows, all samples
  Am <- A[pairs$idx_minus, , drop = FALSE]  # minus rows, all samples

  # a CpG dyad is "callable" in a sample if both + and - are non-NA
  obs     <- !is.na(Ap) & !is.na(Am)           # logical (n_pairs x n_samples)
  obs_int <- obs * 1L                          # 1 where TRUE, 0 where FALSE

  # aggregate by chromosome
  counts_mat <- rowsum(obs_int, group = pairs$Chr, reorder = FALSE)
  # rows = Chr, columns = samples

  df_pairs <- as.data.frame(counts_mat)
  df_pairs$Chr <- rownames(df_pairs)
  df_pairs_long <- tidyr::pivot_longer(
    df_pairs,
    cols      = -Chr,
    names_to  = "Sample",
    values_to = "N_paired_CpGs"
  )

  # order chromosomes nicely
  chr_levels <- paste0("chr", c(1:19, "X", "Y", "M"))
  df_pairs_long$Chr <- factor(df_pairs_long$Chr,
                              levels = intersect(chr_levels, unique(df_pairs_long$Chr)))

  ## ----------------------------
  ## 3b. UNPAIRED CpGs: indices (coordinate-based)
  ## ----------------------------
  pair_keys   <- paste(pairs$Chr, pairs$CpG)
  plus_keys   <- paste(plus_df$Chr, plus_df$CpG)
  minus_keys  <- paste(minus_df$Chr, minus_df$CpG)

  unpaired_plus_idx  <- plus_df$idx_plus[!(plus_keys  %in% pair_keys)]
  unpaired_minus_idx <- minus_df$idx_minus[!(minus_keys %in% pair_keys)]

  # extract unpaired assay submatrices
  A_unpaired_plus  <- A[unpaired_plus_idx,  , drop = FALSE]
  A_unpaired_minus <- A[unpaired_minus_idx, , drop = FALSE]

  ## ----------------------------
  ## 4. Colors for samples (for boxplot)
  ## ----------------------------
  maximally_distinct_colors <- function(n, seed = 1) {
    set.seed(seed)
    hs <- seq(0, 359, by = 1)
    cand <- c(
      grDevices::hcl(h = hs, c = 110, l = 55),
      grDevices::hcl(h = hs, c = 100, l = 70),
      grDevices::hcl(h = hs, c =  90, l = 40),
      "#000000", "#4D4D4D", "#808080", "#B3B3B3"
    )
    cand <- unique(cand)
    rgb <- t(grDevices::col2rgb(cand)) / 255
    start <- which.min(rowSums((rgb - matrix(c(1, 0, 0), nrow(rgb), 3, byrow = TRUE))^2))
    chosen <- start
    while (length(chosen) < n) {
      dmin <- rep(Inf, nrow(rgb))
      for (idx in chosen) {
        d <- rowSums((rgb - matrix(rgb[idx, ], nrow(rgb), 3, byrow = TRUE))^2)
        dmin <- pmin(dmin, d)
      }
      dmin[chosen] <- -Inf
      chosen <- c(chosen, which.max(dmin))
    }
    cand[chosen]
  }

  samples_all   <- sort(unique(df_pairs_long$Sample))
  sample_colors <- setNames(
    maximally_distinct_colors(length(samples_all), seed = 1),
    samples_all
  )

  ## ----------------------------
  ## 5. Boxplot: number of paired CpGs per Chr across samples (with sample colours)
  ## ----------------------------
  p_box <- ggplot2::ggplot(df_pairs_long,
                           ggplot2::aes(x = Chr, y = N_paired_CpGs)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey80", colour = "black") +
    ggplot2::geom_jitter(ggplot2::aes(color = Sample),
                         width = 0.2, alpha = 0.7, size = 1) +
    ggplot2::scale_color_manual(values = sample_colors, breaks = names(sample_colors)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title  = "Callable CpG dyads per chromosome",
      y      = "Number of CpG dyads (non-NA on + and -)",
      x      = "Chromosome",
      color  = "Sample"
    )

  grDevices::pdf(file.path(output_dir, boxplot_filename), width = 10, height = 6)
  print(p_box)
  grDevices::dev.off()

  ## ----------------------------
  ## 6. Per-sample tables:
  ##    - 2×2 paired plus/minus M/U
  ##    - 2×2 unpaired strand (rows = + / -, cols = U / M)
  ## ----------------------------
  classify_state <- function(x) {
    # expects 0 / 100 / NA; anything else becomes NA
    ifelse(
      is.na(x),
      NA_character_,
      ifelse(
        x == 100, "M",
        ifelse(x == 0, "U", NA_character_)
      )
    )
  }

  sample_tables_paired   <- list()
  sample_tables_unpaired <- list()
  sample_names           <- colnames(Ap)

  for (j in seq_along(sample_names)) {
    ## ----- paired states -----
    ps <- classify_state(Ap[, j])  # plus strand on paired CpGs
    ms <- classify_state(Am[, j])  # minus strand on paired CpGs

    valid <- !is.na(ps) & !is.na(ms)
    if (any(valid)) {
      tab <- table(plus = ps[valid], minus = ms[valid])

      lev <- c("U", "M")
      mat <- matrix(0L, nrow = 2, ncol = 2,
                    dimnames = list(plus = lev, minus = lev))
      mat[rownames(tab), colnames(tab)] <- tab
    } else {
      # no fully paired callable CpGs for this sample
      mat <- matrix(0L, nrow = 2, ncol = 2,
                    dimnames = list(plus = c("U","M"), minus = c("U","M")))
    }
    sample_tables_paired[[sample_names[j]]] <- mat

    ## ----- unpaired states -----
    # plus strand only (no matching minus CpG)
    ps_u <- if (length(unpaired_plus_idx)) {
      classify_state(A_unpaired_plus[, j])
    } else {
      character(0)
    }

    # minus strand only (no matching plus CpG)
    ms_u <- if (length(unpaired_minus_idx)) {
      classify_state(A_unpaired_minus[, j])
    } else {
      character(0)
    }

    n_plus_U  <- sum(ps_u == "U", na.rm = TRUE)
    n_plus_M  <- sum(ps_u == "M", na.rm = TRUE)
    n_minus_U <- sum(ms_u == "U", na.rm = TRUE)
    n_minus_M <- sum(ms_u == "M", na.rm = TRUE)

    # 2×2, rows = strand (+ / -), cols = state (U / M)
    unpaired_mat <- matrix(
      c(n_plus_U, n_plus_M,
        n_minus_U, n_minus_M),
      nrow = 2, byrow = TRUE,
      dimnames = list(
        Strand = c("+", "-"),
        State  = c("U", "M")
      )
    )

    sample_tables_unpaired[[sample_names[j]]] <- unpaired_mat
  }

  ## ----------------------------
  ## 7. Save tables to a PDF (one page per sample)
  ##    Layout: sample name, paired 2×2, unpaired 2×2.
  ## ----------------------------
  tables_path <- file.path(output_dir, tables_filename)
  grDevices::pdf(tables_path, width = 6, height = 6)
  for (s in sample_names) {
    paired_tab   <- sample_tables_paired[[s]]
    unpaired_tab <- sample_tables_unpaired[[s]]

    grob_title    <- grid::textGrob(s, gp = grid::gpar(fontface = "bold"))
    grob_paired   <- gridExtra::tableGrob(paired_tab)
    grob_unpaired <- gridExtra::tableGrob(unpaired_tab)
    grob_paired   <- gridExtra::tableGrob(paired_tab)

    # Add a small caption clarifying dimensions of the paired table
    grob_paired <- gridExtra::arrangeGrob(
      grid::textGrob("Paired CpGs: rows = plus strand, cols = minus strand",
                    x = 1, y = 1, just = c("right", "top"),
                    gp = grid::gpar(cex = 0.6)),
      grob_paired,
      ncol = 1,
      heights = c(0.2, 0.8)
    )

    grob_unpaired <- gridExtra::tableGrob(unpaired_tab)

    gridExtra::grid.arrange(
      grob_title,
      grob_paired,
      grob_unpaired,
      ncol    = 1,
      heights = c(0.15, 0.425, 0.425)
    )
  }
  grDevices::dev.off()

  invisible(list(
    pairs_per_chr          = df_pairs_long,          # Chr / Sample / N_paired_CpGs
    tables_paired_per_sample   = sample_tables_paired,   # list of 2×2 (paired)
    tables_unpaired_per_sample = sample_tables_unpaired, # list of 2×2 (unpaired by strand)
    boxplot_pdf            = file.path(output_dir, boxplot_filename),
    tables_pdf             = tables_path
  ))
}
