library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr) # Required for combining plots (ggarrange)
library(scales)
library(HDF5Array)

#' Plot per-sample MePer histograms from a SummarizedExperiment
#'
#' @param se_path Character path to SummarizedExperiment RDS OR a SummarizedExperiment object.
#' @param output_dir Directory to save plots.
#' @param assay_name Assay name containing MePer (default "MePer").
#' @export
plot_all_sample_mepe_se <- function(se_path, output_dir, assay_name = "MePer") {

  # ---- load SE depending on input ----
  se <- if (inherits(se_path, "SummarizedExperiment")) {
    se_path
  } else if (is.character(se_path) && length(se_path) == 1) {
    if (dir.exists(se_path)) {
      HDF5Array::loadHDF5SummarizedExperiment(se_path)   # HDF5SummarizedExperiment directory
    } else if (grepl("\\.rds$", se_path, ignore.case = TRUE)) {
      readRDS(se_path)                                   # regular RDS SE
    } else {
      stop("se_path must be a SummarizedExperiment, an .rds file, or an HDF5SummarizedExperiment directory.")
    }
  } else {
    stop("Invalid se_path.")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  M  <- SummarizedExperiment::assay(se, assay_name)
  ph <- as.data.frame(SummarizedExperiment::colData(se))
  if (!"Group" %in% names(ph)) ph$Group <- NA

  BIN_WIDTH <- 5

  # order samples by Group then name
  sample_order <- rownames(ph)[order(ph$Group, rownames(ph))]
  plot_list <- list()
    # ---- NEW: group-concatenated distributions + violin + KS test ----
  # Collect MePer values pooled within each group (streaming columns from HDF5/DelayedArray).
  # To keep PDFs reasonable, we downsample per group (set to Inf to disable).
  MAX_POINTS_PER_GROUP <- 200000L
  set.seed(1L)

  groups_keep <- c("Control", "MSUS")
  sample_in_groups <- rownames(ph)[ph$Group %in% groups_keep]
  if (length(sample_in_groups) > 0) {
    pooled <- vector("list", length(groups_keep))
    names(pooled) <- groups_keep

    for (s in sample_in_groups) {
      g <- as.character(ph[s, "Group"])
      x <- as.numeric(M[, s])
      x <- x[is.finite(x)]
      if (!length(x)) next
      pooled[[g]] <- c(pooled[[g]], x)
    }

    # Downsample per group (optional)
    pooled_ds <- lapply(pooled, function(x) {
      if (!length(x)) return(x)
      if (is.finite(MAX_POINTS_PER_GROUP) && length(x) > MAX_POINTS_PER_GROUP) {
        sample(x, MAX_POINTS_PER_GROUP)
      } else x
    })

    df_group <- data.frame(
      Group = rep(names(pooled_ds), times = vapply(pooled_ds, length, integer(1))),
      MePer = unlist(pooled_ds, use.names = FALSE)
    )

    # KS test (distributional difference). If either group is empty, skip.
    p_ks <- NA_real_
    if (all(vapply(pooled_ds, length, integer(1)) > 0)) {
      p_ks <- suppressWarnings(stats::ks.test(pooled_ds$Control, pooled_ds$MSUS)$p.value)
    }

    p_ks <- suppressWarnings(stats::ks.test(pooled_ds$Control, pooled_ds$MSUS)$p.value)

    p_group_violin <- ggplot2::ggplot(df_group, ggplot2::aes(x = Group, y = MePer, fill = Group)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
      ggplot2::scale_fill_manual(values = c(Control = "steelblue", MSUS = "firebrick")) +
      ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(
        x = NULL,
        y = "CpG methylation (%)",
        title = "Pooled CpG methylation distributions by group"
      ) +
      ggplot2::annotate(
        "text", x = 1.5, y = 98,
        label = if (is.na(p_ks)) "KS p = NA" else paste0("KS p = ", signif(p_ks, 3))
      )

    grDevices::pdf(file.path(output_dir, "Violin_pooled_CpG_methylation_by_group.pdf"), width = 6, height = 4)
    print(p_group_violin); grDevices::dev.off()
  }
  for (s in sample_order) {
    x <- as.numeric(M[, s])          # works with DelayedArray/HDF5, pulls one column at a time
    x <- x[is.finite(x)]
    if (!length(x)) next

    median_meper <- median(x, na.rm = TRUE)
    df <- data.frame(MePer = x)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = MePer, y = ggplot2::after_stat(density) * BIN_WIDTH * 100)) +
      ggplot2::geom_histogram(binwidth = BIN_WIDTH, fill = "lightblue", color = "black", alpha = 0.8) +
      ggplot2::scale_y_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, 20),
        labels = scales::label_percent(scale = 1)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = s,
        subtitle = paste("Group:", ph[s, "Group"]),
        x = "CpG Methylation (%)",
        y = "Frequency (%)"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 10),
        plot.subtitle = ggplot2::element_text(size = 8),
        axis.text = ggplot2::element_text(size = 7),
        axis.title = ggplot2::element_text(size = 9)
      ) +
      ggplot2::geom_vline(xintercept = median_meper, linetype = "dashed", color = "red")

    grDevices::pdf(file.path(output_dir, paste0("MePer_Distribution_", s, ".pdf")), width = 5, height = 4)
    print(p); grDevices::dev.off()

    plot_list[[s]] <- p
  }

  # mosaic
  if (length(plot_list)) {
    n <- length(plot_list)
    n_cols <- ceiling(sqrt(n))
    n_rows <- ceiling(n / n_cols)

    combined <- ggpubr::ggarrange(
      plotlist = plot_list, ncol = n_cols, nrow = n_rows,
      font.label = list(size = 8)
    )

    grDevices::pdf(file.path(output_dir, "MePer_Distribution_Mosaic_All_Samples.pdf"),
                   width = 2.5 * n_cols, height = 2.5 * n_rows)
    print(combined); grDevices::dev.off()
  }

  # ---- extra summary boxplots (per sample) ----
  # NOTE: colMeans() on an HDF5/DelayedArray can be expensive; this computes by columns safely.
  pct_meth   <- vapply(colnames(M), function(s) mean(as.numeric(M[, s]) >= 80, na.rm = TRUE) * 100, numeric(1))
  pct_unmeth <- vapply(colnames(M), function(s) mean(as.numeric(M[, s]) <= 20, na.rm = TRUE) * 100, numeric(1))

  ### NEW: mean methylation per sample (0–100 scale)
  mean_meper <- vapply(colnames(M), function(s) mean(as.numeric(M[, s]), na.rm = TRUE), numeric(1))

  df_sum <- data.frame(
    Sample = colnames(M),
    Group  = ph[colnames(M), "Group"],
    pct_methylated   = pct_meth,
    pct_unmethylated = pct_unmeth,
    mean_meper       = mean_meper      # NEW
  )

  # % CpGs methylated
  p_meth <- ggplot2::ggplot(df_sum, ggplot2::aes(x = Group, y = pct_methylated)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "% CpGs methylated (MePer ≥ 80)") +
    ggpubr::stat_compare_means(method = "wilcox.test")

  # % CpGs unmethylated
  p_unmeth <- ggplot2::ggplot(df_sum, ggplot2::aes(x = Group, y = pct_unmethylated)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "% CpGs unmethylated (MePer ≤ 20)") +
    ggpubr::stat_compare_means(method = "wilcox.test")

  grDevices::pdf(file.path(output_dir, "Boxplot_pct_methylated_by_group.pdf"), width = 6, height = 4)
  print(p_meth); grDevices::dev.off()

  grDevices::pdf(file.path(output_dir, "Boxplot_pct_unmethylated_by_group.pdf"), width = 6, height = 4)
  print(p_unmeth); grDevices::dev.off()

  ### NEW: boxplot of mean methylation by group
  p_mean <- ggplot2::ggplot(df_sum, ggplot2::aes(x = Group, y = mean_meper)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "Mean CpG methylation (%)") +
    ggpubr::stat_compare_means(method = "wilcox.test")

  grDevices::pdf(file.path(output_dir, "Boxplot_mean_methylation_by_group.pdf"), width = 6, height = 4)
  print(p_mean); grDevices::dev.off()

  invisible(list(sample_plots = plot_list, summary = df_sum))
}

