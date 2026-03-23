library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(scales)
library(forcats)

#' Analyze and Plot CpG Coverage Distribution per Chromosome (memory-safe)
#'
#' Computes per-sample CpG “coverage” (feature present if methylated>0 OR unmethylated>0)
#' and summarizes it by chromosome. Handles two DGEList formats:
#' 1) Separate $unmethylated slot
#' 2) Single $counts matrix with "-Me"/"-Un" column suffixes
#'
#' Outputs (PDF + CSV):
#' - Coverage_Composition_Stacked_PerSample.pdf
#'     Stacked bars of fraction of covered CpGs contributed by each chromosome (per sample).
#' - Coverage_Composition_Stacked_PerSample_RefNCpGNorm.pdf
#'     Stacked bars normalized by reference CpG counts per chromosome (CpGs_Covered/Total_CpGs_In_Ref),
#'     then renormalized within sample to sum to 1.
#' - Coverage_Breadth_Boxplot_PerChr.pdf
#'     Boxplots per chromosome of % dataset-chr CpGs covered per sample.
#' - Coverage_Breadth_Boxplot_PerChr_RefNCpGNorm_PStrand.pdf
#' - Coverage_Breadth_Boxplot_PerChr_RefNCpGNorm_NStrand.pdf
#'     Same as above but % reference-chr CpGs covered, split by Strand from dge$genes$Strand.
#' - Coverage_Breadth_Boxplot_ByGroup_PerChr.pdf + Coverage_Breadth_Wilcoxon_PerChr.csv
#'     Group-wise breadth differences (Wilcoxon per chr, BH-FDR).
#'
#' Chromosome ordering is fixed for stacked bars and boxplots:
#' chr1..chr19, chrX, chrY, chrM (plus any extra contigs appended).
#' Bar fill colors are assigned to maximize contrast between adjacent chromosomes while
#' preserving the above stacking order.
#'
#' @param me_dge DGEList object OR path to an .rds containing the DGEList.
#' @param output_dir Character. Directory where outputs are written.
#' @param nCpGs_per_chr_path Character. CSV with reference CpG counts per chromosome.
#'   Required columns: chromosome, n_cpgs.
#' @return Invisibly returns a list of data.frames used for plotting + Wilcoxon table.
#' @export
analyze_cpg_coverage_by_chr <- function(me_dge,
                                       output_dir,
                                       nCpGs_per_chr_path,
                                       covariate_shape = "Batch",
                                       strand_dup_factor = 2) {

  # ---- Load DGEList (path or object) ----
  if (is.character(me_dge) && length(me_dge) == 1L) {
    if (!file.exists(me_dge)) stop("RDS file not found at: ", me_dge)
    dge <- readRDS(me_dge)
  } else {
    dge <- me_dge
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ---- Basic requirements ----
  if (is.null(dge$genes) || !all(c("Chr", "Strand") %in% colnames(dge$genes))) {
    stop("dge$genes must contain columns: 'Chr' and 'Strand'.")
  }

  chromosomes <- as.character(dge$genes$Chr)
  strands     <- as.character(dge$genes$Strand)

  has_merged <- any(strands == "merged", na.rm = TRUE)
  # ---- helper: natural chromosome ordering (numeric then alpha) ----
  chr_levels_natural <- function(chr_vec) {
    chr_vec <- unique(as.character(chr_vec))
    core <- sub("^chr", "", chr_vec, ignore.case = TRUE)
    is_num <- grepl("^[0-9]+$", core)
    num_val <- suppressWarnings(as.integer(core))
    alpha_val <- tolower(core)

    ord <- order(!is_num, ifelse(is_num, num_val, Inf), ifelse(is_num, "", alpha_val))
    chr_vec[ord]
  }

  # ---- Build coverage matrix (logical) ----
  if (!is.null(dge$unmethylated)) {
    is_covered_mat <- (dge$counts > 0) | (dge$unmethylated > 0)
  } else {
    Me_cols <- grep("-Me$", colnames(dge$counts), value = TRUE)
    if (!length(Me_cols)) stop("No '-Me' columns found and no $unmethylated slot present.")
    Un_cols <- sub("-Me$", "-Un", Me_cols)

    ok <- Un_cols %in% colnames(dge$counts)
    if (!all(ok)) {
      warning("Some matching '-Un' columns are missing; using only complete Me/Un pairs.")
      Me_cols <- Me_cols[ok]; Un_cols <- Un_cols[ok]
    }
    if (!length(Me_cols)) stop("No complete -Me/-Un pairs found.")

    Mat_Me <- dge$counts[, Me_cols, drop = FALSE]
    Mat_Un <- dge$counts[, Un_cols, drop = FALSE]

    is_covered_mat <- (Mat_Me > 0) | (Mat_Un > 0)
    colnames(is_covered_mat) <- sub("-Me$", "", Me_cols)
  }

  if (nrow(is_covered_mat) != length(chromosomes)) {
    stop("Gene annotation (rows) and count matrix dimensions mismatch.")
  }
  n_features <- nrow(is_covered_mat)

  # ---- Aggregate covered CpGs per chromosome per sample ----
  counts_per_chr_sample <- rowsum(is_covered_mat * 1L, group = chromosomes, reorder = FALSE)

  if (has_merged) {
    counts_per_chr_sample_merged <- counts_per_chr_sample
  } else {
    counts_per_chr_sample_plus  <- rowsum(is_covered_mat[strands == "+", , drop = FALSE] * 1L,
                                          group = chromosomes[strands == "+"], reorder = FALSE)
    counts_per_chr_sample_minus <- rowsum(is_covered_mat[strands == "-", , drop = FALSE] * 1L,
                                          group = chromosomes[strands == "-"], reorder = FALSE)
  }

  # Keep overall per-sample totals BEFORE dropping the big matrix
  total_cov_per_sample <- colSums(counts_per_chr_sample)

  rm(is_covered_mat); gc()

  # ---- Convert to long format: Chr / Sample / CpGs_Covered ----
  df_long <- as.data.frame(counts_per_chr_sample) %>%
    tibble::rownames_to_column("Chr") %>%
    tidyr::pivot_longer(-Chr, names_to = "Sample", values_to = "CpGs_Covered")

  if (has_merged) {
    df_long_merged <- as.data.frame(counts_per_chr_sample_merged) %>%
      tibble::rownames_to_column("Chr") %>%
      tidyr::pivot_longer(-Chr, names_to = "Sample", values_to = "CpGs_Covered") %>%
      dplyr::mutate(Strand = "merged")
  } else {
    df_long_plus <- as.data.frame(counts_per_chr_sample_plus) %>%
      tibble::rownames_to_column("Chr") %>%
      tidyr::pivot_longer(-Chr, names_to = "Sample", values_to = "CpGs_Covered") %>%
      dplyr::mutate(Strand = "+")

    df_long_minus <- as.data.frame(counts_per_chr_sample_minus) %>%
      tibble::rownames_to_column("Chr") %>%
      tidyr::pivot_longer(-Chr, names_to = "Sample", values_to = "CpGs_Covered") %>%
      dplyr::mutate(Strand = "-")
}

  # ---- Load reference CpG counts per chromosome ----
  if (!file.exists(nCpGs_per_chr_path)) stop("nCpGs_per_chr_path not found: ", nCpGs_per_chr_path)

  chr_ref_counts <- readr::read_csv(nCpGs_per_chr_path, show_col_types = FALSE) %>%
    dplyr::rename(Chr = chromosome, Total_CpGs_In_Ref = n_cpgs) %>%
    dplyr::mutate(Chr = as.character(Chr))

  ref_denom <- sum(chr_ref_counts$Total_CpGs_In_Ref, na.rm = TRUE)

  # ---- Chromosome order + colors (natural order) ----
  chr_levels <- chr_levels_natural(df_long$Chr)

  df_long$Chr <- factor(df_long$Chr, levels = chr_levels)

  if (has_merged) {
    df_long_merged$Chr <- factor(df_long_merged$Chr, levels = chr_levels)
  } else {
    df_long_plus$Chr  <- factor(df_long_plus$Chr, levels = chr_levels)
    df_long_minus$Chr <- factor(df_long_minus$Chr, levels = chr_levels)
  }

  base_pal <- scales::hue_pal(l = 65, c = 100)(length(chr_levels))
  ord <- c(seq(1, length(chr_levels), by = 2), seq(2, length(chr_levels), by = 2))
  chr_colors <- setNames(base_pal[ord], chr_levels)

  # ---- Plot 1: composition (feature-based) ----
  df_composition <- df_long %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Total_Covered_In_Sample = sum(CpGs_Covered)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Prop_Features_From_Chr = CpGs_Covered / Total_Covered_In_Sample)

  df_composition$Chr <- factor(df_composition$Chr, levels = chr_levels) |> forcats::fct_drop()

  p1 <- ggplot(df_composition, aes(x = Sample, y = Prop_Features_From_Chr, fill = Chr)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = chr_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
    labs(
      title = "Genomic Distribution of Covered CpGs per Sample",
      subtitle = "Per-sample composition based on covered features (CpGs_Covered).",
      y = "Proportion of Covered CpGs",
      x = "Sample",
      fill = "Chr"
    )

  # ---- Plot 1b: composition normalized by reference CpG counts ----
  df_comp_ref_norm <- df_long %>%
    dplyr::left_join(chr_ref_counts, by = "Chr") %>%
    dplyr::filter(!is.na(Total_CpGs_In_Ref)) %>%
    dplyr::mutate(CoverageRate = CpGs_Covered / Total_CpGs_In_Ref) %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Prop_CoverageRate = CoverageRate / sum(CoverageRate)) %>%
    dplyr::ungroup()

  df_comp_ref_norm$Chr <- factor(df_comp_ref_norm$Chr, levels = chr_levels) |> forcats::fct_drop()

  p1_ref_norm <- ggplot(df_comp_ref_norm, aes(x = Sample, y = Prop_CoverageRate, fill = Chr)) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = chr_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
    labs(
      title = "Genomic Distribution of Covered CpGs per Sample (Reference CpG-normalized)",
      subtitle = "Uses CpGs_Covered / Total_CpGs_In_Ref per chr; then renormalizes within sample.",
      y = "Proportion of CoverageRate",
      x = "Sample",
      fill = "Chr"
    )

  # ---- Breadth vs dataset CpG universe (per chromosome) ----
  chr_dataset_counts <- as.data.frame(table(chromosomes), stringsAsFactors = FALSE)
  colnames(chr_dataset_counts) <- c("Chr", "Total_CpGs_In_Dataset")

  df_breadth <- df_long %>%
    dplyr::left_join(chr_dataset_counts, by = "Chr") %>%
    dplyr::mutate(Pct_Universe_Covered = 100 * CpGs_Covered / Total_CpGs_In_Dataset)

  # ---- Breadth vs reference CpG universe (by strand) ----
  if (has_merged) {

    df_breadth_ref_merged <- df_long_merged %>%
      dplyr::left_join(chr_ref_counts, by = "Chr") %>%
      dplyr::filter(!is.na(Total_CpGs_In_Ref)) %>%
      dplyr::mutate(Pct_RefChr_Covered = 100 * CpGs_Covered / Total_CpGs_In_Ref)

  } else {

    df_breadth_ref_plus <- df_long_plus %>%
      dplyr::left_join(chr_ref_counts, by = "Chr") %>%
      dplyr::filter(!is.na(Total_CpGs_In_Ref)) %>%
      dplyr::mutate(Pct_RefChr_Covered = 100 * CpGs_Covered / Total_CpGs_In_Ref)

    df_breadth_ref_minus <- df_long_minus %>%
      dplyr::left_join(chr_ref_counts, by = "Chr") %>%
      dplyr::filter(!is.na(Total_CpGs_In_Ref)) %>%
      dplyr::mutate(Pct_RefChr_Covered = 100 * CpGs_Covered / Total_CpGs_In_Ref)
  }

  # ---- Sample metadata (Group + optional covariate_shape) ----
  if (is.null(dge$samples) || !("Group" %in% colnames(dge$samples))) {
    stop("dge$samples must contain a 'Group' column (e.g., MSUS / Control).")
  }

  sample_meta <- tibble::as_tibble(as.data.frame(dge$samples))
  if (!("Sample" %in% colnames(sample_meta))) {
    sample_meta <- sample_meta %>% tibble::rownames_to_column("Sample")
  }

  sample_meta <- sample_meta %>%
    dplyr::mutate(Sample = as.character(Sample), Group = as.character(Group))

  if (!is.null(covariate_shape)) {
    if (!covariate_shape %in% colnames(sample_meta)) {
      stop("covariate_shape not found in dge$samples: ", covariate_shape)
    }
    sample_meta <- sample_meta %>%
      dplyr::mutate(CovariateShape = as.character(.data[[covariate_shape]]))
  } else {
    sample_meta$CovariateShape <- NA_character_
  }

  sample_meta <- sample_meta %>%
    dplyr::select(Sample, Group, CovariateShape) %>%
    dplyr::distinct()

  df_breadth_g <- df_breadth %>% dplyr::left_join(sample_meta, by = "Sample")
  if (any(is.na(df_breadth_g$Group))) {
    warning("Some samples have no Group after join; dropping them for group plots/tests.")
    df_breadth_g <- df_breadth_g %>% dplyr::filter(!is.na(Group))
  }

  grp_levels <- sort(unique(as.character(df_breadth_g$Group)))
  if (length(grp_levels) != 2) stop("Expected exactly 2 groups, found: ", paste(grp_levels, collapse = ", "))
  df_breadth_g$Group <- factor(df_breadth_g$Group, levels = grp_levels)

  if (any(is.na(df_breadth_g$CovariateShape))) df_breadth_g$CovariateShape[is.na(df_breadth_g$CovariateShape)] <- "NA"
  df_breadth_g$CovariateShape <- factor(df_breadth_g$CovariateShape)

  # ---- sample colors (as before) ----
  if (all(c("Control", "MSUS") %in% levels(df_breadth_g$Group))) {
    light_group <- "Control"; dark_group <- "MSUS"
  } else {
    light_group <- levels(df_breadth_g$Group)[1]; dark_group <- levels(df_breadth_g$Group)[2]
  }

  mix_hex <- function(hex, mix_with = c("white", "black"), amount = 0.45) {
    mix_with <- match.arg(mix_with)
    rgb <- grDevices::col2rgb(hex) / 255
    target <- if (mix_with == "white") c(1, 1, 1) else c(0, 0, 0)
    out <- rgb * (1 - amount) + target * amount
    grDevices::rgb(out[1], out[2], out[3])
  }

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

  samples_all <- sort(unique(df_breadth_g$Sample))
  base_cols <- setNames(maximally_distinct_colors(length(samples_all), seed = 1), samples_all)

  group_by_sample <- df_breadth_g %>% dplyr::distinct(Sample, Group)
  sample_colors <- base_cols
  for (i in seq_len(nrow(group_by_sample))) {
    s <- group_by_sample$Sample[i]
    g <- as.character(group_by_sample$Group[i])
    if (g == light_group) sample_colors[s] <- mix_hex(base_cols[s], "white", amount = 0.45)
    if (g == dark_group)  sample_colors[s] <- mix_hex(base_cols[s], "black", amount = 0.35)
  }

  # ---- Percent axis control for boxplots ----
  y_breaks <- seq(0, 100, by = 25)
  y_limits <- c(0, 102)

  # ---- Boxplot: dataset-universe breadth ----
  p2 <- ggplot(df_breadth, aes(x = Chr, y = Pct_Universe_Covered)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
    geom_jitter(aes(color = Sample), width = 0.2, size = 1, alpha = 0.75) +
    scale_color_manual(values = sample_colors, breaks = names(sample_colors)) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
    labs(
      title = "Breadth of Coverage per Chromosome",
      subtitle = "100 × CpGs_Covered / Total_CpGs_In_Dataset (per chromosome).",
      y = "% of Dataset CpGs Covered",
      x = "Chromosome",
      color = "Sample"
    )

  if (has_merged) {

    p3_ref_breadth_merged <- ggplot(df_breadth_ref_merged, aes(x = Chr, y = Pct_RefChr_Covered)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
      geom_jitter(aes(color = Sample), width = 0.2, size = 1, alpha = 0.75) +
      scale_color_manual(values = sample_colors, breaks = names(sample_colors)) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
      labs(
        title = "Breadth of Coverage per Chromosome (Reference CpG-normalized, merged strand)",
        subtitle = "100 × CpGs_Covered / Total_CpGs_In_Ref (all rows; strands merged).",
        y = "% of Reference Chr CpGs Covered",
        x = "Chromosome",
        color = "Sample"
      )

  } else {

    p3_ref_breadth_plus <- ggplot(df_breadth_ref_plus, aes(x = Chr, y = Pct_RefChr_Covered)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
      geom_jitter(aes(color = Sample), width = 0.2, size = 1, alpha = 0.75) +
      scale_color_manual(values = sample_colors, breaks = names(sample_colors)) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
      labs(
        title = "Breadth of Coverage per Chromosome (Reference CpG-normalized, + strand)",
        subtitle = "100 × CpGs_Covered / Total_CpGs_In_Ref (+ strand rows only).",
        y = "% of Reference Chr CpGs Covered",
        x = "Chromosome",
        color = "Sample"
      )

    p3_ref_breadth_minus <- ggplot(df_breadth_ref_minus, aes(x = Chr, y = Pct_RefChr_Covered)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
      geom_jitter(aes(color = Sample), width = 0.2, size = 1, alpha = 0.75) +
      scale_color_manual(values = sample_colors, breaks = names(sample_colors)) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
      labs(
        title = "Breadth of Coverage per Chromosome (Reference CpG-normalized, - strand)",
        subtitle = "100 × CpGs_Covered / Total_CpGs_In_Ref (- strand rows only).",
        y = "% of Reference Chr CpGs Covered",
        x = "Chromosome",
        color = "Sample"
      )
  }

  # ---- Group-wise per chromosome + Wilcoxon ----
  wilcox_by_chr <- df_breadth_g %>%
    dplyr::group_by(Chr) %>%
    dplyr::summarize(
      n = dplyr::n(),
      p_value = wilcox.test(Pct_Universe_Covered ~ Group, exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

  ann_df <- wilcox_by_chr %>%
    dplyr::mutate(label = paste0("p=", signif(p_value, 3), "\nFDR=", signif(p_adj_BH, 3)))

  jitter_aes <- if (is.null(covariate_shape)) aes(color = Sample) else aes(color = Sample, shape = CovariateShape)

  p2_by_group <- ggplot(df_breadth_g, aes(x = Group, y = Pct_Universe_Covered)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(jitter_aes, width = 0.18, height = 0, size = 1, alpha = 0.75) +
    facet_wrap(~ Chr, scales = "free_y") +
    scale_color_manual(values = sample_colors, breaks = names(sample_colors)) +
    { if (!is.null(covariate_shape)) scale_shape_discrete(name = covariate_shape) } +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right") +
    labs(
      title = "Breadth of Coverage by Group (per Chromosome)",
      subtitle = "Control vs MSUS; Wilcoxon per chr (BH-FDR).",
      x = "Group",
      y = "% of Dataset CpGs Covered",
      color = "Sample"
    ) +
    geom_label(
      data = ann_df,
      aes(x = 1.5, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.1,
      size = 2.6,
      label.size = 0,
      alpha = 0.75
    ) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits)

  # ===========================================================================
  # NEW: Overall coverage breadth (dataset vs reference) by group (CSV + PDF)
  # ===========================================================================

  df_overall <- tibble::tibble(
    Sample = names(total_cov_per_sample),
    CpGs_Covered = as.numeric(total_cov_per_sample)
  ) %>%
    dplyr::left_join(sample_meta, by = "Sample") %>%
    dplyr::filter(!is.na(Group)) %>%
    dplyr::mutate(
      Group = factor(as.character(Group), levels = levels(df_breadth_g$Group)),
      Dataset_Coverage_Pct   = 100 * CpGs_Covered / n_features,
      Reference_Coverage_Pct = 100 * (CpGs_Covered / strand_dup_factor) / ref_denom
    )

  out_overall_csv <- file.path(output_dir, "Coverage_Breadth_Overall_ByGroup.csv")
  write.csv(
    df_overall %>% dplyr::select(Sample, Group, CovariateShape, Dataset_Coverage_Pct, Reference_Coverage_Pct),
    out_overall_csv,
    row.names = FALSE
  )

  p_dataset <- wilcox.test(Dataset_Coverage_Pct ~ Group, data = df_overall, exact = FALSE)$p.value
  p_ref     <- wilcox.test(Reference_Coverage_Pct ~ Group, data = df_overall, exact = FALSE)$p.value

  df_overall_long <- df_overall %>%
    dplyr::select(Sample, Group, CovariateShape, Dataset_Coverage_Pct, Reference_Coverage_Pct) %>%
    tidyr::pivot_longer(
      cols = c(Dataset_Coverage_Pct, Reference_Coverage_Pct),
      names_to = "Metric",
      values_to = "Coverage_Pct"
    ) %>%
    dplyr::mutate(
      Metric = dplyr::recode(
        Metric,
        Dataset_Coverage_Pct   = "Coverage vs dataset CpG universe",
        Reference_Coverage_Pct = "Coverage vs reference CpG universe"
      )
    )

  ann_overall <- tibble::tibble(
    Metric = c("Coverage vs dataset CpG universe", "Coverage vs reference CpG universe"),
    label  = c(paste0("Wilcoxon p=", signif(p_dataset, 3)),
               paste0("Wilcoxon p=", signif(p_ref, 3)))
  )

  p_overall <- ggplot(df_overall_long, aes(x = Group, y = Coverage_Pct)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(
      if (is.null(covariate_shape)) aes(color = Sample) else aes(color = Sample, shape = CovariateShape),
      width = 0.12, height = 0, size = 2, alpha = 0.8
    ) +
    facet_wrap(~ Metric, scales = "free_y") +
    scale_color_manual(values = sample_colors, breaks = names(sample_colors)) +
    { if (!is.null(covariate_shape)) scale_shape_discrete(name = covariate_shape) } +
    theme_bw() +
    labs(
      title = "Overall coverage breadth by group",
      x = "Group",
      y = "% CpGs Covered",
      color = "Sample"
    ) +
    geom_label(
      data = ann_overall,
      aes(x = 1.5, y = Inf, label = label),
      inherit.aes = FALSE,
      vjust = 1.2,
      size = 3.0,
      label.size = 0,
      alpha = 0.75
    )

  ggsave(file.path(output_dir, "Coverage_Breadth_Overall_Boxplot_ByGroup.pdf"),
         plot = p_overall, width = 10, height = 4.5)

  # ---- Save existing outputs ----
  ggsave(file.path(output_dir, "Coverage_Composition_Stacked_PerSample.pdf"),
         plot = p1, width = 12, height = 8)

  ggsave(file.path(output_dir, "Coverage_Composition_Stacked_PerSample_RefNCpGNorm.pdf"),
         plot = p1_ref_norm, width = 12, height = 8)

  ggsave(file.path(output_dir, "Coverage_Breadth_Boxplot_PerChr_DataNCpGNorm.pdf"),
         plot = p2, width = 12, height = 6)

  if (has_merged) {
    ggsave(file.path(output_dir, "Coverage_Breadth_Boxplot_PerChr_RefNCpGNorm_MStrand.pdf"),
          plot = p3_ref_breadth_merged, width = 12, height = 6)
  } else {
    ggsave(file.path(output_dir, "Coverage_Breadth_Boxplot_PerChr_RefNCpGNorm_PStrand.pdf"),
          plot = p3_ref_breadth_plus, width = 12, height = 6)

    ggsave(file.path(output_dir, "Coverage_Breadth_Boxplot_PerChr_RefNCpGNorm_NStrand.pdf"),
          plot = p3_ref_breadth_minus, width = 12, height = 6)
  }

  ggsave(file.path(output_dir, "Coverage_Breadth_Boxplot_ByGroup_PerChr.pdf"),
         plot = p2_by_group, width = 14, height = 9)

  write.csv(wilcox_by_chr,
            file = file.path(output_dir, "Coverage_Breadth_Wilcoxon_PerChr.csv"),
            row.names = FALSE)

  message("Plots saved to ", output_dir)

ref_tbl <- if (has_merged) df_breadth_ref_merged else df_breadth_ref_plus

invisible(list(
  composition = df_composition,
  composition_ref_norm = df_comp_ref_norm,
  breadth_dataset = df_breadth,
  breadth_reference = ref_tbl,
  breadth_reference_plus = if (!has_merged) df_breadth_ref_plus else NULL,
  breadth_reference_minus = if (!has_merged) df_breadth_ref_minus else NULL,
  breadth_by_group = df_breadth_g,
  wilcoxon_by_chr = wilcox_by_chr,
  overall = df_overall
))
}
