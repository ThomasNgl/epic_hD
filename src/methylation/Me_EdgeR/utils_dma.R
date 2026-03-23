library(readr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(edgeR)
library(SummarizedExperiment)
library(pheatmap)
library(gprofiler2)
library(stringr)
library(fdrtool)
library(dplyr)
library(tidyr)
#########################################################


##########################################
dge_prep <- function(se,
                     model,
                     assay_name = "counts",
                     norm_method = "TMM",
                     plot_bcv = TRUE,
                     disp_robust = FALSE) {

  counts <- SummarizedExperiment::assay(se, assay_name)
  if (!is.matrix(counts)) counts <- as.matrix(counts)

  # ---- If negative values: edgeR can't run. Return limma-friendly EList. ----
  if (any(counts < 0, na.rm = TRUE)) {
    message(
      "dge_prep(): assay '", assay_name,
      "' has negative values -> returning limma-style EList ",
      "(use limma::lmFit directly; no voom/edgeR)."
    )
    el <- list(
      E = counts,
      targets = as.data.frame(SummarizedExperiment::colData(se))
    )
    class(el) <- "EList"
    return(el)
  }

  # Use lib.size from colData if present (common in your objects)
  lib_size <- NULL
  cd <- SummarizedExperiment::colData(se)
  if ("lib.size" %in% colnames(cd)) {
    lib_size <- cd[["lib.size"]]
    if (length(lib_size) != ncol(counts)) lib_size <- NULL
  }

  dge <- edgeR::DGEList(counts = counts, lib.size = lib_size)

  if (!is.null(norm_method) && !identical(norm_method, "none")) {
    dge <- edgeR::calcNormFactors(dge, method = norm_method)
  }

  dge <- edgeR::estimateDisp(dge, design = model, robust = disp_robust)

  if (isTRUE(plot_bcv)) edgeR::plotBCV(dge)

  dge
}

##########################################
makeCorrectedAssay <- function(se, assay_name, fit, model, VarExp) {

  logCPM <- get_logcpm(se, assay_name)
  coefs <- fit$coefficients  # genes × coefficients
  
  if (is.null(coefs)) stop("fit has no coefficients. Run glmFit or glmQLFit first.")
  
  vars_in_model <- colnames(model)
  coef_names <- colnames(coefs)
  
  # Identify covariates to remove (exclude intercept and main variable)
  main_terms <- grep(VarExp, vars_in_model, value = TRUE)
  covars_to_remove <- setdiff(vars_in_model, c("(Intercept)", main_terms))
  
  if (length(covars_to_remove) == 0)
    stop("No covariates to correct for.")
  
  # Only keep terms present in both
  covars_to_remove <- intersect(covars_to_remove, coef_names)
  if (length(covars_to_remove) == 0)
    stop("None of the covariates to remove are in coefficient matrix.")
  
  cat("Removing covariates:", paste(covars_to_remove, collapse = ", "), "\n")
  
  # Compute unwanted component
  correction <- model[, covars_to_remove, drop = FALSE] %*%
    t(coefs[, covars_to_remove, drop = FALSE])
  
  # FIX: transpose correction to match logCPM orientation
  corrected_mat <- logCPM - t(correction)
  
  SummarizedExperiment::assay(se, "covcorrected") <- corrected_mat
  return(se)
}

##########################################
# Perform the dma
run_dma <- function(se, 
                    assay_name = 'counts',
                    row_name = 'gene_id',
                    VarExp = "Group", 
                    ExpGroup = "MSUS",
                    RefGroup = "Control",
                    nSV = 0, 
                    extra_formula = NULL,
                    scalefor = NULL,
                    CatVars = NULL,
                    VarNested = NULL,
                    disp_robust = FALSE,
                    plot_bcv = FALSE,
                    trend = TRUE,
                    engine = 'edgeR_lrt',
                    norm_method = 'TMM',
                    save_csv  = NULL
                  ) {
  
  # 🔹 Ensure grouping variable is a factor + set reference
  colData(se)[[VarExp]] <- factor(colData(se)[[VarExp]])
  colData(se)[[VarExp]] <- stats::relevel(colData(se)[[VarExp]], ref = RefGroup)
  
  # makes model 
  res_model <- model_prep(
    se = se,
    VarExp = VarExp,
    RefGroup = RefGroup,
    nSV = nSV,
    extra_formula = extra_formula,
    scalefor = scalefor,
    CatVars = CatVars,
    VarNested = VarNested
  )
  se <- res_model$se
  formula <- res_model$formula
  model <- res_model$model
  use_block <- res_model$use_block
  
  # Create DGEList
  dge <- dge_prep(se = se,
                  assay_name = assay_name,
                  model = model,
                  norm_method =norm_method,
                  plot_bcv = plot_bcv,
                  disp_robust = disp_robust
  )
  coef_name <- paste0(VarExp, ExpGroup)  # e.g., "GroupMSUS"
  
  # If coef was dropped by rank reduction, fail early with advice
  if (!coef_name %in% colnames(model)) {
    stop("The contrast '", coef_name, "' is not estimable with the current design. ",
         "If you intended to control for a nested factor, include it in extra_formula and set VarNested to its name.")
  }
  
  # ---------- Fit ----------
  if (use_block) {
    message("Use the limma voom technique because of nested variable")
    nested <- SummarizedExperiment::colData(se)[[VarNested]]
    v <- limma::voom(dge, model)
    corfit <- limma::duplicateCorrelation(v, model, block = nested)
    v <- limma::voom(dge, model, block = nested, correlation = corfit$consensus.correlation)
    fit <- limma::lmFit(v, model, block = nested, correlation = corfit$consensus.correlation)
    fit2 <- limma::eBayes(fit, trend = trend, robust = disp_robust)
    res <- limma::topTable(fit2, coef = coef_name, number = nrow(se), adjust.method = "BH")
    res <- res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    colnames(res) <- c("logFC", "avgExpr", "stat", "pval", "padj", "B")
  } else {
      if (grepl("limma", engine, ignore.case = TRUE)) {

        if (grepl("voom", engine, ignore.case = TRUE)) {
          v <- limma::voom(dge, design = model, plot = FALSE)
        } else {
          v <- dge
          attr(v, "design") <- model   # keep model attached even when skipping voom
        }

        fit <- limma::lmFit(v, design = model)
  
      fit2 <- limma::eBayes(fit, trend = trend, robust = disp_robust)
      res <- limma::topTable(
        fit2,
        coef          = coef_name,       # cd[[VarExp]]ExpGroup vs RefGroup
        number        = nrow(se),
        adjust.method = "BH"
      )
      res <- res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
      colnames(res) <- c("logFC", "avgExpr", "stat", "pval", "padj", "B")
    }
    if (engine == "edgeR_ql") {
      cat("Use the ql estimate")
      fit <- edgeR::glmQLFit(dge, model, robust = TRUE)
      qlf <- edgeR::glmQLFTest(fit, coef = coef_name)
      tab <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
      res <- data.frame(
        logFC   = tab$logFC,
        avgExpr = tab$logCPM,
        stat    = tab$F,
        pval    = tab$PValue,
        padj    = tab$FDR,
        B       = NA_real_
      )
      
    } else if (engine == "edgeR_lrt")  {
      fit <- edgeR::glmFit(dge, model)
      lrt <- edgeR::glmLRT(fit, coef = coef_name)
      tab <- lrt$table
      res <- data.frame(
        logFC   = tab$logFC,
        avgExpr = tab$logCPM,
        stat    = tab$LR,
        pval    = tab$PValue,
        padj    = p.adjust(tab$PValue, method = "BH"),
        B       = NA_real_
      )
    }
  }
  if (!is.null(row_name)) {
    rn <- as.character(SummarizedExperiment::rowData(se)[[row_name]])
    print(table(duplicated(rowData(se)[[row_name]])))

  } else {
    rn <- rownames(se)
  }
  
  # make duplicates unique: foo, foo_1, foo_2, ...
  rownames(res) <- make.unique(rn, sep = "_")
  
  formula_terms <- attr(terms(formula), "term.labels")
  makecorrected <- any(!formula_terms %in% VarExp)
  if(makecorrected){
  se <- makeCorrectedAssay(
    se     = se,
    assay_name = assay_name,
    fit    = fit,
    model  = model,
    VarExp = VarExp
  )}
  
  # Attach results + summary
  SummarizedExperiment::rowData(se)$dma <- res
  #if (!is.null(row_name)) {
  #  rownames(res) <- SummarizedExperiment::rowData(se)[[row_name]]
  #}
  ntranscripts <- nrow(res)
  nSignP <- sum(res$pval < 0.05, na.rm = TRUE)
  nSignFDR <- sum(res$padj < 0.05, na.rm = TRUE)
  cat(paste0(
    "\nDMA was performed using ", ntranscripts, " transcripts and using ", nSV, " SVs and engine " , engine, ".\n",
    nSignP, " transcripts had a p-value < 0.05. ",
    nSignFDR, " transcripts had an FDR-adjusted p-value < 0.05\n"
  ))
  res <- res[order(res$padj, res$pval, na.last = TRUE), ]
  if (!is.null(save_csv)) {
    dir.create(dirname(save_csv), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(res, save_csv, row.names = TRUE)
    message("Results written to: ", save_csv)
  }
  return(list(
    se = se,
    formula = formula
  ))}

##########################################
run_dma_2 <- function(se, 
                            assay_name = 'counts',
                            row_name = 'gene_id',
                            nSV, 
                            extra_formula = NULL,
                            VarExp = "Group", 
                            ExpGroup = "MSUS", 
                            RefGroup = "Control",
                            VarNested = NULL,
                            scalefor = NULL,
                            engine = 'edgeR_lrt',
                            norm_method = 'TMM',
                            CatVars = c("Group"),
                            disp_robust = TRUE,
                            trend = TRUE,
                            save_csv = NULL
) {
  
  # Run DmA
  dma_res <- run_dma(
    se           = se,
    assay_name = assay_name,
    row_name = row_name,
    VarExp       = VarExp,
    ExpGroup     = ExpGroup,
    RefGroup     = RefGroup,
    nSV          = nSV,
    extra_formula = extra_formula,
    CatVars      = CatVars,
    scalefor = scalefor,
    VarNested = VarNested,
    disp_robust  = disp_robust,
    trend = trend,
    plot_bcv     = TRUE,
    engine = engine,
    save_csv = save_csv
  )
  
  # Extract SummarizedExperiment
  se_out <- dma_res$se
  print(rowData(se_out)$dma)

  # Choose statistic name depending on QL
  feature_name <- if (engine == "edgeR_ql") {
    "F"
  } else if (engine == "edgeR_lrt") {
    "LR"
  } else if (engine == "limma") {
    "t"
  } else {
    stop("Unknown engine: ", engine)
  }
  
  # Extract p-value data frame
  dma_tbl <- as.data.frame(SummarizedExperiment::rowData(se_out)$dma)
  
    method_label <- if (engine == "edgeR_ql") {
      "QL"
    } else if (engine == "edgeR_lrt") {
      "LRT"
    } else if (engine == "limma") {
      "limma"
    } else {
      stop("Unknown engine: ", engine)
    }
    
  df <- tibble(
    pval           = dma_tbl$pval,
    pval_bh        = dma_tbl$padj,
    Method         = method_label
  )
  
  list(se = se_out, df = df)
}

#########################
# ───────────────────────────────
# Generic helper functions
# ───────────────────────────────

# 1️⃣ Sample annotation builder (general SE utility)
makeAnnotation <- function(se, VarExp, Covariates = NULL) {
  # Combine main and additional covariates, ensuring uniqueness
  vars <- c(setdiff(Covariates, VarExp), VarExp)  
  anno <- as.data.frame(SummarizedExperiment::colData(se)[, vars, drop = FALSE])
  
  # Order samples for consistent layout
  anno$order_num <- suppressWarnings(as.numeric(gsub("[[:alpha:]]", "", rownames(anno))))
  anno <- anno[order(anno[[VarExp]], anno$order_num), , drop = FALSE]
  anno$order_num <- NULL
  
  return(anno)
}

# 3️⃣ Heatmap drawer (used for any matrix)
# at top of your file (for ComplexHeatmap numeric scales)
if (!requireNamespace("circlize", quietly = TRUE)) {
  stop("Package 'circlize' is required for ComplexHeatmap numeric annotations.")
}
makeHeatmap <- function(
    mat,
    anno,
    VarExp,
    RefGroup,
    ExpGroup,
    title,
    save_path = NULL,
    width = 15,
    height = 10,
    show_rownames = TRUE,
    fontsize_row = 10,
    hide_row_tree = TRUE,
    engine = c("pheatmap","ComplexHeatmap"),
    show_colnames = TRUE,
    fontsize_colnames = 8,
    colnames_rot = 90,
    row_hclust = NULL               # NEW: precomputed row clustering
) {
  engine <- match.arg(engine)
  mat <- mat[, rownames(anno), drop = FALSE]
  
  # --- Build annotation color maps (as before) ---
  ann_colors_pheatmap <- list()
  ann_colors_ch       <- list()
  
  palette_names <- c("Set1", "Dark2", "Accent", "Pastel1", "Paired", "Set3")
  palette_index <- 1
  gradient_list <- list(
    c("navy","white","firebrick3"),
    c("darkgreen","white","violet"),
    c("darkorange","white","purple4"),
    c("black","white","darkred")
  )
  
  for (var in colnames(anno)) {
    vals <- anno[[var]]
    
    is_cont <- is.numeric(vals) || suppressWarnings(!any(is.na(as.numeric(as.character(vals)))))
    
    if (is_cont) {
      num_vals <- suppressWarnings(as.numeric(as.character(vals)))
      anno[[var]] <- num_vals
      
      rng <- range(num_vals, na.rm = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
        rng <- c(num_vals[is.finite(num_vals)][1] - 1, num_vals[is.finite(num_vals)][1] + 1)
      }
      mid <- stats::median(num_vals, na.rm = TRUE)
      gpal <- gradient_list[[((palette_index - 1) %% length(gradient_list)) + 1]]
      ann_colors_ch[[var]] <- circlize::colorRamp2(c(rng[1], mid, rng[2]), gpal)
      palette_index <- palette_index + 1
      
    } else {
      fac <- as.factor(vals)
      anno[[var]] <- fac
      levs <- levels(fac)
      
      if (var == VarExp) {
        pal_vec <- setNames(c("#989898", "#be0b34"), c(RefGroup, ExpGroup))
      } else {
        pal_name <- palette_names[(palette_index - 1) %% length(palette_names) + 1]
        if (length(levs) >= 3) {
          max_n <- RColorBrewer::brewer.pal.info[pal_name, "max"]
          ncols <- min(length(levs), max_n)
          pal   <- RColorBrewer::brewer.pal(ncols, pal_name)[seq_along(levs)]
        } else if (length(levs) == 2) {
          pal <- RColorBrewer::brewer.pal(3, pal_name)[1:2]
        } else {
          pal <- "#aaaaaa"
        }
        names(pal) <- levs
        pal_vec <- pal
        palette_index <- palette_index + 1
      }
      
      ann_colors_pheatmap[[var]] <- pal_vec
      ann_colors_ch[[var]]       <- pal_vec
    }
  }
  
  # If row_hclust is provided, use it; otherwise let the engine cluster
  cluster_rows_arg <- if (is.null(row_hclust)) TRUE else row_hclust
  
  # --- pheatmap path ---
  if (engine == "pheatmap") {
    hm <- pheatmap::pheatmap(
      mat,
      annotation_col    = anno,
      annotation_colors = ann_colors_pheatmap,
      cluster_cols      = FALSE,
      cluster_rows      = cluster_rows_arg,
      treeheight_row    = if (hide_row_tree) 0 else 50,
      show_rownames     = show_rownames,
      show_colnames     = show_colnames,
      fontsize_row      = fontsize_row,
      fontsize_col      = fontsize_colnames,
      color             = colorRampPalette(c("navy","white","firebrick3"))(50),
      scale             = "row",
      main              = title,
      fontsize          = 15,
      border_color      = NA,
      silent            = FALSE,
      angle_col         = colnames_rot
    )
    if (!is.null(save_path)) {
      dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
      pdf(save_path, width = width, height = height)
      print(hm)
      dev.off()
    } else {
      print(hm)
    }
    return(invisible(hm))
  }
  
  # --- ComplexHeatmap path ---
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df   = anno,
    col  = ann_colors_ch,
    which = "column",
    annotation_name_side = "left"
  )
  
  # z-score by row like pheatmap(scale="row")
  zfun <- function(m) { m2 <- t(scale(t(m))); m2[is.nan(m2)] <- 0; m2 }
  zmat <- zfun(as.matrix(mat))
  
  ht <- ComplexHeatmap::Heatmap(
    zmat,
    name = "Z-score",
    col = colorRampPalette(c("navy","white","firebrick3"))(50),
    cluster_rows = cluster_rows_arg,
    cluster_columns = FALSE,
    show_row_names = show_rownames,
    row_names_gp   = grid::gpar(fontsize = fontsize_row),
    show_column_names = show_colnames,
    column_names_side = "bottom",
    column_names_rot  = colnames_rot,
    column_names_gp   = grid::gpar(fontsize = fontsize_colnames),
    top_annotation    = ha,
    column_title      = title,
    border            = FALSE,
    show_row_dend     = !hide_row_tree       # hide or show dendrogram itself
  )
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(save_path, width = width, height = height)
    ComplexHeatmap::draw(ht, heatmap_legend_side = "left", annotation_legend_side = "bottom")
    grDevices::dev.off()
  } else {
    ComplexHeatmap::draw(ht, heatmap_legend_side = "left", annotation_legend_side = "bottom")
  }
  invisible(ht)
}

# ───────────────────────────────
# Heatmap-specific core function
# ───────────────────────────────
plot_heatmap <- function(
    se,
    assay_name      = 'counts',
    log_assay_name  = 'logcpm',
    outputdir       = NULL,
    RefGroup        = "Control",
    ExpGroup        = "MSUS",
    VarExp          = "Group",
    Covariates      = c("Group"),
    fdr_threshold   = 0.05,
    pval_threshold  = 0.05,
    sva_corrected   = FALSE,
    cov_corrected   = FALSE,
    convert_ids     = TRUE,
    title           = "Heatmap",
    fontsize_colnames = 10,
    fontsize_row      = 10,
    engine          = c("pheatmap"),
    result_name     = "dma",          # which result in rowData(se)
    label_column    = "gene_symbol"   # which column in that result to use as labels
) {
  engine <- match.arg(engine)
  
  # ---------------------------
  # Sample annotation
  # ---------------------------
  anno <- makeAnnotation(se, VarExp, Covariates)
  
  # ---------------------------
  # Extract result data.frame from rowData(se)
  # ---------------------------
  res_df_all <- SummarizedExperiment::rowData(se)[[result_name]]
  if (is.null(res_df_all)) {
    stop("rowData(se)[[", result_name, "]] is NULL; did you attach results as '", result_name, "'?")
  }
  
  # rownames(res_df_all) are technical IDs (feature_id / assay rownames)
  feature_ids <- rownames(res_df_all)
  
  # master DE table with explicit feature_id + label column
  label_vec <- if (label_column %in% names(res_df_all)) as.character(res_df_all[[label_column]]) else rep(NA_character_, length(feature_ids))
  label_vec[is.na(label_vec) | label_vec == ""] <- feature_ids

  res_heat <- data.frame(
    feature_id  = feature_ids,
    label       = label_vec,
    logFC       = res_df_all$logFC,
    avgExp      = res_df_all$avgExpr,
    stat        = res_df_all$stat,
    pval        = res_df_all$pval,
    padj        = res_df_all$padj,
    stringsAsFactors = FALSE,
    row.names   = feature_ids
  )
  # Order by p-value so "top" means smallest p
  res_heat <- res_heat[order(res_heat$pval), , drop = FALSE]
  
  # Base assay used for hierarchical clustering
  base_assay_for_clustering <- if (sva_corrected) "corrected" else log_assay_name
  
  # ---------------------------
  # Helpers
  # ---------------------------
  
  # Numeric matrix from an assay, indexed by technical IDs, in the given order
  get_assay_matrix <- function(feature_ids_subset, assay_name) {
    mat_all <- SummarizedExperiment::assay(se, assay_name)
    if (!all(feature_ids_subset %in% rownames(mat_all))) {
      missing <- setdiff(feature_ids_subset, rownames(mat_all))
      if (length(missing) > 0L) {
        warning("Assay '", assay_name, "' is missing ", length(missing),
                " features; they will be dropped.")
      }
      feature_ids_subset <- intersect(feature_ids_subset, rownames(mat_all))
    }
    mat_all[feature_ids_subset, , drop = FALSE]
  }
  
  # Matrix for plotting: same data as assay, but rownames = labels (gene symbols, etc.)
  build_matrix_with_labels <- function(feature_ids_subset, assay_name, res_heat) {
    if (length(feature_ids_subset) == 0L) {
      return(NULL)
    }
    mat <- get_assay_matrix(feature_ids_subset, assay_name)
    if (nrow(mat) == 0L) return(NULL)
    
    # derive labels from res_heat$label
    labels <- res_heat$label[match(rownames(mat), res_heat$feature_id)]
    labels[is.na(labels) | labels == ""] <- rownames(mat)  # fallback to IDs
    rownames(mat) <- labels
    mat
  }
  
  # Emit all requested heatmaps for a given ordered feature set
  emitHeatmaps <- function(feature_ids_ordered,
                           filter_label,
                           suffix_tag = "",
                           show_rownames,
                           engine) {
    
    # 1) SVA-corrected
    if (sva_corrected) {
      mat_corr <- build_matrix_with_labels(feature_ids_ordered, "corrected", res_heat)
      if (!is.null(mat_corr) && nrow(mat_corr) >= 2) {
        makeHeatmap(
          mat_corr,
          anno, VarExp, RefGroup, ExpGroup,
          fontsize_colnames = fontsize_colnames,
          fontsize_row      = fontsize_row,
          title             = sprintf("%s — %s%s (SVA corrected)", title, filter_label, suffix_tag),
          save_path         = if (!is.null(outputdir)) file.path(outputdir, sprintf("Heatmap%s%s.sva_corrected.pdf", filter_label, suffix_tag)) else NULL,
          show_rownames     = show_rownames,
          engine            = engine,
          row_hclust        = FALSE,   # <- already ordered; do not re-cluster
          hide_row_tree     = TRUE
        )
      } else {
        message("Skipping SVA heatmap for ", filter_label, suffix_tag, ": < 2 rows.")
      }
    }
    
    # 2) Covariate-corrected
    if (cov_corrected) {
      mat_cov <- build_matrix_with_labels(feature_ids_ordered, "covcorrected", res_heat)
      if (!is.null(mat_cov) && nrow(mat_cov) >= 2) {
        makeHeatmap(
          mat_cov,
          anno, VarExp, RefGroup, ExpGroup,
          fontsize_colnames = fontsize_colnames,
          fontsize_row      = fontsize_row,
          title             = sprintf("%s — %s%s (cov corrected)", title, filter_label, suffix_tag),
          save_path         = if (!is.null(outputdir)) file.path(outputdir, sprintf("Heatmap%s%s.cov_corrected.pdf", filter_label, suffix_tag)) else NULL,
          show_rownames     = show_rownames,
          engine            = engine,
          row_hclust        = FALSE,
          hide_row_tree     = TRUE
        )
      } else {
        message("Skipping cov_corrected heatmap for ", filter_label, suffix_tag, ": < 2 rows.")
      }
    }
    
    # 3) logCPM (or whatever log_assay_name is)
    mat_log <- build_matrix_with_labels(feature_ids_ordered, log_assay_name, res_heat)
    if (!is.null(mat_log) && nrow(mat_log) >= 2) {
      makeHeatmap(
        mat_log,
        anno, VarExp, RefGroup, ExpGroup,
        fontsize_colnames = fontsize_colnames,
        fontsize_row      = fontsize_row,
        title             = sprintf("%s — %s%s (%s)", title, filter_label, suffix_tag, log_assay_name),
        save_path         = if (!is.null(outputdir)) file.path(outputdir, sprintf("Heatmap%s%s.%s.pdf", filter_label, suffix_tag, log_assay_name)) else NULL,
        show_rownames     = show_rownames,
        engine            = engine,
        row_hclust        = FALSE,
        hide_row_tree     = TRUE
      )
    } else {
      message("Skipping ", log_assay_name, " heatmap for ", filter_label, suffix_tag, ": < 2 rows.")
    }
  }
  
  # ---------------------------
  # Core per-filter logic
  # ---------------------------
  makeFilteredHeatmap <- function(filter_col,
                                  threshold,
                                  filter_label,
                                  show_rownames_full,
                                  engine,
                                  top20 = FALSE) {
    
    # select features by padj / pval, using technical IDs
    sel       <- res_heat[[filter_col]] < threshold
    res_sub   <- res_heat[sel, , drop = FALSE]
    feat_all  <- res_sub$feature_id
    
    if (length(feat_all) < 2) {
      message("No heatmap for ", filter_label, ": only ", length(feat_all),
              " features pass ", filter_col, " < ", threshold)
      return(NULL)
    }
    
    # --- split into up/down (and neutral) by logFC ---
    up_idx      <- which(res_sub$logFC > 0)
    down_idx    <- which(res_sub$logFC < 0)
    neutral_idx <- setdiff(seq_len(nrow(res_sub)), c(up_idx, down_idx))
    
    feat_up      <- res_sub$feature_id[up_idx]
    feat_down    <- res_sub$feature_id[down_idx]
    feat_neutral <- res_sub$feature_id[neutral_idx]  # rarely used, but we keep it
    
    # helper to cluster a subset with z-normalisation
    cluster_subset <- function(feat_subset) {
      if (length(feat_subset) < 2) return(character(0))
      base_mat <- get_assay_matrix(feat_subset, base_assay_for_clustering)
      if (nrow(base_mat) < 2) return(character(0))
      # z-score rows for clustering (to match the visual scale)
      zmat <- t(scale(t(as.matrix(base_mat))))
      zmat[is.nan(zmat)] <- 0
      hc <- stats::hclust(stats::dist(zmat))
      feat_subset[hc$order]
    }
    
    order_up      <- cluster_subset(feat_up)
    order_neutral <- cluster_subset(feat_neutral)
    order_down    <- cluster_subset(feat_down)
    
    # concatenated order: all MSUS-up at top, then neutral, then down
    feat_ordered_all <- c(order_up, order_neutral, order_down)
    
    if (length(feat_ordered_all) < 2) {
      message("Not enough rows after clustering for ", filter_label)
      return(NULL)
    }
    
    # 1) Full heatmaps in this order
    emitHeatmaps(
      feature_ids_ordered = feat_ordered_all,
      filter_label        = filter_label,
      suffix_tag          = "",
      show_rownames       = show_rownames_full,
      engine              = engine
    )
    
    # 2) Optional "top20 / duplicate" behavior
    if (top20) {
      n_sub <- nrow(res_sub)
      n_top <- min(20L, n_sub)
      # top 20 by p-value (res_heat is already ordered by pval)
      feat_top_by_pval <- res_sub$feature_id[seq_len(n_top)]
      # keep them in the same up→neutral→down clustered order
      feat_top_ordered <- intersect(feat_ordered_all, feat_top_by_pval)
      
      if (length(feat_top_ordered) >= 2) {
        emitHeatmaps(
          feature_ids_ordered = feat_top_ordered,
          filter_label        = filter_label,
          suffix_tag          = ".top20",
          show_rownames       = TRUE,     # always show names for top 20
          engine              = engine
        )
      } else {
        # ≤ 20 genes overall: re-plot full set with rownames
        emitHeatmaps(
          feature_ids_ordered = feat_ordered_all,
          filter_label        = filter_label,
          suffix_tag          = ".withnames",
          show_rownames       = TRUE,
          engine              = engine
        )
      }
    }
  }
  
  # ---------------------------
  # Run both filters
  # ---------------------------
  # FDR: typically show rownames
  makeFilteredHeatmap(
    filter_col         = "padj",
    threshold          = fdr_threshold,
    filter_label       = "FDR",
    show_rownames_full = TRUE,
    engine             = engine,
    top20              = FALSE
  )
  
  # P-value: full heatmap without names, plus top20 logic
  makeFilteredHeatmap(
    filter_col         = "pval",
    threshold          = pval_threshold,
    filter_label       = "Pval",
    show_rownames_full = FALSE,
    engine             = engine,
    top20              = TRUE
  )
  
  return(res_heat)
}



##########################################
#MA matrix
plot_ma <- function(
    se,
    outputdir      = NULL,
    title          = "MA Plot",
    fdrtool        = FALSE,
    padj_threshold = 0.05
) {
  # Expression (A) from logCPM
  assay_names   <- SummarizedExperiment::assayNames(se)
  
  has_counts   <- "counts" %in% assay_names
  if (has_counts){
  counts <- SummarizedExperiment::assay(se, "counts")
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  }
  else{logcpm <- SummarizedExperiment::assay(se, "expr")}
  A <- rowMeans(logcpm)  # average expression per gene
  
  # Pull DE results
  rd <- SummarizedExperiment::rowData(se)
  
  has_fdrtool <- !is.null(rd$fdrtool)
  has_dma     <- !is.null(rd$dma)
  
  if (fdrtool && has_fdrtool) {
    res <- rd$fdrtool
    M    <- res$logFC
    padj <- res$qval_fdrtool
    src  <- "fdrtool"
  } else if (has_dma) {
    res <- rd$dma
    M    <- res$logFC
    padj <- res$padj
    src  <- "DMA"
    if (fdrtool && !has_fdrtool) {
      message("fdrtool results not found; falling back to 'dma'.")
    }
  } else {
    stop("No DE results found in rowData(se)$fdrtool or rowData(se)$dma.")
  }
  
  # Align lengths defensively (should already match)
  n <- min(length(A), length(M), length(padj))
  A    <- A[seq_len(n)]
  M    <- M[seq_len(n)]
  padj <- padj[seq_len(n)]
  
  sig <- is.finite(padj) & (padj < padj_threshold)
  
  df <- data.frame(
    A = A,
    M = M,
    group = factor(ifelse(sig, "FDR < 0.05", "Not significant"),
                   levels = c("FDR < 0.05", "Not significant"))
  )
  
  # Draw non-sig first, then sig on top
  df_sig  <- df[df$group == "FDR < 0.05", , drop = FALSE]
  df_nsig <- df[df$group == "Not significant", , drop = FALSE]
  
  subtitle_src <- if (src == "fdrtool") "Significance by fdrtool q-values" else "Significance by DMA padj"
  legend_name  <- NULL
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = A, y = M)) +
    ggplot2::geom_point(data = df_nsig, ggplot2::aes(color = group), size = 1.2, alpha = 0.6, na.rm = TRUE) +
    ggplot2::geom_point(data = df_sig,  ggplot2::aes(color = group), size = 1.6, alpha = 0.9, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_color_manual(
      values = c("FDR < 0.05" = "red", "Not significant" = "black"),
      breaks = c("FDR < 0.05", "Not significant"),
      name = legend_name
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle_src,
      x = "Average expression (logCPM)",
      y = "Log2 fold change"
    ) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "top")
  
  if (!is.null(outputdir)) {
    dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(file.path(outputdir, "MAPlot.pdf"))
    print(p)
    grDevices::dev.off()
  } else {
    print(p)
  }
  
  invisible(p)
}


###################################
# GO analysis wrapper (revised)
run_go_analysis <- function(se,
                            outputdir = NULL,
                            organism = "mmusculus",
                            sources = "GO:BP",
                            fdr_threshold = 0.05,
                            pval_threshold = 0.01,
                            logfc_fdr = 1,
                            logfc_pval = 2,
                            fdrtool = FALSE,
                            logic = c("and", "or"),
                            save_path = NULL) {
  
  logic <- match.arg(logic)
  
  # -----------------------------
  # Pull DE results (DMA or fdrtool)
  # -----------------------------
  rd <- SummarizedExperiment::rowData(se)
  
  has_fdrtool <- !is.null(rd$fdrtool)
  has_dma     <- !is.null(rd$dma)
  
  if (fdrtool && has_fdrtool) {
    res <- as.data.frame(rd$fdrtool)
    padj_col <- "qval_fdrtool"  # adjusted p-values from fdrtool
    pval_col <- "pval"          # raw p-values in fdrtool table
    src      <- "fdrtool"
  } else if (has_dma) {
    res <- as.data.frame(rd$dma)
    padj_col <- "padj"          # adjusted p-values from DMA
    pval_col <- "pval"          # raw p-values from DMA
    src      <- "DMA"
    if (fdrtool && !has_fdrtool) {
      message("fdrtool results not found; falling back to 'dma'.")
    }
  } else {
    stop("No DE results found in rowData(se)$fdrtool or rowData(se)$dma.")
  }
  
  if (!("logFC" %in% colnames(res))) {
    stop("Column 'logFC' not found in the selected DE result (", src, ").")
  }
  if (!(padj_col %in% colnames(res))) {
    stop("Adjusted p-value column '", padj_col, "' not found in the selected DE result (", src, ").")
  }
  if (!(pval_col %in% colnames(res))) {
    stop("P-value column '", pval_col, "' not found in the selected DE result (", src, ").")
  }
  
  # -----------------------------
  # Detect significant genes
  # -----------------------------
  nSignFDR <- sum(res[[padj_col]] < fdr_threshold, na.rm = TRUE)
  
  fdr_sig <- res[[padj_col]] < fdr_threshold
  p_sig   <- res[[pval_col]] < pval_threshold
  
  # --- Select transcripts ---
  if (nSignFDR > 5) {
    # Use FDR thresholds
    if (logic == "and") {
      transcripts_up   <- rownames(res[fdr_sig & res$logFC >  logfc_fdr, , drop = FALSE])
      transcripts_down <- rownames(res[fdr_sig & res$logFC < -logfc_fdr, , drop = FALSE])
    } else { # logic == "or"
      up_mask   <- (fdr_sig & res$logFC > 0)  | (res$logFC >  logfc_fdr)
      down_mask <- (fdr_sig & res$logFC < 0)  | (res$logFC < -logfc_fdr)
      transcripts_up   <- rownames(res[up_mask,   , drop = FALSE])
      transcripts_down <- rownames(res[down_mask, , drop = FALSE])
    }
  } else {
    # Use p-value thresholds
    if (logic == "and") {
      transcripts_up   <- rownames(res[p_sig & res$logFC >  logfc_pval, , drop = FALSE])
      transcripts_down <- rownames(res[p_sig & res$logFC < -logfc_pval, , drop = FALSE])
    } else { # logic == "or"
      up_mask   <- (p_sig & res$logFC > 0)  | (res$logFC >  logfc_pval)
      down_mask <- (p_sig & res$logFC < 0)  | (res$logFC < -logfc_pval)
      transcripts_up   <- rownames(res[up_mask,   , drop = FALSE])
      transcripts_down <- rownames(res[down_mask, , drop = FALSE])
    }
  }
  
  # -----------------------------
  # Helper to run g:Profiler
  # -----------------------------
  run_gost <- function(ids) {
    if (length(ids) == 0) return(NULL)
    gprofiler2::gost(
      query = ids,
      organism = organism,
      ordered_query = FALSE,
      multi_query = FALSE,
      significant = TRUE,
      exclude_iea = FALSE,
      measure_underrepresentation = FALSE,
      evcodes = TRUE,
      user_threshold = fdr_threshold,
      correction_method = "fdr",
      domain_scope = "annotated",
      custom_bg = NULL,
      numeric_ns = "",
      sources = sources,
      as_short_link = FALSE,
      highlight = TRUE
    )
  }
  
  # -----------------------------
  # Run enrichment
  # -----------------------------
  genes_up_symb      <- convertGeneNames(transcripts_up,   organism = organism)
  genes_down_symb    <- convertGeneNames(transcripts_down, organism = organism)
  transcripts_merged <- unique(c(transcripts_up, transcripts_down))
  genes_merged_symb  <- unique(c(genes_up_symb, genes_down_symb))
  
  go_up     <- run_gost(genes_up_symb)
  go_down   <- run_gost(genes_down_symb)
  go_merged <- run_gost(genes_merged_symb)
  print("unrec genes")
  print(go_merged$query_metadata$genes_without_hits)
  print("voila")
  
  # Save into SE metadata
  S4Vectors::metadata(se)$go_up     <- go_up
  S4Vectors::metadata(se)$go_down   <- go_down
  S4Vectors::metadata(se)$go_merged <- go_merged
  
  # -----------------------------
  # Helper to plot & save GO tables
  # -----------------------------
  plot_and_save <- function(go_res, title, filename, save_path = NULL, outputdir = NULL, max_terms = 50) {
    if (is.null(go_res) || is.null(go_res$result) || nrow(go_res$result) == 0) {
      message(title, ": no enriched terms found.")
      return(NULL)
    }
    
    ord <- order(go_res$result$p_value, na.last = NA)
    res_sorted  <- go_res$result[ord, , drop = FALSE]
    res_limited <- utils::head(res_sorted, max_terms)
    
    # Console preview
    print(utils::head(res_sorted[, c("term_name","p_value","term_size","intersection_size")], 10))
    if (nrow(res_sorted) > max_terms) {
      message("Plot limited to top ", max_terms, " terms. Full table still available.")
    }
    
    # Prepare trimmed object for plotting
    go_res_plot <- go_res
    go_res_plot$result <- res_limited
    
    gg <- gprofiler2::publish_gosttable(
      go_res_plot,
      use_colors   = TRUE,
      show_columns = c("term_name","p_value","term_size","intersection_size"),
      filename     = NULL,
      ggplot       = TRUE
    ) + ggplot2::ggtitle(title)
    
    # Compute the path (outputdir OR explicit save_path)
    final_path <- NULL
    if (!is.null(save_path)) {
      final_path <- save_path
    } else if (!is.null(outputdir)) {
      final_path <- file.path(outputdir, filename)
    }
    
    # Save via pdf()/print()/dev.off()
    if (!is.null(final_path)) {
      dir.create(dirname(final_path), showWarnings = FALSE, recursive = TRUE)
      fig_height <- max(5, min(20, nrow(res_limited) / 3 + 2))
      grDevices::pdf(final_path, width = 12, height = fig_height)
      print(gg)
      grDevices::dev.off()
      
      # Flatten list columns before writing CSV
      flat_res_sorted <- res_sorted
      flat_res_sorted[] <- lapply(flat_res_sorted, function(x) {
        if (is.list(x)) vapply(x, function(el)
          paste(el, collapse = ","), character(1L)) else x
      })
      
      utils::write.csv(
        flat_res_sorted,
        file = sub("\\.pdf$", ".csv", final_path),
        row.names = FALSE
      )
    } else {
      print(gg)
    }
    
    list(result = res_sorted, plot = gg)
  }
  
  # -----------------------------
  # Make GO plots + tables
  # -----------------------------
  res_up <- plot_and_save(
    go_res   = go_up,
    title    = "GO analysis of upregulated transcripts",
    filename = "GO.up.pdf",
    save_path = if (!is.null(outputdir)) file.path(outputdir, "GO.up.pdf") else NULL,
    outputdir = NULL
  )
  
  res_down <- plot_and_save(
    go_res   = go_down,
    title    = "GO analysis of downregulated transcripts",
    filename = "GO.down.pdf",
    save_path = if (!is.null(outputdir)) file.path(outputdir, "GO.down.pdf") else NULL,
    outputdir = NULL
  )
  
  res_merged <- plot_and_save(
    go_res   = go_merged,
    title    = "GO analysis of merged transcripts",
    filename = "GO.merged.pdf",
    save_path = if (!is.null(outputdir)) file.path(outputdir, "GO.merged.pdf") else NULL,
    outputdir = NULL
  )
  
  # -----------------------------
  # Save DEG list as CSV with pval, padj, logFC
  # -----------------------------
  if (!is.null(save_path)) {
    # Ensure .csv extension
    if (!grepl("\\.csv$", save_path, ignore.case = TRUE)) {
      save_path <- paste0(save_path, ".csv")
    }
    
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    
    # Subset DE table for merged transcripts
    deg_df <- res[transcripts_merged, c(pval_col, padj_col, "logFC"), drop = FALSE]
    
    # Standardize column names to pval / padj / logFC in the output
    colnames(deg_df) <- c("pval", "padj", "logFC")
    
    # Build output table: transcript ID, gene symbol, pval, padj, logFC
    out_df <- data.frame(
      deg_id = rownames(deg_df),
      deg    = convertGeneNames(rownames(deg_df), organism = organism),
      pval   = deg_df$pval,
      padj   = deg_df$padj,
      logFC  = deg_df$logFC,
      row.names = NULL
    )
    
    # Order by padj (increasing) and format numeric columns to 3 decimals
    out_df <- out_df[order(out_df$padj, decreasing = FALSE), , drop = FALSE]
    num_cols <- c("pval", "padj", "logFC")
    out_df[num_cols] <- lapply(out_df[num_cols], function(x) sprintf("%.3f", as.numeric(x)))
    
    utils::write.csv(
      out_df,
      file = save_path,
      row.names = FALSE
    )
  }
  
  invisible(list(up = res_up, down = res_down, merged = res_merged))
  return(list(deg_id = transcripts_merged, deg = genes_merged_symb))
}

###################
# Summary
get_gene_ranges_table <- function(
    deg_id_list,
    genes_bed = "/mnt/thomas/matriline/results/GRCm39/annotations/mm39_anno.genes.primary.bed",
    chrom_sizes = "/mnt/thomas/matriline/results/GRCm39/annotations/mm39_anno.chrom.sizes.primary",
    output_tsv = NULL
) {
  # Vérifications de base
  stopifnot(is.character(deg_id_list), length(deg_id_list) > 0)
  if (!file.exists(genes_bed)) stop("Genes BED not found: ", genes_bed)
  
  # 1️⃣ Conversion ENSMUSG -> gene symbol (tu as déjà convertGeneNames)
  gene_symbols <- convertGeneNames(deg_id_list, convert_ids = TRUE, organism = "mmusculus")
  
  # 2️⃣ Lecture du fichier BED d’annotation (colonnes: chr, start, end, symbol, score, strand)
  g <- utils::read.table(genes_bed, sep = "\t", header = FALSE, quote = "",
                         comment.char = "", stringsAsFactors = FALSE)
  colnames(g)[1:6] <- c("chr","start","end","gene_symbol","score","strand")
  
  # 3️⃣ Filtrage sur les symboles correspondants
  keep <- g$gene_symbol %in% unique(gene_symbols)
  g_sel <- g[keep, c("chr","start","end","strand","gene_symbol"), drop = FALSE]
  
  # Si aucun match
  if (!nrow(g_sel)) {
    warning("Aucun gène correspondant trouvé dans le fichier BED.")
    return(g_sel)
  }
  
  # 4️⃣ Calcul de la taille du gène
  g_sel$size <- g_sel$end - g_sel$start
  
  # 5️⃣ Ajout de la colonne gene_id (à partir de ta liste originale)
  symbol_to_id <- setNames(deg_id_list, gene_symbols)
  g_sel$gene_id <- symbol_to_id[g_sel$gene_symbol]
  
  # 6️⃣ Tri par chromosome et position de départ
  if (!is.null(chrom_sizes) && file.exists(chrom_sizes)) {
    chr_levels <- read.table(chrom_sizes, sep = "\t", header = FALSE, stringsAsFactors = FALSE)[[1]]
  } else {
    chr_levels <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM", "chrMT")
  }
  chr_rank <- match(g_sel$chr, chr_levels)
  chr_rank[is.na(chr_rank)] <- max(chr_rank, na.rm = TRUE) + as.integer(factor(g_sel$chr[is.na(chr_rank)]))
  g_sel <- g_sel[order(chr_rank, g_sel$start), , drop = FALSE]
  
  # 7️⃣ Réorganisation des colonnes finales
  out <- g_sel[, c("chr","start","end","strand","size","gene_id","gene_symbol")]
  
  # 8️⃣ Écriture optionnelle du TSV
  if (!is.null(output_tsv)) {
    utils::write.table(out, file = output_tsv, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = TRUE)
    message("Saved: ", output_tsv)
  }
  
  return(out)
}


##################################################################
# Take a csv file containing a column of feature and a column of log FC and pval 
# Returns a subaset of the csv 
filter_edgeR_results <- function(
    file,
    padj_cutoff = NULL,
    pval_cutoff = NULL,
    logFC_min   = NULL,
    logFC_max   = NULL,
    connector   = c("AND", "OR")
) {
  connector <- match.arg(connector)
  
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  if (names(df)[1] == "" || is.na(names(df)[1])) {
    names(df)[1] <- "Top.annotation.name"
  }
  conds <- list()
  
  if (!is.null(padj_cutoff) && "padj" %in% names(df)) {
    conds$padj <- df$padj < padj_cutoff
  }
  
  if (!is.null(pval_cutoff) && "pval" %in% names(df)) {
    conds$pval <- df$pval < pval_cutoff
  }
  
  if (!is.null(logFC_min) && "logFC" %in% names(df)) {
    conds$logFC_min <- df$logFC >= logFC_min
  }
  
  if (!is.null(logFC_max) && "logFC" %in% names(df)) {
    conds$logFC_max <- df$logFC <= logFC_max
  }
  
  if (length(conds) == 0L) {
    warning("No filters provided (all arguments NULL). Returning full data frame.")
    return(df)
  }
  
  cond_mat <- do.call(cbind, conds)
  
  keep <- if (connector == "AND") {
    apply(cond_mat, 1, all)
  } else {
    apply(cond_mat, 1, any)
  }
  
  df[keep, , drop = FALSE]
}

