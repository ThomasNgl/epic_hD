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

#######################
# Define the full and the null model
set_model_formulas <- function(
    se,
    VarExp,
    CovarBio = NULL,
    CovarTec = NULL
) {
  # Helper to create formula from vector of terms
  create_formula <- function(terms) {
    if (length(terms) == 0) return(stats::as.formula("~ 1"))
    stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  }
  
  # Helper to generate all non-empty subsets of a character vector
  all_subsets <- function(x) {
    if (length(x) == 0) return(list())
    out <- list()
    idx <- 1
    for (k in seq_along(x)) {
      cmb <- utils::combn(x, k, simplify = FALSE)
      for (v in cmb) {
        out[[idx]] <- v
        idx <- idx + 1
      }
    }
    out
  }
  
  # Ensure VarExp & CovarBio/CovarTec exist
  cd <- SummarizedExperiment::colData(se)
  used_vars <- unique(c(VarExp, CovarBio, CovarTec))
  missing_vars <- setdiff(used_vars, colnames(cd))
  if (length(missing_vars)) {
    stop("Variables not found in colData(se): ", paste(missing_vars, collapse = ", "))
  }
  
  used_vars <- used_vars[!is.na(used_vars) & nzchar(used_vars)]
  
  # Generate all non-empty subsets
  subsets <- all_subsets(used_vars)
  
  # Create formulas list (named like "A", "A+B", "A+B+C", etc.)
  formulas <- list()
  formulas[["null"]] <- create_formula(character(0))
  if (length(subsets)) {
    for (s in subsets) {
      nm <- paste(s, collapse = "+")
      formulas[[nm]] <- create_formula(s)
    }
  }
  # Return dictionary of formulas
  formulas
}

##########################################
model_prep <- function(
    se,
    VarExp = "Group",
    RefGroup = "Control",
    nSV = 0,
    VarNested = NULL,
    extra_formula = NULL,
    scalefor = NULL,
    CatVars = NULL
) {
  # ---- Setup ----
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  
  # ---- 1. Handle extra formula terms ----
  extra_terms <- character(0)
  if (!is.null(extra_formula)) {
    if (is.character(extra_formula)) {
      # Tolerate "x + y" or "~ x + y"
      if (!grepl("~", extra_formula, fixed = TRUE))
        extra_formula <- paste("~", extra_formula)
      extra_formula <- stats::as.formula(extra_formula)
    }
    extra_terms <- attr(stats::terms(extra_formula), "term.labels")
  }
  
  # ---- 2. Handle nested variable ----
  use_block <- FALSE
  if (!is.null(VarNested)) {
    if (identical(VarNested, VarExp)) {
      warning("VarNested equals VarExp (", VarNested, "); ignoring VarNested.")
      VarNested <- NULL
    } else {
      pat <- paste0("(^|:)", VarNested, "(:|$)")
      in_formula <- any(extra_terms %in% VarNested | grepl(pat, extra_terms))
      if (in_formula) {
        message("Using '", VarNested, "' as blocking factor and removing it (and its interactions) from fixed effects.")
        extra_terms <- extra_terms[!(extra_terms %in% VarNested | grepl(pat, extra_terms))]
        use_block <- TRUE
      }
    }
  }
  
  # ---- 3. Construct term vector ----
  sv_terms <- if (nSV > 0) paste0("SV", seq_len(nSV)) else character(0)
  term_vec <- unique(c(VarExp, extra_terms, sv_terms))
  if (!is.null(VarNested)) term_vec <- setdiff(term_vec, VarNested)
  
  # ---- 4. Scale selected numeric variables ----
  if (!is.null(scalefor)) {
    scalefor <- intersect(scalefor, colnames(cd))
    for (v in scalefor) {
      if (!is.numeric(cd[[v]]))
        stop("Trying to scale non-numeric variable: ", v)
      cd[[paste0(v, ".scaled")]] <- as.numeric(scale(cd[[v]]))
    }
    # Replace scaled names in formula terms
    term_vec[term_vec %in% scalefor] <- paste0(term_vec[term_vec %in% scalefor], ".scaled")
  }
  
  # ---- 5. Convert categorical variables to factors ----
  cats <- unique(c(CatVars, VarExp))
  cats <- cats[cats %in% colnames(cd)]
  for (nm in cats) {
    if (is.character(cd[[nm]]) || is.logical(cd[[nm]]))
      cd[[nm]] <- factor(cd[[nm]])
    if (is.factor(cd[[nm]]))
      cd[[nm]] <- droplevels(cd[[nm]])
  }
  
  # ---- 6. Set reference level for VarExp ----
  if (!is.null(RefGroup) && VarExp %in% colnames(cd)) {
    if (is.factor(cd[[VarExp]])) {
      if (!(RefGroup %in% levels(cd[[VarExp]]))) {
        warning("RefGroup '", RefGroup, "' not found in levels of ", VarExp, 
                "; available levels are: ", paste(levels(cd[[VarExp]]), collapse = ", "))
      } else {
        cd[[VarExp]] <- stats::relevel(cd[[VarExp]], ref = RefGroup)
        message("Reference level for ", VarExp, " set to '", RefGroup, "'.")
      }
    } else {
      warning("Variable ", VarExp, " is not a factor; cannot set reference level.")
    }
  }
  
  # ---- 7. Construct formula and model matrix ----
  formula <- stats::reformulate(term_vec)
  cat("\nThe applied formula is:", deparse(formula), "\n")
  
  model <- stats::model.matrix(formula, data = cd)
  
  # ---- 8. Ensure full rank ----
  qrM <- qr(model)
  if (qrM$rank < ncol(model)) {
    dropped <- colnames(model)[setdiff(seq_len(ncol(model)), qrM$pivot[seq_len(qrM$rank)])]
    message("Design not full rank. Dropping aliased columns: ", paste(dropped, collapse = ", "))
    model <- model[, qrM$pivot[seq_len(qrM$rank)], drop = FALSE]
  }
  
  # ---- 9. Write back scaled/factored data ----
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)
  
  # ---- 10. Return both ----
  return(list(
    se = se,
    model = model,
    formula = formula,
    use_block = use_block
  ))
}

##########################################
#TODO
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
                    trend = TRUE,
                    engine = 'edgeR_lrt',
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

  coef_name <- paste0(VarExp, ExpGroup)  # e.g., "GroupMSUS"
  
  # If coef was dropped by rank reduction, fail early with advice
  if (!coef_name %in% colnames(model)) {
    stop("The contrast '", coef_name, "' is not estimable with the current design. ",
         "If you intended to control for a nested factor, include it in extra_formula and set VarNested to its name.")
  }

  # Get matrix features*samples
  dge <- SummarizedExperiment::assay(se, assay_name)

  # ---------- Fit ----------
  if (use_block) {
    message("Use the limma voom technique because of nested variable")
    nested <- SummarizedExperiment::colData(se)[[VarNested]]
    v <- limma::voom(dge, model)
    corfit <- limma::duplicateCorrelation(v, model, block = nested)
    v <- limma::voom(dge, model, block = nested, correlation = corfit$consensus.correlation)
    fit_ <- limma::lmFit(v, model, block = nested, correlation = corfit$consensus.correlation)
    fit <- limma::eBayes(fit_, trend = trend, robust = disp_robust)
    res <- limma::topTable(fit, coef = coef_name, number = nrow(se), adjust.method = "BH")
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

      fit_ <- limma::lmFit(v, design = model)
  
      fit <- limma::eBayes(fit_, trend = trend, robust = disp_robust)
      res <- limma::topTable(
        fit,
        coef          = coef_name,       # cd[[VarExp]]ExpGroup vs RefGroup
        number        = nrow(se),
        adjust.method = "BH"
      )
      res <- res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
      colnames(res) <- c("logFC", "avgExpr", "stat", "pval", "padj", "B")
    }
    if (engine == "edgeR_ql") {
      cat("Use the ql estimate")
      fit <- edgeR::glmQLFit(dge, model, robust = disp_robust)
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
  # Get feature name  
  if (!is.null(row_name)) {
    rn <- as.character(SummarizedExperiment::rowData(se)[[row_name]])
    print(table(duplicated(rowData(se)[[row_name]])))

  } else {
    rn <- rownames(se)
  }

  rownames(res) <- make.unique(rn, sep = "_")
  SummarizedExperiment::rowData(se)$dma <- res

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

  se_out <- se

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
  
  list(se = se_out, res = res, df = df)
}

##########################################
#pval dist
plot_pval_distribution <- function(
    se,
    save_path   = NULL,
    line_at     = 0.05,
    alpha       = 0.05,
    title       = "P value distribution",
    # what to report in the corner: BH, FDRtool, or both
    binwidth    = 0.05
) {
  
  # Grab the fdrtool table from rowData
  dma_tbl <- as.data.frame(SummarizedExperiment::rowData(se)$dma)
  
  # Safety checks
  needed_cols <- c("pval", "padj")
  missing_cols <- setdiff(needed_cols, names(dma_tbl))
  if (length(missing_cols) > 0) {
    stop("Missing columns in rowData(se)$fdrtool: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Choose which p-values to plot on x-axis
  pvals <- dma_tbl$pval
  pvals <- pvals[!is.na(pvals)]
  
  # Counts for BH and FDRtool
  n_sig_bh      <- sum(dma_tbl$padj        < alpha, na.rm = TRUE)
  
  # Build label based on sig_type
  label <- paste0("BH FDR < ", alpha, ": ", n_sig_bh)

  
  
  # Compute histogram counts to place annotation nicely
  breaks  <- seq(0, 1, by = binwidth)
  h       <- hist(pvals, breaks = breaks, plot = FALSE)
  ymax    <- if (length(h$counts) > 0) max(h$counts, na.rm = TRUE) else 0
  y_annot <- ymax * 1.15
  
  pvals_long <- data.frame(pval = pvals)
  
  pval_plot <- ggplot(pvals_long, aes(x = pval)) +
    geom_histogram(
      color  = "#e9ecef",
      alpha  = 0.9,
      breaks = breaks
    ) +
    geom_vline(xintercept = line_at, linetype = 2) +
    annotate(
      "text",
      x      = 0.98,
      y      = y_annot,
      label  = label,
      hjust  = 1,
      vjust  = 0,
      size   = 3
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    coord_cartesian(ylim = c(0, ymax * 1.35)) +
    labs(
      x     = "Raw p-values",
      y     = "Number of features",
      title = title
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 9, face = "bold", margin = margin(b = 3)),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 7)
    )
  
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    pdf(save_path)
    print(pval_plot)
    dev.off()
  } else {
    print(pval_plot)
  }
  pval_plot
}

##########################################
plot_volcano <- function(
  se,
  VarExp = "Group",
  ExpGroup = "MSUS",
  RefGroup = "Control",
  title_fc  = "Volcano: log2 fold change",
  title_pmd = "Volcano: methylation difference",
  convert = TRUE,
  convert_ids = TRUE,
  organism = "mmusculus",
  FCcutoff = 1.5,
  pmd_cutoff = 10,            # in percentage points
  pCutoff = 0.05,
  save_path = NULL
) {
  grp <- SummarizedExperiment::colData(se)[[VarExp]]
  nCtr <- sum(grp == RefGroup, na.rm = TRUE)
  nExp <- sum(grp == ExpGroup, na.rm = TRUE)

  tmp <- SummarizedExperiment::rowData(se)$dma

  gene_symbol <- if ("gene_symbol" %in% names(tmp)) as.character(tmp[["gene_symbol"]]) else rep(NA_character_, nrow(tmp))
  gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- as.character(tmp[["feature_id"]])

  res_vol <- data.frame(
    feature_id    = as.character(tmp[["feature_id"]]),
    gene_symbol   = gene_symbol,
    logFC         = tmp[["logFC"]],
    per_meth_diff = tmp[["per_meth_diff"]],
    pval          = tmp[["pval"]],
    padj          = tmp[["padj"]],
    row.names     = as.character(tmp[["feature_id"]])
  )
  # labels
  res_vol$Gene <- if (convert) {
    convertGeneNames(res_vol$feature_id, convert_ids = convert_ids, organism = organism)
  } else {
    res_vol$gene_symbol
  }

  # caption + optional FDR-derived hline
  nSignP   <- sum(res_vol$pval < pCutoff, na.rm = TRUE)
  nSignFDR <- sum(res_vol$padj < 0.05,   na.rm = TRUE)

  if (nSignFDR > 0) {
    hlinetype  <- "solid"
    thresh_fdr <- max(res_vol$pval[res_vol$padj < 0.05], na.rm = TRUE)
  } else {
    hlinetype  <- NULL
    thresh_fdr <- NULL
  }

  caption <- paste0(
    "n features = ", nrow(res_vol),
    "\nN p-value < ", pCutoff, " = ", nSignP,
    "\nN FDR-adjusted p-value < 0.05 = ", nSignFDR,
    "\nSample size: ", RefGroup, " = ", nCtr, " and ", ExpGroup, " = ", nExp
  )

  # label sets (separate, so each plot labels what matters for that x-variable)
  sel_fc  <- which(abs(res_vol$logFC)         > FCcutoff | res_vol$padj < 0.05)
  sel_pmd <- which(abs(res_vol$per_meth_diff) > pmd_cutoff | res_vol$padj < 0.05)

  # symmetric xlim for methylation difference, adaptive but at least ±50

  rng_pmd <- range(res_vol$per_meth_diff, na.rm = TRUE)
  lim_pmd <- c(rng_pmd[1] - 2, rng_pmd[2] + 2)

  VolPlot_fc <- EnhancedVolcano::EnhancedVolcano(
    res_vol,
    lab            = res_vol$Gene,
    selectLab      = res_vol$Gene[sel_fc],
    x              = "logFC",
    y              = "pval",
    xlab           = "log2 fold change",
    pCutoff        = pCutoff,
    FCcutoff       = FCcutoff,
    hline          = thresh_fdr,
    hlineType      = hlinetype,
    legendPosition = "bottom",
    pointSize      = 2,
    labSize        = 3.5,
    legendIconSize = 4.0,
    xlim           = c(min(res_vol$logFC, na.rm = TRUE) - 0.5,
                       max(res_vol$logFC, na.rm = TRUE) + 0.5),
    ylim           = c(0, max(-log10(res_vol$pval), na.rm = TRUE) + 1),
    title          = title_fc,
    subtitle       = "",
    caption        = caption,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    drawConnectors  = TRUE,
    widthConnectors = 0.5,
    boxedLabels     = FALSE
  )

  VolPlot_pmd <- EnhancedVolcano::EnhancedVolcano(
    res_vol,
    lab            = res_vol$Gene,
    selectLab      = res_vol$Gene[sel_pmd],
    x              = "per_meth_diff",
    y              = "pval",
    xlab           = "Methylation difference (percentage points)",
    pCutoff        = pCutoff,
    FCcutoff       = pmd_cutoff,
    hline          = thresh_fdr,
    hlineType      = hlinetype,
    legendPosition = "bottom",
    pointSize      = 2,
    labSize        = 3.5,
    legendIconSize = 4.0,
    xlim           = c(-lim_pmd, lim_pmd),
    ylim           = c(0, max(-log10(res_vol$pval), na.rm = TRUE) + 1),
    title          = title_pmd,
    subtitle       = "",
    caption        = caption,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    drawConnectors  = TRUE,
    widthConnectors = 0.5,
    boxedLabels     = FALSE
  )

  if (!is.null(save_path)) {
    dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(save_path, onefile = TRUE, width = 10, height = 8)
    print(VolPlot_fc)
    print(VolPlot_pmd)
    grDevices::dev.off()
  }

  print(VolPlot_fc)
  print(VolPlot_pmd)

  invisible(list(volcano_logFC = VolPlot_fc, volcano_pmd = VolPlot_pmd, data = res_vol))
}

