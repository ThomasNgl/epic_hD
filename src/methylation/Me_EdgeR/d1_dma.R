suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(edgeR)
  library(limma)
  library(sva)
})

run_sva_glm <- function(
  se_path,
  dge_rds,
  output_dir,
  m_assay = "M_from_dge",
  group_var = "Group",
  contrast_levels = c("MSUS", "Control"),
  covariates = NULL,
  nsv_override = NULL
) {
  if (!exists("modelMatrixMeth", mode = "function")) {
    stop("modelMatrixMeth() not found on search path. Source/load it before calling run_sva_glm().")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("start")

  ## ---- load inputs ----
  se <- if (inherits(se_path, "SummarizedExperiment")) {
    se_path
  } else if (dir.exists(se_path)) {
    HDF5Array::loadHDF5SummarizedExperiment(se_path)
  } else if (file.exists(se_path) && grepl("\\.rds$", se_path, ignore.case = TRUE)) {
    readRDS(se_path)
  } else {
    stop("se_path must be a SummarizedExperiment, an HDF5SummarizedExperiment directory, or an .rds file: ", se_path)
  }

  dge <- readRDS(dge_rds)


  targets <- as.data.frame(colData(se))
  if (!group_var %in% names(targets)) stop("group_var not found in colData(se): ", group_var)

  M <- assay(se, m_assay)
  if (!is.matrix(M)) M <- as.matrix(M)

  # -------------------------
  # Align on BASE sample names (M columns), then expand to -Me/-Un for DGE
  # -------------------------
  base_M <- colnames(M)
  if (is.null(base_M)) stop("M assay has no colnames; expected base sample IDs.")

  dge_cols <- colnames(dge)
  if (is.null(dge_cols)) stop("DGEList has no colnames.")

  dge_base <- sub("-(Me|Un)$", "", dge_cols)
  has_meun <- grepl("-(Me|Un)$", dge_cols)
  if (!all(has_meun)) stop("Not all DGE colnames end with -Me/-Un; cannot align automatically.")

  common_base <- intersect(base_M, unique(dge_base))
  if (!length(common_base)) {
    stop("No overlap between M colnames (base samples) and DGE base names (from -Me/-Un stripping).")
  }

  # keep M in its current order; drop anything not in DGE
  keep_base <- base_M[base_M %in% common_base]
  M <- M[, keep_base, drop = FALSE]

  # align targets to base sample order: prefer rownames(targets), else a 'Sample' column, else fail
  if (!is.null(rownames(targets)) && all(keep_base %in% rownames(targets))) {
    targets <- targets[keep_base, , drop = FALSE]
  } else if ("Sample" %in% names(targets) && all(keep_base %in% as.character(targets$Sample))) {
    targets <- targets[match(keep_base, as.character(targets$Sample)), , drop = FALSE]
    rownames(targets) <- keep_base
  } else {
    stop("Cannot align colData to M colnames. Need rownames(colData) == base samples or a 'Sample' column containing them.")
  }

  # subset/reorder DGE columns to match keep_base expanded to Me/Un
  dge_cols_needed <- as.vector(rbind(paste0(keep_base, "-Me"), paste0(keep_base, "-Un")))
  if (!all(dge_cols_needed %in% dge_cols)) {
    miss <- setdiff(dge_cols_needed, dge_cols)
    stop("Missing expected DGE columns (first few): ", paste(head(miss, 10), collapse = ", "))
  }
  dge <- dge[, dge_cols_needed, keep.lib.sizes = FALSE]

  # -------------------------
  # SVA on M (base samples)
  # -------------------------
  G <- factor(targets[[group_var]])
  if (!all(contrast_levels %in% levels(G))) {
    stop("contrast_levels not present in ", group_var, ": ",
         paste(setdiff(contrast_levels, levels(G)), collapse = ", "))
  }
  G <- relevel(G, ref = contrast_levels[2])
  targets$Group <- G  # standardize name for formulas

  cov_df <- if (is.null(covariates)) {
    targets[, "Group", drop = FALSE]
  } else {
    missing_cov <- setdiff(covariates, names(targets))
    if (length(missing_cov)) stop("Missing covariates in colData(se): ", paste(missing_cov, collapse = ", "))
    targets[, c("Group", covariates), drop = FALSE]
  }

  mod  <- model.matrix(~ Group + ., data = cov_df)
  cov0 <- cov_df[, setdiff(colnames(cov_df), "Group"), drop = FALSE]
  mod0 <- if (ncol(cov0) == 0) model.matrix(~ 1, data = targets) else model.matrix(~ 1 + ., data = cov0)

  message("Running SVA on assay: ", m_assay)
  svobj <- sva::sva(M, mod, mod0)

  nSV <- svobj$n.sv
  if (!is.null(nsv_override)) nSV <- min(nSV, as.integer(nsv_override))
  message("Using ", nSV, " surrogate variables.")

  sv_df <- as.data.frame(svobj$sv)
  if (nSV == 0) {
    sv_df <- sv_df[, 0, drop = FALSE]
  } else {
    sv_df <- sv_df[, seq_len(nSV), drop = FALSE]
    colnames(sv_df) <- paste0("SV", seq_len(nSV))
  }
  rownames(sv_df) <- rownames(targets)

  targets_with_sv <- as.data.frame(cbind(cov_df, sv_df))  # base samples, rownames preserved

  # -------------------------
  # edgeR GLM QLF (Me/Un counts) using modelMatrixMeth expansion
  # -------------------------
  design_terms <- c("Group", covariates, colnames(sv_df))
  design_terms <- design_terms[!is.na(design_terms) & nzchar(design_terms)]

  design_svs <- model.matrix(
    as.formula(paste("~0 +", paste(design_terms, collapse = " + "))),
    data = targets_with_sv
  )
  design_svs_matmeth <- modelMatrixMeth(design_svs)

  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design_svs_matmeth, robust = TRUE)
  fit_svs <- glmQLFit(dge, design_svs_matmeth, robust = TRUE)

  c1 <- paste0("Group", make.names(contrast_levels[1]))
  c0 <- paste0("Group", make.names(contrast_levels[2]))
  if (!all(c(c1, c0) %in% colnames(design_svs))) {
    stop("Could not find group columns in design. Expected: ", c1, " and ", c0,
         ". Found: ", paste(colnames(design_svs), collapse = ", "))
  }

  contr <- makeContrasts(contrasts = paste0(c1, " - ", c0), levels = design_svs_matmeth)
  glm_svs <- glmQLFTest(fit_svs, contrast = contr)

  glm_edgeR_svs_path <- file.path(output_dir, "edger_glm_svs.rds")
  saveRDS(glm_svs, glm_edgeR_svs_path)

  # -------------------------
  # limma on M with SVs (base samples)
  # -------------------------
  design_limma <- model.matrix(
    as.formula(paste("~", paste(design_terms, collapse = " + "))),
    data = targets_with_sv
  )
  fit_limma <- eBayes(lmFit(M, design_limma))

  glm_limma_svs_path <- file.path(output_dir, "limma_glm_svs.rds")
  saveRDS(fit_limma, glm_limma_svs_path)

  # -------------------------
  # Add SVs to SE colData and save
  # -------------------------
  cd <- as.data.frame(colData(se))
  if (!is.null(rownames(cd)) && all(rownames(sv_df) %in% rownames(cd))) {
    cd <- cd[rownames(sv_df), , drop = FALSE]
  } else if ("Sample" %in% names(cd) && all(rownames(sv_df) %in% as.character(cd$Sample))) {
    cd <- cd[match(rownames(sv_df), as.character(cd$Sample)), , drop = FALSE]
    rownames(cd) <- rownames(sv_df)
  } else {
    stop("Cannot align SE colData for SV insertion (need rownames or Sample column matching base samples).")
  }

  colData(se) <- S4Vectors::DataFrame(as.data.frame(cbind(cd, sv_df)))

  # save as HDF5SummarizedExperiment directory (NOT .rds)
  se_svs_dir <- file.path(output_dir, "se_svs")
  if (dir.exists(se_svs_dir)) unlink(se_svs_dir, recursive = TRUE, force = TRUE)
  HDF5Array::saveHDF5SummarizedExperiment(se, dir = se_svs_dir, replace = TRUE)

  saveRDS(
    list(svobj = svobj, targets_with_sv = targets_with_sv),
    file.path(output_dir, "sva_details.rds")
  )

  invisible(list(
    nSV = nSV,
    paths = list(
      edger_glm_svs = glm_edgeR_svs_path,
      limma_glm_svs = glm_limma_svs_path,
      se_svs = se_svs_dir
    )
  ))

}
