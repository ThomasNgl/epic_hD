svacor_from_colData <- function(SE, form, form0 = ~1,
                                assayName = 1,
                                sv_pattern = "^SV[0-9]+$",
                                regressOutNull = TRUE,
                                corrected_assay = "corrected") {
  stopifnot(is(SE, "SummarizedExperiment"))

  CD <- as.data.frame(colData(SE))

  # SVs must already be in colData
  sv_cols <- grep(sv_pattern, names(CD), value = TRUE)
  if (length(sv_cols) == 0) stop("No SV columns found in colData (pattern: ", sv_pattern, ").")

  en <- as.matrix(assays(SE)[[assayName]])

  # Ensure sample order matches
  if (!identical(colnames(en), rownames(CD))) {
    CD <- CD[colnames(en), , drop = FALSE]
  }

  mm  <- model.matrix(form,  data = CD)
  mm0 <- model.matrix(form0, data = CD)
  sv  <- as.matrix(CD[, sv_cols, drop = FALSE])

  X <- cbind(mm, sv)

  # Fit coefficients (more stable than solve(t(X)%*%X))
  B <- qr.solve(X, t(en))  # dim: ncol(X) x nrow(en)

  # Which columns to regress out?
  if (regressOutNull) {
    keep_terms <- setdiff(colnames(mm), colnames(mm0))     # variables of interest
    cn <- setdiff(colnames(X), keep_terms)                 # remove everything else
  } else {
    cn <- setdiff(colnames(X), colnames(mm))               # remove only SVs
  }
  cn <- setdiff(cn, "(Intercept)")

  encor <- en
  if (length(cn) > 0) {
    encor <- en - t(X[, cn, drop = FALSE] %*% B[cn, , drop = FALSE])
  }

  assays(SE)[[corrected_assay]] <- encor
  SE
}
