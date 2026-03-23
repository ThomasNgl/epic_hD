suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(HDF5Array)
})

#' Remove samples from a MePer SummarizedExperiment and optionally binarize (HDF5-friendly)
#'
#' - Input: SE object, HDF5SummarizedExperiment directory, or .rds SE
#' - Optional: drop samples
#' - Optional: create a binarized assay (0 / 100 / NA) with user-defined thresholds
#' - Output: HDF5SummarizedExperiment directory: file.path(output_dir, output_name)
#'
#' @param se_or_path SummarizedExperiment object, HDF5 SE directory, or path to SE RDS.
#' @param samples_to_remove Character vector of sample IDs to drop (colnames of SE).
#' @param output_dir Directory to save output SE (as HDF5SummarizedExperiment).
#' @param output_name Directory name for saved SE (not a .rds filename).
#' @param assay_name Name of input assay (default "MePer").
#' @param binarize FALSE or numeric vector c(low, high).
#'        FALSE: no binarization.
#'        c(low, high): NA stays NA; <low -> 0; >high -> 100; [low, high] -> NA.
#' @param bin_assay_name Name for the binarized assay if created (default "MePerBin").
#' @return SummarizedExperiment (invisibly); written to disk as HDF5SummarizedExperiment.
clean_mepe_se <- function(se_or_path,
                          samples_to_remove = NULL,
                          output_dir,
                          output_name,
                          assay_name   = "MePer",
                          binarize     = FALSE,          # FALSE or c(low, high)
                          bin_assay_name = "MePerBin") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ## ---- 1) Load SE ----
  se <- if (inherits(se_or_path, "SummarizedExperiment")) {
    se_or_path
  } else if (dir.exists(se_or_path)) {
    HDF5Array::loadHDF5SummarizedExperiment(se_or_path)
  } else {
    readRDS(se_or_path)
  }

  ## ---- 2) Remove samples (if requested) ----
  if (!is.null(samples_to_remove) && length(samples_to_remove)) {
    se <- se[, !(colnames(se) %in% samples_to_remove), drop = FALSE]
  }

  ## ---- 3) Bring the main assay fully into memory, clean infinities ----
  A <- as.matrix(SummarizedExperiment::assay(se, assay_name))
  A[!is.finite(A)] <- NA_real_
  SummarizedExperiment::assay(se, assay_name) <- A

  ## ---- 4) Optional binarization ----
  if (!identical(binarize, FALSE)) {
    if (!is.numeric(binarize) || length(binarize) != 2) {
      stop("'binarize' must be FALSE or a numeric vector c(low, high).")
    }
    low  <- binarize[1]
    high <- binarize[2]
    if (!is.finite(low) || !is.finite(high) || low >= high) {
      stop("Invalid binarize thresholds: require finite low < high.")
    }

    B <- A
    B[B < low]             <- 0
    B[B > high]            <- 100
    B[B >= low & B <= high] <- NA_real_

    SummarizedExperiment::assay(se, bin_assay_name) <- B
  }

  ## ---- 5) Save as HDF5SummarizedExperiment (NOT saveRDS) ----
  out_dir <- file.path(output_dir, output_name)
  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }

  HDF5Array::saveHDF5SummarizedExperiment(se, dir = out_dir, replace = TRUE)

  message("Saved HDF5SummarizedExperiment to: ", out_dir)
  invisible(se)
}
