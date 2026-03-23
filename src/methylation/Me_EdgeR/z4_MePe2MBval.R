suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(HDF5Array)
})

#' Add Bval and Mval assays to an existing MePer SummarizedExperiment (HDF5-friendly)
#'
#' - Input: SE object, HDF5SummarizedExperiment directory, or .rds SE
#' - Uses existing assay (default "MePer", 0–100)
#' - Computes:
#'     * Bval = MePer / 100 (0–1)
#'     * Mval = log2(Bval / (1 - Bval)), with clamping at [eps, 1-eps]
#' - Output: new HDF5SummarizedExperiment in `file.path(output_dir, file_path_sans_ext(output_name))`
#'           with three assays: MePer, Bval, Mval
#'
#' @param se_input SummarizedExperiment object, HDF5 SE directory, or path to SE .rds.
#' @param output_dir Directory where the new SE will be written.
#' @param output_name Base name for the output (can be with or without .rds; only the stem is used).
#' @param assay_name Name of the existing Me-percentage assay (default "MePer").
#' @param beta_assay_name Name of the beta assay to create (default "Bval").
#' @param m_assay_name Name of the M-value assay to create (default "Mval").
#' @param chunk_nrows Chunk size for HDF5 writing (per row block).
#'
#' @return SummarizedExperiment (invisibly); written to disk as HDF5SummarizedExperiment.
add_beta_M_to_mepe_se <- function(se_input,
                                  output_dir,
                                  output_name,
                                  assay_name      = "MePer",
                                  beta_assay_name = "Bval",
                                  m_assay_name    = "Mval",
                                  chunk_nrows     = 50000) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ## ---- 1) Load SE ----
  se <- if (inherits(se_input, "SummarizedExperiment")) {
    se_input
  } else if (is.character(se_input) && length(se_input) == 1L) {
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

  if (!(assay_name %in% assayNames(se))) {
    stop("Assay '", assay_name, "' not found in SummarizedExperiment.")
  }

  ## ---- 2) Extract MePer and compute Bval and Mval ----
  meper <- SummarizedExperiment::assay(se, assay_name)
  if (!is.matrix(meper)) meper <- as.matrix(meper)

  # MePer expected 0–100; allow NAs and infinities (convert to NA)
  meper[!is.finite(meper)] <- NA_real_

  # Beta 0–1
  beta <- meper / 100
  beta[!is.finite(beta)] <- NA_real_

  # M-values from beta
  eps <- 1e-6
  beta_clamped <- pmin(pmax(beta, eps), 1 - eps)
  Mval <- log2(beta_clamped / (1 - beta_clamped))

  ## ---- 3) Write all three assays as HDF5-backed arrays ----
  assays_list <- S4Vectors::SimpleList()

  base_stub <- file.path(output_dir, tools::file_path_sans_ext(output_name))

  # HDF5 file paths
  h5_meper <- paste0(base_stub, "_", assay_name,      ".h5")
  h5_beta  <- paste0(base_stub, "_", beta_assay_name, ".h5")
  h5_M     <- paste0(base_stub, "_", m_assay_name,    ".h5")

  # Remove old files if present
  for (f in c(h5_meper, h5_beta, h5_M)) {
    if (file.exists(f)) unlink(f)
  }

  chunkdim <- c(min(nrow(meper), chunk_nrows), min(ncol(meper), 1))

  meper_h5 <- HDF5Array::writeHDF5Array(
    meper,
    filepath      = h5_meper,
    name          = assay_name,
    with.dimnames = TRUE,
    chunkdim      = chunkdim,
    verbose       = FALSE
  )
  beta_h5 <- HDF5Array::writeHDF5Array(
    beta,
    filepath      = h5_beta,
    name          = beta_assay_name,
    with.dimnames = TRUE,
    chunkdim      = chunkdim,
    verbose       = FALSE
  )
  M_h5 <- HDF5Array::writeHDF5Array(
    Mval,
    filepath      = h5_M,
    name          = m_assay_name,
    with.dimnames = TRUE,
    chunkdim      = chunkdim,
    verbose       = FALSE
  )

  assays_list[[assay_name]]      <- meper_h5
  assays_list[[beta_assay_name]] <- beta_h5
  assays_list[[m_assay_name]]    <- M_h5

  rm(meper, beta, beta_clamped, Mval); gc()

  ## ---- 4) Build new SE with same rowData/colData, new assays ----
  new_se <- SummarizedExperiment::SummarizedExperiment(
    assays  = assays_list,
    rowData = SummarizedExperiment::rowData(se),
    colData = SummarizedExperiment::colData(se)
  )

  ## ---- 5) Save as HDF5SummarizedExperiment ----
  se_dir <- file.path(output_dir, tools::file_path_sans_ext(output_name))
  if (dir.exists(se_dir)) {
    unlink(se_dir, recursive = TRUE, force = TRUE)
  }

  HDF5Array::saveHDF5SummarizedExperiment(new_se, dir = se_dir, replace = TRUE)
  message("Saved HDF5SummarizedExperiment with assays [",
          paste(names(assays_list), collapse = ", "),
          "] to: ", se_dir)
  message("Load with: se <- HDF5Array::loadHDF5SummarizedExperiment('", se_dir, "')")

  invisible(new_se)
}
