suppressPackageStartupMessages({
  library(edgeR)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(HDF5Array)
})

#' Convert Me/Un edgeR DGEList to an HDF5-backed SummarizedExperiment.
#'
#' @description
#' Reads an edgeR DGEList saved as an RDS. The DGEList is assumed to contain:
#' (i) `counts` with paired columns per biological sample suffixed by `-Me` and `-Un`,
#' (ii) `samples` with one row per `counts` column, and
#' (iii) optional `genes` with one row per feature (CpG).
#'
#' The function computes percent methylation (MePer), and optionally Beta values
#' and M-values. It then saves an HDF5SummarizedExperiment to disk as a directory
#' containing `assays.h5` and `se.rds`. No standalone assay `.h5` files are created.
#'
#' @param me_dge_path Path to the input `.rds` containing an edgeR `DGEList`.
#' @param output_dir Directory where the HDF5SummarizedExperiment directory will be written.
#' @param output_name Output basename (with or without `.rds` extension). The final
#'   directory will be `file.path(output_dir, <stub>)` where `<stub>` is `output_name`
#'   without extension.
#' @param assay_name Name for the percent methylation assay (default `"MePer"`).
#' @param compute_beta_M Logical; if TRUE, compute Beta and M assays (default TRUE).
#' @param beta_assay_name Name for the Beta assay (default `"Bval"`).
#' @param m_assay_name Name for the M-value assay (default `"Mval"`).
#' @param chunk_nrows Number of rows per HDF5 chunk when saving assays (default 50000).
#'
#' @details
#' Pairing logic:
#' - Methylated columns are those ending in `-Me`.
#' - The corresponding unmethylated column is obtained by replacing `-Me` with `-Un`.
#' - Only complete pairs are used.
#'
#' Assays:
#' - MePer: `100 * Me / (Me + Un)` with non-finite values set to NA.
#' - Beta: `MePer / 100` with non-finite values set to NA.
#' - Mval: `log2(beta / (1 - beta))` after clamping beta to `[1e-6, 1-1e-6]`.
#'
#' colData:
#' A single row per biological sample is derived from `y$samples` by stripping `-(Me|Un)$`
#' from rownames to form sample IDs, keeping the first occurrence per sample, and
#' re-ordering to match assay column order.
#'
#' @return Invisibly returns the in-memory `SummarizedExperiment` object.
#'   The on-disk representation is written to `file.path(output_dir, <stub>)`.
#'
#' @examples
#' \dontrun{
#' me_dge2mepe_se(
#'   me_dge_path = "F1_yall_me_dge_sgCpGs_clean.rds",
#'   output_dir  = "02_sgCpGs_clean",
#'   output_name = "F1_mepe_se_sgCpGs_clean.rds"
#' )
#' se <- HDF5Array::loadHDF5SummarizedExperiment("02_sgCpGs_clean/F1_mepe_se_sgCpGs_clean")
#' }
me_dge2mepe_se <- function(me_dge_path,
                           output_dir,
                           output_name,
                           assay_name       = "MePer",
                           compute_beta_M   = TRUE,
                           beta_assay_name  = "Bval",
                           m_assay_name     = "Mval",
                           chunk_nrows      = 50000) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load edgeR DGEList (expects $counts, $samples, and optional $genes).
  y <- readRDS(me_dge_path)

  # ---- 1) Identify complete Me/Un column pairs and compute percent methylation ----
  Me_cols <- grep("-Me$", colnames(y$counts), value = TRUE)
  Un_cols <- sub("-Me$", "-Un", Me_cols)

  ok <- Un_cols %in% colnames(y$counts)
  Me_cols <- Me_cols[ok]
  Un_cols <- Un_cols[ok]
  if (!length(Me_cols)) stop("No complete -Me/-Un pairs found in y$counts.")

  # Matrices are (features x technical columns). After pairing, columns correspond to biological samples.
  Me <- y$counts[, Me_cols, drop = FALSE]
  Un <- y$counts[, Un_cols, drop = FALSE]

  meper <- 100 * Me / (Me + Un)
  meper[!is.finite(meper)] <- NA_real_
  colnames(meper) <- sub("-Me$", "", Me_cols)

  # ---- 2) Build one-row-per-biological-sample phenotype table from y$samples ----
  samp <- as.data.frame(y$samples)
  samp$.raw <- rownames(y$samples)
  samp$Sample <- sub("-(Me|Un)$", "", samp$.raw)

  # Keep first row per biological sample, then align to assay column order.
  pheno <- samp[!duplicated(samp$Sample), , drop = FALSE]
  rownames(pheno) <- pheno$Sample
  pheno$.raw <- NULL
  pheno$Sample <- NULL

  pheno <- pheno[colnames(meper), , drop = FALSE]
  pheno$Sample <- rownames(pheno)

  # ---- 3) Construct assays (in-memory) ----
  assays_list <- S4Vectors::SimpleList()
  assays_list[[assay_name]] <- meper

  if (compute_beta_M) {
    beta <- meper / 100
    beta[!is.finite(beta)] <- NA_real_

    eps <- 1e-6
    beta_clamped <- pmin(pmax(beta, eps), 1 - eps)
    Mval <- log2(beta_clamped / (1 - beta_clamped))

    assays_list[[beta_assay_name]] <- beta
    assays_list[[m_assay_name]]    <- Mval
  }

  # ---- 4) Build SummarizedExperiment ----
  se <- SummarizedExperiment(
    assays  = assays_list,
    rowData = if (!is.null(y$genes)) S4Vectors::DataFrame(y$genes) else NULL,
    colData = S4Vectors::DataFrame(pheno)
  )

  # ---- 5) Save as HDF5SummarizedExperiment (single assays.h5 inside se_dir) ----
  stub <- tools::file_path_sans_ext(output_name)
  se_dir <- file.path(output_dir, stub)
  if (dir.exists(se_dir)) unlink(se_dir, recursive = TRUE, force = TRUE)

  # Store assays as chunks of (chunk_nrows x 1 column) to keep I/O bounded.
  HDF5Array::saveHDF5SummarizedExperiment(
    se,
    dir      = se_dir,
    replace  = TRUE,
    chunkdim = c(chunk_nrows, 1)
  )

  message("Saved HDF5SummarizedExperiment to: ", se_dir)
  message("Load with: se <- HDF5Array::loadHDF5SummarizedExperiment('", se_dir, "')")

  invisible(se)
}
