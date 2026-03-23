suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(HDF5Array)
})

#' Merge plus/minus strand CpG rows into one CpG row with a single MePer assay
#'
#' - Input: SummarizedExperiment, HDF5SummarizedExperiment directory, or .rds SE
#' - Logic (CpG coordinate, always in plus-strand notation):
#'     * plus strand:  CpG = Locus            (Strand == "+")
#'     * minus strand: CpG = Locus - 1        (Strand == "-")
#' - All CpGs are kept:
#'     * "paired"     : both plus and minus rows exist
#'     * "plus_only"  : only plus row exists
#'     * "minus_only" : only minus row exists
#' - For each CpG and sample:
#'     * both NA                  -> NA#'     * only one strand non-NA   -> that value
#'     * both strands non-NA      -> mean(plus, minus)
#'
#' - Output: HDF5SummarizedExperiment directory with a single assay "MePer",
#'           rowData columns: Chr, Locus (CpG coordinate), Strand="merged",
#'           Source ∈ { "paired", "plus_only", "minus_only" }.
#'
#' @param se_input SummarizedExperiment object, HDF5 SE directory, or path to SE .rds.
#' @param output_dir Directory where the merged SE will be written.
#' @param output_name Name of the output HDF5SummarizedExperiment directory (not a .rds file).
#' @param assay_name Name of the assay to merge (default "MePer").
#'
#' @return Merged SummarizedExperiment (invisibly); written to disk as HDF5SummarizedExperiment.
merge_cpg_strands_meper <- function(se_input,
                                    output_dir,
                                    output_name,
                                    assay_name = "MePer") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ## 1) Load SE -------------------------------------------------------------
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

  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  if (!all(c("Chr", "Locus", "Strand") %in% colnames(rd))) {
    stop("rowData(se) must contain columns: Chr, Locus, Strand.")
  }

  A <- SummarizedExperiment::assay(se, assay_name)
  if (!is.matrix(A)) A <- as.matrix(A)

  ## 2) Build plus/minus index tables and CpG coordinates -------------------
  plus_idx  <- which(rd$Strand == "+")
  minus_idx <- which(rd$Strand == "-")

  if (!length(plus_idx) && !length(minus_idx)) {
    stop("merge_cpg_strands_meper: no '+' or '-' strand rows found in rowData(se)$Strand.")
  }

  plus_df <- data.frame(
    Chr      = as.character(rd$Chr[plus_idx]),
    CpG      = rd$Locus[plus_idx],     # plus-strand CpG coordinate
    idx_plus = plus_idx,
    stringsAsFactors = FALSE
  )

  minus_df <- data.frame(
    Chr       = as.character(rd$Chr[minus_idx]),
    CpG       = rd$Locus[minus_idx] - 1L,  # convert to plus-strand CpG coordinate
    idx_minus = minus_idx,
    stringsAsFactors = FALSE
  )

  ## (optional but safer) enforce uniqueness per (Chr, CpG, strand) ----------
  if (nrow(plus_df)) {
    dup_plus <- duplicated(plus_df[c("Chr", "CpG")])
    if (any(dup_plus)) {
      stop("Multiple '+' rows per (Chr,CpG); cannot merge safely.")
    }
  }
  if (nrow(minus_df)) {
    dup_minus <- duplicated(minus_df[c("Chr", "CpG")])
    if (any(dup_minus)) {
      stop("Multiple '-' rows per (Chr,CpG); cannot merge safely.")
    }
  }

  ## 3) Union of CpGs: keep all (paired + plus_only + minus_only) -----------
  # Base R merge with all=TRUE gives full outer join.
  pairs <- merge(plus_df, minus_df, by = c("Chr", "CpG"), all = TRUE)

  # After merge, we expect columns: Chr, CpG, idx_plus, idx_minus.
  if (!nrow(pairs)) {
    stop("merge_cpg_strands_meper: no CpGs found after merging plus/minus tables.")
  }

  # Source classification per CpG
  has_plus  <- !is.na(pairs$idx_plus)
  has_minus <- !is.na(pairs$idx_minus)
  Source <- ifelse(
    has_plus & has_minus, "paired",
    ifelse(has_plus, "plus_only", "minus_only")
  )

  ## 4) Build Ap (plus) and Am (minus) matrices aligned to 'pairs' ----------
  n_cpg    <- nrow(pairs)
  n_sample <- ncol(A)

  Ap <- matrix(NA_real_, nrow = n_cpg, ncol = n_sample,
               dimnames = list(NULL, colnames(A)))
  Am <- matrix(NA_real_, nrow = n_cpg, ncol = n_sample,
               dimnames = list(NULL, colnames(A)))

  idxp_non_na <- which(!is.na(pairs$idx_plus))
  if (length(idxp_non_na)) {
    Ap[idxp_non_na, ] <- A[pairs$idx_plus[idxp_non_na], , drop = FALSE]
  }

  idxm_non_na <- which(!is.na(pairs$idx_minus))
  if (length(idxm_non_na)) {
    Am[idxm_non_na, ] <- A[pairs$idx_minus[idxm_non_na], , drop = FALSE]
  }

  ## 5) Merge per CpG and sample -------------------------------------------
  MePer_merged <- (Ap + Am) / 2

  only_minus <- is.na(Ap) & !is.na(Am)
  only_plus  <- !is.na(Ap) & is.na(Am)

  MePer_merged[only_minus] <- Am[only_minus]
  MePer_merged[only_plus]  <- Ap[only_plus]
  # both NA -> remain NA

  rownames(MePer_merged) <- paste0(pairs$Chr, "-", pairs$CpG)
  colnames(MePer_merged) <- colnames(A)

  ## 6) Build new rowData and SE -------------------------------------------
  new_rowData <- S4Vectors::DataFrame(
    Chr    = pairs$Chr,
    Locus  = pairs$CpG,              # unified CpG coordinate in plus-strand notation
    Strand = rep("merged", n_cpg),
    Source = Source                  # "paired", "plus_only", "minus_only"
  )
  rownames(new_rowData) <- rownames(MePer_merged)

  new_assays <- S4Vectors::SimpleList(MePer = MePer_merged)

  new_se <- SummarizedExperiment::SummarizedExperiment(
    assays  = new_assays,
    rowData = new_rowData,
    colData = SummarizedExperiment::colData(se)
  )

  ## 7) Save as HDF5SummarizedExperiment -----------------------------------
  out_dir <- file.path(output_dir, output_name)
  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }

  HDF5Array::saveHDF5SummarizedExperiment(new_se, dir = out_dir, replace = TRUE)
  message("Merged CpG SE (MePer only, plus-strand notation, all CpGs) saved to: ", out_dir)

  invisible(new_se)
}
