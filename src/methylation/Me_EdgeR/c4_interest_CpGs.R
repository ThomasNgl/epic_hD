suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(HDF5Array)
  library(S4Vectors)
  library(edgeR)
})

remove_cst_mepe_and_dge <- function(se_in_dir,
                                   dge_in_rds,
                                   out_dir,
                                   meper_assay_name = "MePer",
                                   tol = 0.1,
                                   chunk_nrows = 200000,
                                   keep_rows_with_na = TRUE,
                                   se_out_stub = "mepe_se_sgCpGs_interest",
                                   dge_out_name = "yall_me_dge_sgCpGs_interest.rds") {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ## ---- helpers ----
  .dge_cov_bio <- function(dge) {
    cn <- colnames(dge$counts)
    me_cols <- grep("-Me$", cn, value = TRUE)
    un_cols <- sub("-Me$", "-Un", me_cols)
    ok <- un_cols %in% cn
    me_cols <- me_cols[ok]; un_cols <- un_cols[ok]
    if (!length(me_cols)) stop("No complete -Me/-Un pairs found in dge$counts.")

    base <- sub("-Me$", "", me_cols)
    me_map <- setNames(me_cols, base)
    un_map <- setNames(un_cols, base)

    cov <- dge$counts[, me_map, drop = FALSE] + dge$counts[, un_map, drop = FALSE]
    colnames(cov) <- base
    cov
  }

  ## ---- reporting init ----
  report_prefix <- file.path(out_dir, paste0(tools::file_path_sans_ext(dge_out_name), "_report"))
  report_summary_path <- paste0(report_prefix, "_summary.txt")
  report_removed_tsv  <- paste0(report_prefix, "_removed_features.tsv")
  .report_has_header <- FALSE

  writeLines(c(
    paste0("interest pipeline report: ", format(Sys.time())),
    paste0("SE input dir: ", se_in_dir),
    paste0("DGE input rds: ", dge_in_rds),
    paste0("MePer assay: ", meper_assay_name),
    paste0("Constant tolerance (max-min <= tol): ", tol),
    paste0("Chunk rows: ", chunk_nrows),
    ""
  ), con = report_summary_path)

  .report_summary <- function(line) write(line, file = report_summary_path, append = TRUE)

  .append_removed <- function(step, dge_obj, rm_ids) {
    if (!length(rm_ids)) return(invisible(NULL))

    cov_bio <- .dge_cov_bio(dge_obj)
    rm_mask <- rownames(dge_obj$counts) %in% rm_ids
    if (!any(rm_mask)) return(invisible(NULL))

    out <- as.data.frame(cov_bio[rm_mask, , drop = FALSE], check.names = FALSE)
    out <- cbind(Step = step, Feature = rownames(dge_obj$counts)[rm_mask], out)

    if (!is.null(dge_obj$genes)) {
      gdf <- as.data.frame(dge_obj$genes[rm_mask, , drop = FALSE])
      keep_cols <- intersect(c("Chr", "Locus", "Strand"), colnames(gdf))
      if (length(keep_cols)) {
        out <- cbind(out[, 1:2, drop = FALSE], gdf[, keep_cols, drop = FALSE], out[, -(1:2), drop = FALSE])
      }
    }

    write.table(
      out, file = report_removed_tsv,
      sep = "\t", quote = FALSE, row.names = FALSE,
      col.names = !.report_has_header,
      append = .report_has_header
    )
    .report_has_header <<- TRUE
    invisible(NULL)
  }

  ## ---- load inputs ----
  se_in  <- HDF5Array::loadHDF5SummarizedExperiment(se_in_dir)
  dge_in <- readRDS(dge_in_rds)

  .report_summary(paste0("SE dims: ", nrow(se_in), " x ", ncol(se_in)))
  .report_summary(paste0("DGE dims: ", nrow(dge_in$counts), " x ", ncol(dge_in$counts)))

  ## ---- Step 1: remove (near-)constant MePer rows ----
  A <- SummarizedExperiment::assay(se_in, meper_assay_name)
  n <- nrow(A)
  keep <- rep(TRUE, n)

  for (start in seq.int(1L, n, by = chunk_nrows)) {
    end <- min(n, start + chunk_nrows - 1L)
    blk <- as.matrix(A[start:end, , drop = FALSE])

    if (keep_rows_with_na) {
      has_na <- rowSums(is.na(blk)) > 0
    } else {
      has_na <- rep(FALSE, nrow(blk))
    }

    rmin <- apply(blk, 1, min, na.rm = TRUE)
    rmax <- apply(blk, 1, max, na.rm = TRUE)
    rng  <- rmax - rmin

    # keep rows with NA (optional), else keep only if range > tol
    keep_blk <- has_na | (is.finite(rng) & (rng > tol))
    keep[start:end] <- keep_blk
  }

  rm_ids_const <- rownames(se_in)[!keep]
  .report_summary(paste0("Step=se_const_meper; removed rows=", length(rm_ids_const)))
  .append_removed("se_const_meper", dge_in, rm_ids_const)

  se_var <- se_in[keep, ]

  ## ---- Step 2: align SE and DGE to intersection of remaining CpGs ----
  ids_se  <- rownames(se_var)
  ids_dge <- rownames(dge_in$counts)
  ids_keep <- intersect(ids_se, ids_dge)

  rm_from_se  <- setdiff(ids_se, ids_keep)
  rm_from_dge <- setdiff(ids_dge, ids_keep)

  .report_summary(paste0("Step=align_intersect; drop_from_se=", length(rm_from_se),
                        "; drop_from_dge=", length(rm_from_dge)))
  .append_removed("align_drop_from_se",  dge_in, rm_from_se)
  .append_removed("align_drop_from_dge", dge_in, rm_from_dge)

  se_out  <- se_var[ids_keep, ]
  dge_out <- dge_in[ids_keep, , keep.lib.sizes = FALSE]

  ## ---- save outputs ----
  se_dir_name <- tools::file_path_sans_ext(se_out_stub)  # tolerates callers passing ".rds"
  se_out_dir  <- file.path(out_dir, se_dir_name)
  if (dir.exists(se_out_dir)) unlink(se_out_dir, recursive = TRUE, force = TRUE)

  HDF5Array::saveHDF5SummarizedExperiment(se_out, dir = se_out_dir, replace = TRUE)
  dge_out_path <- file.path(out_dir, dge_out_name)
  saveRDS(dge_out, dge_out_path)

  .report_summary(paste0("Saved SE dir: ", se_out_dir))
  .report_summary(paste0("Saved DGE rds: ", dge_out_path))

  invisible(list(
    se_out_dir = se_out_dir,
    dge_out_path = dge_out_path,
    report_summary = report_summary_path,
    report_removed_features = report_removed_tsv
  ))
}
