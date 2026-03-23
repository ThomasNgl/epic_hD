suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(HDF5Array)
  library(S4Vectors)
})

make_strict_mepe_and_dge <- function(se_in_dir,
                                     dge_in_rds,
                                     out_dir,
                                     meper_assay_name = "MePer",
                                     add_M_from_dge = TRUE,
                                     se_out_stub = "mepe_se_sgCpGs_strict",
                                     dge_out_name = "yall_me_dge_sgCpGs_strict.rds",
                                     se_plain_rds_path = NULL,
                                     m_pseudocount = 2) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## ----------------------------
  ## Helpers: DGE coverage per biological sample (Me+Un)
  ## ----------------------------
  .dge_cov_bio <- function(dge) {
    cn <- colnames(dge$counts)
    me_cols <- grep("-Me$", cn, value = TRUE)
    un_cols <- sub("-Me$", "-Un", me_cols)
    ok <- un_cols %in% cn
    me_cols <- me_cols[ok]
    un_cols <- un_cols[ok]
    if (!length(me_cols)) stop("No complete -Me/-Un pairs found in dge$counts.")
    
    # map by base sample name
    base <- sub("-Me$", "", me_cols)
    me_map <- setNames(me_cols, base)
    un_map <- setNames(un_cols, base)
    
    cov <- dge$counts[, me_map, drop = FALSE] + dge$counts[, un_map, drop = FALSE]
    colnames(cov) <- base
    cov
  }
  
  ## ----------------------------
  ## Reporting init
  ## ----------------------------
  report_prefix <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(dge_out_name), "_report")
  )
  report_summary_path  <- paste0(report_prefix, "_summary.txt")
  report_removed_tsv <- paste0(report_prefix, "_removed_features.tsv")
  .report_has_header <- FALSE
  
  writeLines(c(
    paste0("strict pipeline report: ", format(Sys.time())),
    paste0("SE input dir: ", se_in_dir),
    paste0("DGE input rds: ", dge_in_rds),
    paste0("Output dir: ", out_dir),
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
  
  ## ----------------------------
  ## Load inputs
  ## ----------------------------
  se_in <- HDF5Array::loadHDF5SummarizedExperiment(se_in_dir)
  dge_in <- readRDS(dge_in_rds)
  
  .report_summary(paste0("SE dims: ", nrow(se_in), " x ", ncol(se_in)))
  .report_summary(paste0("DGE dims: ", nrow(dge_in$counts), " x ", ncol(dge_in$counts)))
  
  ## ----------------------------
  ## Step 1: strict SE = no NA in MePer
  ## ----------------------------
  A_meper <- SummarizedExperiment::assay(se_in, meper_assay_name)
  keep_se <- rowSums(is.na(A_meper)) == 0
  rm_se_ids <- rownames(se_in)[!keep_se]
  
  .report_summary(paste0("Step=se_noNA; removed rows=", length(rm_se_ids)))
  .append_removed("se_noNA", dge_in, rm_se_ids)
  
  se_strict <- se_in[keep_se, ]
  
  ## ----------------------------
  ## Step 2: strict DGE = coverage >=1 in every biological sample
  ## ----------------------------
  cov_bio <- .dge_cov_bio(dge_in)
  n_bio <- ncol(cov_bio)
  keep_dge_cov <- rowSums(cov_bio >= 1) == n_bio
  rm_dge_cov_ids <- rownames(dge_in$counts)[!keep_dge_cov]
  
  .report_summary(paste0("Step=dge_allSamplesCov; removed rows=", length(rm_dge_cov_ids)))
  .append_removed("dge_allSamplesCov", dge_in, rm_dge_cov_ids)
  
  dge_cov_strict <- dge_in[keep_dge_cov, , keep.lib.sizes = FALSE]
  
  ## ----------------------------
  ## Step 3: align SE and DGE to the intersection of remaining CpGs
  ## ----------------------------
  ids_se  <- rownames(se_strict)
  ids_dge <- rownames(dge_cov_strict$counts)
  ids_keep <- intersect(ids_se, ids_dge)
  
  rm_align_from_se  <- setdiff(ids_se,  ids_keep)
  rm_align_from_dge <- setdiff(ids_dge, ids_keep)
  
  .report_summary(paste0("Step=align_intersect; removing from SE=", length(rm_align_from_se),
                         "; removing from DGE=", length(rm_align_from_dge)))
  .append_removed("align_drop_from_se",  dge_cov_strict, rm_align_from_se)
  .append_removed("align_drop_from_dge", dge_cov_strict, rm_align_from_dge)
  
  se_strict <- se_strict[ids_keep, ]
  dge_strict <- dge_cov_strict[ids_keep, , keep.lib.sizes = FALSE]
  
  ## ----------------------------
  ## Optional: add M_from_dge assay to SE
  ## ----------------------------
  if (add_M_from_dge) {
    cn <- colnames(dge_strict$counts)
    me_cols <- grep("-Me$", cn, value = TRUE)
    un_cols <- sub("-Me$", "-Un", me_cols)
    ok <- un_cols %in% cn
    me_cols <- me_cols[ok]
    un_cols <- un_cols[ok]
    base <- sub("-Me$", "", me_cols)
    
    Me <- dge_strict$counts[, me_cols, drop = FALSE]
    Un <- dge_strict$counts[, un_cols, drop = FALSE]
    M_from_dge <- log2(Me + m_pseudocount) - log2(Un + m_pseudocount)
    colnames(M_from_dge) <- base
    
    # align columns to SE colnames
    M_from_dge <- M_from_dge[, colnames(se_strict), drop = FALSE]
    
    SummarizedExperiment::assay(se_strict, "M_from_dge") <- M_from_dge
    .report_summary("Added assay: M_from_dge")
  }
  
  ## ----------------------------
  ## Save outputs
  ## ----------------------------
  se_out_dir <- file.path(out_dir, se_out_stub)
  if (dir.exists(se_out_dir)) unlink(se_out_dir, recursive = TRUE, force = TRUE)
  
  HDF5Array::saveHDF5SummarizedExperiment(se_strict, dir = se_out_dir, replace = TRUE)
  .report_summary(paste0("Saved HDF5SummarizedExperiment: ", se_out_dir))
  
  dge_out_path <- file.path(out_dir, dge_out_name)
  saveRDS(dge_strict, dge_out_path)
  .report_summary(paste0("Saved strict DGE: ", dge_out_path))
  
  # Optional: save plain RDS of SE (fully realized in memory)
  if (!is.null(se_plain_rds_path)) {
    dir.create(dirname(se_plain_rds_path), recursive = TRUE, showWarnings = FALSE)
    
    assay_list_inmem <- lapply(SummarizedExperiment::assays(se_strict), function(a) as.matrix(a))
    se_inmem <- SummarizedExperiment(
      assays  = assay_list_inmem,
      rowData = S4Vectors::DataFrame(as.data.frame(SummarizedExperiment::rowData(se_strict))),
      colData = S4Vectors::DataFrame(as.data.frame(SummarizedExperiment::colData(se_strict)))
    )
    
    saveRDS(se_inmem, se_plain_rds_path)
    .report_summary(paste0("Saved plain SE RDS: ", se_plain_rds_path))
  }
  
  invisible(list(
    se_out_dir = se_out_dir,
    dge_out_path = dge_out_path,
    report_summary = report_summary_path,
    report_removed_features = report_removed_tsv
  ))
}
