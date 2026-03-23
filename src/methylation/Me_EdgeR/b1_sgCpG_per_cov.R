sum_MeUn_by_cov <- function(dge, cov = "Female_id", suffix_re = "-(Me|Un)$") {
  stopifnot(inherits(dge, "DGEList"))
  if (!cov %in% colnames(dge$samples)) stop("cov not found in dge$samples: ", cov)

  cn <- colnames(dge$counts)

  # Parse status from column names (expects -Me / -Un by default)
  status <- sub(paste0("^.*", suffix_re), "\\1", cn)
  if (any(status == cn)) stop("Could not parse Me/Un from colnames using suffix_re = ", suffix_re)

  cov_val <- as.character(dge$samples[[cov]])
  key <- paste(cov_val, status, sep = "-")  # e.g. "14-Me"

  # Sum counts by key
  idx_list <- split(seq_along(key), key)
  new_counts <- matrix(0L, nrow = nrow(dge$counts), ncol = length(idx_list),
                       dimnames = list(rownames(dge$counts), names(idx_list)))
  for (j in seq_along(idx_list)) {
    ii <- idx_list[[j]]
    new_counts[, j] <- if (length(ii) == 1L) dge$counts[, ii] else rowSums(dge$counts[, ii, drop = FALSE])
  }

  # Build samples table
  parts <- do.call(rbind, strsplit(colnames(new_counts), "-", fixed = TRUE))
  new_samples <- data.frame(
    row.names    = colnames(new_counts),
    MeUn         = parts[, 2],
    stringsAsFactors = FALSE
  )
  new_samples[[cov]] <- suppressWarnings(type.convert(parts[, 1], as.is = TRUE))
  new_samples$lib.size <- colSums(new_counts)
  new_samples$norm.factors <- 1

  # Carry Group if consistent within cov
  if ("Group" %in% colnames(dge$samples)) {
    cov_to_group <- tapply(as.character(dge$samples$Group), cov_val, function(x) {
      ux <- unique(x); if (length(ux) == 1) ux else NA_character_
    })
    new_samples$Group <- unname(cov_to_group[as.character(new_samples[[cov]])])
  }

  edgeR::DGEList(counts = new_counts, samples = new_samples, genes = dge$genes)
}


sum_MeUn_by_cov_save <- function(dge_or_path,
                                cov = "Female_id",
                                output_dir,
                                output_name = NULL,
                                suffix_re = "-(Me|Un)$",
                                recalc_norm = TRUE,
                                overwrite = TRUE) {

  # Load if a path was provided
  if (is.character(dge_or_path) && length(dge_or_path) == 1L) {
    if (!file.exists(dge_or_path)) stop("File does not exist: ", dge_or_path)
    dge <- readRDS(dge_or_path)

    if (is.null(output_name)) {
      base <- tools::file_path_sans_ext(basename(dge_or_path))
      output_name <- paste0(base, "_sumBy_", cov, ".rds")
    }
  } else {
    dge <- dge_or_path
    if (is.null(output_name)) output_name <- paste0("dge_sumBy_", cov, ".rds")
  }

  if (!inherits(dge, "DGEList")) stop("Input must be a DGEList or a path to an RDS containing a DGEList.")

  out <- sum_MeUn_by_cov(dge, cov = cov, suffix_re = suffix_re)
  if (recalc_norm) out <- edgeR::calcNormFactors(out)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(output_dir, output_name)

  if (file.exists(out_path) && !overwrite) stop("Output exists and overwrite=FALSE: ", out_path)
  saveRDS(out, out_path)

  invisible(list(dge = out, path = out_path))
}
