zero_counts_where_mepe_NA <- function(mepe,
                                      dge,
                                      assay_name  = "MePer",
                                      report_path = NULL) {
  stopifnot(inherits(mepe, "SummarizedExperiment"))

  # 1) Extract assay and make sure it's a matrix
  A <- SummarizedExperiment::assay(mepe, assay_name)
  if (!is.matrix(A)) A <- as.matrix(A)

  # 2) Sanity check: loci alignment between mepe and dge
  rd_mepe <- as.data.frame(SummarizedExperiment::rowData(mepe))
  rd_dge  <- as.data.frame(dge$genes)

  if (!all(c("Chr", "Locus") %in% colnames(rd_mepe)) ||
      !all(c("Chr", "Locus") %in% colnames(rd_dge))) {
    stop("Both mepe$rowData and dge$genes must have Chr and Locus.")
  }

  if (nrow(rd_mepe) != nrow(rd_dge) ||
      any(rd_mepe$Chr != rd_dge$Chr) ||
      any(rd_mepe$Locus != rd_dge$Locus)) {
    stop("Row order / loci do not match between mepe and dge.")
  }

  # 3) Loop over samples (columns in mepe)
  sample_names <- colnames(A)
  report <- data.frame(
    Sample                      = sample_names,
    N_NA                        = integer(length(sample_names)),
    N_sites_with_NA             = integer(length(sample_names)),
    N_Me_zeroed_from_nonzero    = integer(length(sample_names)),
    N_Un_zeroed_from_nonzero    = integer(length(sample_names)),
    stringsAsFactors            = FALSE
  )

  for (j in seq_along(sample_names)) {
    s <- sample_names[j]
    mask_NA <- is.na(A[, j])

    me_col  <- paste0(s, "-Me")
    un_col  <- paste0(s, "-Un")

    if (!all(c(me_col, un_col) %in% colnames(dge$counts))) {
      stop("Missing Me/Un columns in dge$counts for sample: ", s)
    }

    old_me <- dge$counts[, me_col]
    old_un <- dge$counts[, un_col]

    report$N_NA[j]            <- sum(mask_NA)
    report$N_sites_with_NA[j] <- sum(mask_NA)
    report$N_Me_zeroed_from_nonzero[j] <- sum(mask_NA & old_me != 0)
    report$N_Un_zeroed_from_nonzero[j] <- sum(mask_NA & old_un != 0)

    # Set counts to 0 where mepe assay is NA
    dge$counts[mask_NA, me_col] <- 0L
    dge$counts[mask_NA, un_col] <- 0L
  }

  # 4) Optional: recompute library sizes to reflect the zeroed counts
  dge$samples$lib.size <- colSums(dge$counts)

  # 5) Optionally write report to CSV
  if (!is.null(report_path)) {
    dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(report, report_path, row.names = FALSE)
  }

  list(
    dge    = dge,    # updated DGEList
    report = report  # per-sample summary of what was zeroed
  )
}
