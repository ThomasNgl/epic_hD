
#' Prepare and clean a DGEList for methylation analysis
#'
#' Loads a DGEList from an RDS file (or uses a provided DGEList), removes
#' all-zero rows, optionally removes specified samples, then filters rows by:
#' (i) minimum total reads per group, (ii) maximum total reads per group,
#' and (iii) minimum number of samples per group with >=1 read.
#'
#' Totals are computed on coverage = methylated + unmethylated (if `unmethylated`
#' exists as a matrix in the object).
#'
#' @param me_dge_path Either a DGEList object or a character path to an .rds file.
#' @param samples_to_remove Character vector of sample IDs to remove (matched to
#'   `y$samples$sample_clean` if present; otherwise matched by column name patterns).
#' @param chr2rm Character vector of chromosomes to remove (e.g. c("chrM","2") or c("chr2","chrM")).
#' @param min_reads Numeric. Keep rows where each group has >= min_reads total reads
#'   across its samples. Set NULL to skip.
#' @param max_reads Numeric. Keep rows where each group has <= max_reads total reads
#'   across its samples. Set NULL to skip.
#' @param min_samples Integer. Keep rows where each group has >= min_samples samples
#'   with >=1 read. Set NULL to skip.
#' @param output_dir Character. Directory to save the cleaned object.
#' @param output_name Character. Filename for the saved .rds.
#' @param groups Character length-2. Group names to enforce filters on
#'   (default c("Control","MSUS")).
#' @return The cleaned DGEList (invisibly). Also saved to disk.
#' @export
clean_dge <- function(me_dge_path,
                                     samples_to_remove = NULL,
                                     chr2rm = NULL,
                                     min_reads = NULL,
                                     max_reads = NULL,
                                     min_samples = NULL,
                                     output_dir,
                                     output_name,
                                     groups = c("Control", "MSUS"),
                                     chr_col = 'Chr') {

  .msg_step <- function(label, before, after) {
    message(label, ": kept ", after, " / ", before, " (removed ", before - after, ")")
  }

  # Coverage per COLUMN (for all-zero row removal)
  .total_cov_colwise <- function(y) {
    if (!is.null(y$unmethylated)) {
      y$counts + y$unmethylated
    } else {
      y$counts
    }
  }

  # Coverage per BIOLOGICAL sample (Me + Un), and BIO group indices
  .bio_cov_and_idx <- function(y, groups) {

    # Case A: separate matrices (counts + unmethylated) already represent biological samples
    if (!is.null(y$unmethylated)) {
      cov_bio <- y$counts + y$unmethylated

      if (!("Group" %in% colnames(y$samples))) stop("y$samples must contain 'Group'.")
      bio_group <- as.character(y$samples$Group)
      names(bio_group) <- colnames(y$counts)

      idx <- lapply(groups, function(g) which(bio_group == g))
      names(idx) <- groups

      return(list(cov_bio = cov_bio, bio_group = bio_group, idx = idx))
    }

    # Case B: single matrix with -Me/-Un columns (2 columns per biological sample)
    cn <- colnames(y$counts)
    me_cols <- grep("-Me$", cn, value = TRUE)
    un_cols <- grep("-Un$", cn, value = TRUE)

    if (length(me_cols) == 0 || length(un_cols) == 0) {
      stop("No y$unmethylated slot and could not find both '-Me' and '-Un' columns in y$counts.")
    }

    bases <- sort(unique(c(sub("-Me$", "", me_cols), sub("-Un$", "", un_cols))))

    me_map <- setNames(me_cols, sub("-Me$", "", me_cols))
    un_map <- setNames(un_cols, sub("-Un$", "", un_cols))

    cov_bio <- matrix(0,
                      nrow = nrow(y$counts),
                      ncol = length(bases),
                      dimnames = list(rownames(y$counts), bases))

    for (b in bases) {
      if (!is.na(me_map[b])) cov_bio[, b] <- cov_bio[, b] + y$counts[, me_map[b]]
      if (!is.na(un_map[b])) cov_bio[, b] <- cov_bio[, b] + y$counts[, un_map[b]]
    }

    # In your object, y$samples has a 'Sample' column = base sample id (same for -Me/-Un)
    smp <- as.data.frame(y$samples)
    smp$colname <- rownames(smp)

    if ("Sample" %in% colnames(smp)) {
      smp$base <- as.character(smp$Sample)
    } else {
      smp$base <- sub("-(Me|Un)$", "", smp$colname)
    }

    if (!("Group" %in% colnames(smp))) stop("y$samples must contain 'Group'.")
    smp$Group <- as.character(smp$Group)

    grp_by_base <- dplyr::distinct(smp, base, Group)

    # Each biological sample must map to exactly one group
    if (any(duplicated(grp_by_base$base))) {
      stop("A biological sample maps to multiple Group values in y$samples (check metadata consistency).")
    }

    bio_group <- setNames(grp_by_base$Group, grp_by_base$base)
    bio_group <- bio_group[bases]  # align order to cov_bio columns

    idx <- lapply(groups, function(g) which(bio_group == g))
    names(idx) <- groups

    list(cov_bio = cov_bio, bio_group = bio_group, idx = idx)
  }
  ## ----------------------------
  ## Reporting
  ## ----------------------------
  report_prefix <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(output_name), "_report")
  )
  report_summary_path  <- paste0(report_prefix, "_summary.txt")
  report_features_path <- paste0(report_prefix, "_removed_features.tsv")

  .report_has_header <- FALSE

  dir.create(dirname(report_summary_path), recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    paste0("clean_dge report: ", format(Sys.time())),
    paste0("Input: ", if (is.character(me_dge_path)) me_dge_path else "<DGEList object>"),
    paste0("Output: ", file.path(output_dir, output_name)),
    ""
  ), con = report_summary_path)

  .report_summary <- function(line) {
    write(line, file = report_summary_path, append = TRUE)
  }

  # Append removed features for one filtering step.
  # Writes biological-sample coverage (Me+Un) per removed feature.
  .report_removed_features <- function(step, y_obj, rm_mask) {
    if (!any(rm_mask)) return(invisible(NULL))

    bio <- .bio_cov_and_idx(y_obj, groups = groups)
    cov_bio <- bio$cov_bio

    cov_rm <- cov_bio[rm_mask, , drop = FALSE]
    feat_id <- rownames(y_obj$counts)[rm_mask]

    out <- as.data.frame(cov_rm, check.names = FALSE)
    out <- cbind(Step = step, Feature = feat_id, out)

    if (!is.null(y_obj$genes)) {
      g <- as.data.frame(y_obj$genes[rm_mask, , drop = FALSE])
      keep_cols <- intersect(c(chr_col, "Chr", "Locus", "Strand"), colnames(g))
      if (length(keep_cols)) {
        out <- cbind(out[, 1:2, drop = FALSE], g[, keep_cols, drop = FALSE], out[, -(1:2), drop = FALSE])
      }
    }

    write.table(
      out,
      file = report_features_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = !.report_has_header,
      append = .report_has_header
    )
    .report_has_header <<- TRUE

    invisible(NULL)
  }



  # --- 0) load / use DGEList ---
  if (inherits(me_dge_path, "DGEList")) {
    y <- me_dge_path
    message("Using provided DGEList object.")
  } else if (is.character(me_dge_path) && length(me_dge_path) == 1) {
    if (!file.exists(me_dge_path)) stop("File not found: ", me_dge_path)
    y <- readRDS(me_dge_path)
    message("Loaded DGEList from: ", me_dge_path, " (regions=", nrow(y$counts),
            ", samples=", ncol(y$counts), ")")
  } else {
    stop("'me_dge_path' must be a DGEList or a path to an RDS file.")
  }

  if (is.null(y$samples) || !("Group" %in% colnames(y$samples))) {
    stop("y$samples must contain a 'Group' column.")
  }

# --- 0b) remove chromosomes listed in chr2rm ---
if (!is.null(chr2rm) && length(chr2rm) > 0) {
  stopifnot(!is.null(y$genes), chr_col %in% colnames(y$genes))
  before <- nrow(y$counts)

  chr_rm <- ifelse(grepl("^chr", chr2rm), chr2rm, paste0("chr", chr2rm))

  gene_chr <- as.character(y$genes[[chr_col]])
  keep <- !(gene_chr %in% chr_rm)

  .report_summary(paste0("Removed chromosomes: ", paste(chr_rm, collapse = ",")))
  .report_summary(paste0("Step=chr_remove; removing rows=", sum(!keep)))
  .report_removed_features("chr_remove", y, rm_mask = !keep)
  
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (!is.null(y$unmethylated)) y$unmethylated <- y$unmethylated[keep, , drop = FALSE]

  .msg_step(paste0("Remove chromosomes (", chr_col, "): ", paste(chr_rm, collapse = ",")),
            before, nrow(y$counts))
}


  # --- 1) remove all-zero rows (using coverage per COLUMN) ---
  before <- nrow(y$counts)
  cov_col <- .total_cov_colwise(y)
  keep <- rowSums(cov_col) > 0
  .report_summary(paste0("Step=all_zero; removing rows=", sum(!keep)))
  .report_removed_features("all_zero", y, rm_mask = !keep)
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (!is.null(y$unmethylated)) y$unmethylated <- y$unmethylated[keep, , drop = FALSE]
  .msg_step("Filter all-zero rows", before, nrow(y$counts))

  # --- 2) remove specified samples (remove both -Me and -Un via y$samples$Sample when available) ---
  if (!is.null(samples_to_remove) && length(samples_to_remove) > 0) {
    before <- ncol(y$counts)

    if ("Sample" %in% colnames(y$samples)) {
      keep_cols <- !(as.character(y$samples$Sample) %in% samples_to_remove)
    } else if ("sample_clean" %in% colnames(y$samples)) {
      keep_cols <- !(y$samples$sample_clean %in% samples_to_remove)
    } else {
      pattern <- paste0("(", paste(samples_to_remove, collapse="|"), ")(-Me|-Un)?$")
      keep_cols <- !grepl(pattern, colnames(y$counts))
    }

    removed_cols <- colnames(y$counts)[!keep_cols]
    .report_summary(paste0("Removed columns (technical): ", paste(removed_cols, collapse = ",")))

    if ("Sample" %in% colnames(y$samples)) {
      removed_bio <- unique(as.character(y$samples$Sample[!keep_cols]))
      .report_summary(paste0("Removed samples (biological): ", paste(removed_bio, collapse = ",")))
    }

    y <- y[, keep_cols, keep.lib.sizes = FALSE]
    if (!is.null(y$unmethylated)) y$unmethylated <- y$unmethylated[, keep_cols, drop = FALSE]
    .msg_step("Remove samples", before, ncol(y$counts))
  } else {
    message("Remove samples: skipped (samples_to_remove is NULL/empty).")
  }

  # --- Build BIO coverage + BIO group indices (after sample removal) ---
  bio_info <- .bio_cov_and_idx(y, groups)
  cov_bio <- bio_info$cov_bio
  idx_bio <- bio_info$idx

  if (any(vapply(idx_bio, length, 0L) == 0L)) {
    stop("Missing group samples after filtering. Biological-sample counts per group: ",
         paste(names(idx_bio), vapply(idx_bio, length, 0L), sep="=", collapse=", "))
  }

  # --- 3a) min_reads per group (BIOLOGICAL) ---
  if (!is.null(min_reads)) {
    before <- nrow(y$counts)

    g1 <- rowSums(cov_bio[, idx_bio[[groups[1]]], drop = FALSE])
    g2 <- rowSums(cov_bio[, idx_bio[[groups[2]]], drop = FALSE])
    keep <- (g1 >= min_reads) & (g2 >= min_reads)

    .report_summary(paste0("Step=min_reads_", min_reads, "; removing rows=", sum(!keep)))
    .report_removed_features(paste0("min_reads_", min_reads), y, rm_mask = !keep)
    
    y <- y[keep, , keep.lib.sizes = FALSE]
    if (!is.null(y$unmethylated)) y$unmethylated <- y$unmethylated[keep, , drop = FALSE]
    .msg_step(paste0("Filter min_reads (", min_reads, ") [BIO]"), before, nrow(y$counts))

    # refresh bio matrices after row filtering
    bio_info <- .bio_cov_and_idx(y, groups); cov_bio <- bio_info$cov_bio; idx_bio <- bio_info$idx
  } else message("Filter min_reads: skipped (min_reads is NULL).")

  # --- 3b) max_reads per group (BIOLOGICAL) ---
  if (!is.null(max_reads)) {
    before <- nrow(y$counts)

    g1 <- rowSums(cov_bio[, idx_bio[[groups[1]]], drop = FALSE])
    g2 <- rowSums(cov_bio[, idx_bio[[groups[2]]], drop = FALSE])
    keep <- (g1 <= max_reads) & (g2 <= max_reads)
    
    .report_summary(paste0("Step=max_reads_", max_reads, "; removing rows=", sum(!keep)))
    .report_removed_features(paste0("max_reads_", max_reads), y, rm_mask = !keep)
    
    y <- y[keep, , keep.lib.sizes = FALSE]
    if (!is.null(y$unmethylated)) y$unmethylated <- y$unmethylated[keep, , drop = FALSE]
    .msg_step(paste0("Filter max_reads (", max_reads, ") [BIO]"), before, nrow(y$counts))

    bio_info <- .bio_cov_and_idx(y, groups); cov_bio <- bio_info$cov_bio; idx_bio <- bio_info$idx
  } else message("Filter max_reads: skipped (max_reads is NULL).")

  # --- 3c) min_samples with >=1 read per group (BIOLOGICAL) ---
  if (!is.null(min_samples)) {
    before <- nrow(y$counts)

    k1 <- rowSums(cov_bio[, idx_bio[[groups[1]]], drop = FALSE] > 0)
    k2 <- rowSums(cov_bio[, idx_bio[[groups[2]]], drop = FALSE] > 0)
    keep <- (k1 >= min_samples) & (k2 >= min_samples)

    .report_summary(paste0("Step=min_samples_", min_samples, "; removing rows=", sum(!keep)))
    .report_removed_features(paste0("min_samples_", min_samples), y, rm_mask = !keep)

    y <- y[keep, , keep.lib.sizes = FALSE]
    if (!is.null(y$unmethylated)) y$unmethylated <- y$unmethylated[keep, , drop = FALSE]
    .msg_step(paste0("Filter min_samples (", min_samples, ") [BIO]"), before, nrow(y$counts))

    bio_info <- .bio_cov_and_idx(y, groups); cov_bio <- bio_info$cov_bio; idx_bio <- bio_info$idx
  } else message("Filter min_samples: skipped (min_samples is NULL).")

  message("Final DGEList: regions=", nrow(y$counts), ", samples=", ncol(y$counts))

  # --- 4) save ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_path <- file.path(output_dir, output_name)
  saveRDS(y, output_path)
  message("Saved cleaned DGEList to: ", output_path)

  invisible(y)
}

# usage:
# out <- sum_meun_by_meta(dge_clean, "Female_id")
# out$Me  # features x Female_id
# out$Un
