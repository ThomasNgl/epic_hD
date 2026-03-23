merge_cpg_strands_dge <- function(dge_input,
                                  output_rds = NULL,
                                  chr_col   = "Chr",
                                  locus_col = "Locus",
                                  strand_col = "Strand") {
  # ---- 1) Load DGEList ----
  dge <- if (is.character(dge_input) && length(dge_input) == 1L) {
    if (!file.exists(dge_input)) stop("dge_input file does not exist: ", dge_input)
    readRDS(dge_input)
  } else {
    dge_input
  }

  if (is.null(dge$counts) || is.null(dge$genes)) {
    stop("Input must be a DGEList-like object with $counts and $genes.")
  }

  genes <- as.data.frame(dge$genes)
  if (!all(c(chr_col, locus_col, strand_col) %in% colnames(genes))) {
    stop("dge$genes must contain columns: ", paste(c(chr_col, locus_col, strand_col), collapse = ", "))
  }

  genes[[chr_col]]   <- as.character(genes[[chr_col]])
  genes[[locus_col]] <- as.integer(genes[[locus_col]])

  # ---- 2) Build plus / minus tables with unified CpG coordinate (plus notation) ----
  strand <- as.character(genes[[strand_col]])
  plus_idx  <- which(strand == "+")
  minus_idx <- which(strand == "-")

  if (!length(plus_idx) && !length(minus_idx)) {
    stop("No '+' or '-' strand rows found in dge$genes[[\"", strand_col, "\"]].")
  }

  plus_df <- data.frame(
    Chr      = genes[[chr_col]][plus_idx],
    CpG      = genes[[locus_col]][plus_idx],      # plus-strand CpG coordinate
    idx_plus = plus_idx,
    stringsAsFactors = FALSE
  )

  minus_df <- data.frame(
    Chr       = genes[[chr_col]][minus_idx],
    CpG       = genes[[locus_col]][minus_idx] - 1L,  # convert to plus coordinate
    idx_minus = minus_idx,
    stringsAsFactors = FALSE
  )

  # Optional safety: no duplicates per (Chr,CpG, strand)
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

  # ---- 3) Full outer join to keep all CpGs (paired + plus_only + minus_only) ----
  pairs <- merge(plus_df, minus_df, by = c("Chr", "CpG"), all = TRUE)
  if (!nrow(pairs)) stop("No CpGs found after merging plus/minus tables.")

  has_plus  <- !is.na(pairs$idx_plus)
  has_minus <- !is.na(pairs$idx_minus)

  Source <- ifelse(
    has_plus & has_minus, "paired",
    ifelse(has_plus, "plus_only", "minus_only")
  )

  # ---- 4) Build new counts by summing Me/Un across strands ----
  counts <- dge$counts
  n_cpg  <- nrow(pairs)
  n_col  <- ncol(counts)

  # Start with all zeros; we'll add plus and minus contributions
  new_counts <- matrix(0L, nrow = n_cpg, ncol = n_col,
                       dimnames = list(NULL, colnames(counts)))

  # Add plus-strand counts where present
  idxp_non_na <- which(has_plus)
  if (length(idxp_non_na)) {
    new_counts[idxp_non_na, ] <- new_counts[idxp_non_na, ] +
      counts[pairs$idx_plus[idxp_non_na], , drop = FALSE]
  }

  # Add minus-strand counts where present
  idxm_non_na <- which(has_minus)
  if (length(idxm_non_na)) {
    new_counts[idxm_non_na, ] <- new_counts[idxm_non_na, ] +
      counts[pairs$idx_minus[idxm_non_na], , drop = FALSE]
  }

  # Row names = plus-strand CpG notation
  rownames(new_counts) <- paste0(pairs$Chr, "-", pairs$CpG)

  # ---- 5) Build new genes table in plus-strand notation ----
  genes_new <- data.frame(
    Chr    = pairs$Chr,
    Locus  = pairs$CpG,        # unified CpG coordinate (plus-notation)
    Strand = "merged",
    Source = Source,
    stringsAsFactors = FALSE
  )
  rownames(genes_new) <- rownames(new_counts)

  # ---- 6) Assemble new DGEList and recompute lib.size ----
  dge_merged <- dge
  dge_merged$counts <- new_counts
  dge_merged$genes  <- genes_new

  # Recompute library sizes
  if (!is.null(dge_merged$samples)) {
    dge_merged$samples$lib.size <- colSums(new_counts)
  }

  # ---- 7) Save if requested ----
  if (!is.null(output_rds)) {
    dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(dge_merged, output_rds)
  }

  invisible(dge_merged)
}
