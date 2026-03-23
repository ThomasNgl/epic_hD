#' Filter a methylation DGEList to CpG loci defined by per-chromosome CpG-location CSVs
#'
#' This function subsets a DGEList (e.g., from readBismark2DGE) to retain only rows whose
#' genomic locus matches CpG positions provided in `CpGslocation_chr*.csv` files.
#'
#' Strand assignment rule (relative to the CSV position list for each chromosome):
#' - "+" (plus):  DGE locus exactly matches a CSV position p
#' - "-" (minus): DGE locus matches (p + 1)
#' - "ambig":     DGE locus matches both a position p and some (p' + 1) simultaneously
#'                (rare; typically indicates overlapping lists or duplicated/shifted entries)
#'
#' Outputs are always written:
#' - A per-chromosome mapping summary CSV
#' - The filtered DGEList saved as an RDS
#'
#' @param me_dge Either:
#'   (1) a DGEList-like object with a `$genes` data.frame containing columns `Chr` and `Locus`, or
#'   (2) a character scalar path to an `.rds` file that contains such an object.
#'
#' @param CpGs_loc_dir Character scalar. Directory containing per-chromosome CpG location CSV files.
#'   Files are expected to be named like `paste0(file_prefix, chr, ".csv")`,
#'   e.g. `CpGslocation_chrM.csv`.
#'
#' @param ChrNames Character vector of chromosome names to iterate over.
#'   Default: `chr1..chr19, chrX, chrY, chrM`.
#'
#' @param file_prefix Character scalar. Prefix for CpG location CSV filenames.
#'   Default: `"CpGslocation_"`.
#'
#' @param position_col Character scalar. Name of the CpG position column in each CSV.
#'   Default: `"position"`.
#'
#' @param summary_out_csv Character scalar. Output path for the per-chromosome summary CSV.
#'
#' @param dge_out_rds Character scalar. Output path for the filtered DGEList `.rds`.
#'
#' @return The filtered DGEList (same class as input) with an added column
#'   `me_dge_cpg$genes$Strand` containing `"+"`, `"-"`, or `"ambig"`.
#'
#' @details
#' Summary CSV fields (per chromosome):
#' - n_csv_positions: number of unique CpG positions in the CSV for that chromosome
#' - n_csv_positions_mapped: number of CSV positions p for which p OR (p+1) exists in the DGE loci
#' - pct_csv_positions_mapped: percentage of CSV positions mapped
#' - n_dge_rows_chr: number of DGE rows on that chromosome
#' - n_dge_rows_mapped / n_dge_rows_unmapped: counts of mapped / unmapped DGE rows
#' - pct_dge_rows_mapped: percent of DGE rows mapped on that chromosome
#' - n_plus / n_minus / n_ambig: strand counts among mapped rows
#' - pct_plus_among_mapped / pct_minus_among_mapped: strand percentages among mapped rows
#'
#' The DGEList is subset using `keep.lib.sizes = FALSE` (library sizes recomputed), as requested.
#'
#' @examples
#' \dontrun{
#' me_dge_cpg <- get_CpGs_dge(
#'   me_dge = me_dge,
#'   CpGs_loc_dir = "/mnt/thomas/matriline/results/GRCm39/bismark_genome/CpGs_location",
#'   summary_out_csv = "/mnt/thomas/matriline/results/GRCm39/bismark_genome/mapping_summary_per_chr.csv",
#'   dge_out_rds = "/mnt/thomas/matriline/results/GRCm39/bismark_genome/me_dge_cpg.rds"
#' )
#'
#' # Or, load DGEList from an .rds path:
#' me_dge_cpg <- get_CpGs_dge(
#'   me_dge = "/path/to/me_dge.rds",
#'   CpGs_loc_dir = "/mnt/thomas/matriline/results/GRCm39/bismark_genome/CpGs_location",
#'   summary_out_csv = "/tmp/mapping_summary_per_chr.csv",
#'   dge_out_rds = "/tmp/me_dge_cpg.rds"
#' )
#' }
get_CpGs_dge <- function(me_dge,
                        CpGs_loc_dir,
                        ChrNames = paste0("chr", c(1:19, "X", "Y", "M")),
                        file_prefix = "CpGslocation_",
                        position_col = "position",
                        summary_out_csv,
                        dge_out_rds) {

  # ---- Load DGEList if input is a path ----
  # If `me_dge` is a character scalar, interpret it as an .rds file path.
  if (is.character(me_dge) && length(me_dge) == 1L) {
    if (!file.exists(me_dge)) stop("me_dge path does not exist: ", me_dge)
    me_dge <- readRDS(me_dge)
  }

  # ---- Basic checks ----
  if (missing(CpGs_loc_dir) || !dir.exists(CpGs_loc_dir)) {
    stop("CpGs_loc_dir does not exist: ", CpGs_loc_dir)
  }
  if (missing(summary_out_csv) || !nzchar(summary_out_csv)) {
    stop("summary_out_csv is required.")
  }
  if (missing(dge_out_rds) || !nzchar(dge_out_rds)) {
    stop("dge_out_rds is required.")
  }
  if (is.null(me_dge$genes)) {
    stop("me_dge$genes is NULL; expected genes data.frame with Chr and Locus.")
  }
  if (!all(c("Chr", "Locus") %in% colnames(me_dge$genes))) {
    stop("me_dge$genes must contain columns: Chr and Locus.")
  }

  # Work on a local copy of genes, normalize types.
  genes <- me_dge$genes
  genes$Chr <- as.character(genes$Chr)
  genes$Locus <- suppressWarnings(as.integer(genes$Locus))

  # `keep_all` marks which original rows to keep in the final filtered DGEList.
  keep_all <- rep(FALSE, nrow(genes))

  # Strand vector aligned to the original genes rows; filled only for kept rows.
  strand_all <- rep(NA_character_, nrow(genes))

  # Per-chromosome summary rows will be accumulated here.
  summary_list <- vector("list", length(ChrNames))
  names(summary_list) <- ChrNames

  # ---- Iterate across chromosomes ----
  for (chr in ChrNames) {

    # Expected CSV path for this chromosome
    csv_path <- file.path(CpGs_loc_dir, paste0(file_prefix, chr, ".csv"))
    if (!file.exists(csv_path)) {
      message("Missing CSV: ", csv_path, " (skipping)")
      next
    }

    # Read CpG positions and validate schema
    cpg_df <- read.csv(csv_path, stringsAsFactors = FALSE)
    if (!(position_col %in% names(cpg_df))) {
      stop("CSV ", csv_path, " does not contain column: ", position_col)
    }

    # Unique integer positions p (NA removed)
    pos <- unique(suppressWarnings(as.integer(cpg_df[[position_col]])))
    pos <- pos[!is.na(pos)]
    pos_plus1 <- pos + 1L

    # Indices of DGE rows for this chromosome
    idx_chr <- which(genes$Chr == chr)
    locus_chr <- genes$Locus[idx_chr]

    # If no rows in DGE for this chromosome, still record a summary row.
    if (length(idx_chr) == 0) {
      summary_list[[chr]] <- data.frame(
        chromosome = chr,
        n_csv_positions = length(pos),
        n_csv_positions_mapped = 0L,
        pct_csv_positions_mapped = if (length(pos) == 0) NA_real_ else 0,
        n_dge_rows_chr = 0L,
        n_dge_rows_mapped = 0L,
        n_dge_rows_unmapped = 0L,
        pct_dge_rows_mapped = NA_real_,
        n_plus = 0L, n_minus = 0L, n_ambig = 0L,
        pct_plus_among_mapped = NA_real_,
        pct_minus_among_mapped = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Row-level matches
    hit_plus  <- locus_chr %in% pos       # locus == p
    hit_minus <- locus_chr %in% pos_plus1 # locus == p+1
    hit_any   <- hit_plus | hit_minus     # keep if either

    # Strand assignment for each DGE row on this chromosome
    strand_chr <- rep(NA_character_, length(idx_chr))
    strand_chr[ hit_plus & !hit_minus] <- "+"
    strand_chr[!hit_plus &  hit_minus] <- "-"
    strand_chr[ hit_plus &  hit_minus] <- "ambig"

    # Update global keep and strand for rows on this chromosome
    keep_all[idx_chr] <- hit_any
    strand_all[idx_chr[hit_any]] <- strand_chr[hit_any]

    # Compute how many CSV positions are "mapped" (p or p+1 appears in DGE loci)
    locus_set <- unique(locus_chr)
    n_csv_positions <- length(pos)
    n_csv_positions_mapped <- sum((pos %in% locus_set) | ((pos + 1L) %in% locus_set))
    pct_csv_positions_mapped <- if (n_csv_positions == 0) NA_real_ else 100 * n_csv_positions_mapped / n_csv_positions

    # DGE mapping statistics
    n_dge_rows_chr <- length(idx_chr)
    n_dge_rows_mapped <- sum(hit_any)
    n_dge_rows_unmapped <- n_dge_rows_chr - n_dge_rows_mapped
    pct_dge_rows_mapped <- if (n_dge_rows_chr == 0) NA_real_ else 100 * n_dge_rows_mapped / n_dge_rows_chr

    # Strand breakdown among mapped rows
    n_plus  <- sum(strand_chr == "+", na.rm = TRUE)
    n_minus <- sum(strand_chr == "-", na.rm = TRUE)
    n_ambig <- sum(strand_chr == "ambig", na.rm = TRUE)

    denom_mapped <- n_plus + n_minus + n_ambig
    pct_plus_among_mapped  <- if (denom_mapped == 0) NA_real_ else 100 * n_plus  / denom_mapped
    pct_minus_among_mapped <- if (denom_mapped == 0) NA_real_ else 100 * n_minus / denom_mapped

    summary_list[[chr]] <- data.frame(
      chromosome = chr,
      n_csv_positions = n_csv_positions,
      n_csv_positions_mapped = n_csv_positions_mapped,
      pct_csv_positions_mapped = pct_csv_positions_mapped,
      n_dge_rows_chr = n_dge_rows_chr,
      n_dge_rows_mapped = n_dge_rows_mapped,
      n_dge_rows_unmapped = n_dge_rows_unmapped,
      pct_dge_rows_mapped = pct_dge_rows_mapped,
      n_plus = n_plus, n_minus = n_minus, n_ambig = n_ambig,
      pct_plus_among_mapped = pct_plus_among_mapped,
      pct_minus_among_mapped = pct_minus_among_mapped,
      stringsAsFactors = FALSE
    )
  }

  # Combine summaries (dropping skipped chromosomes without entries)
  summary_df <- do.call(rbind, summary_list)
  summary_df <- summary_df[!is.na(summary_df$chromosome), , drop = FALSE]

  # ---- Create filtered DGEList ----
  # keep.lib.sizes = FALSE ensures library sizes are recalculated.
  me_dge_cpg <- me_dge[keep_all, , keep.lib.sizes = FALSE]

  # Add strand annotation to genes in the filtered object.
  me_dge_cpg$genes$Strand <- strand_all[keep_all]

  # ---- Write outputs (always) ----
  dir.create(dirname(summary_out_csv), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(dge_out_rds), recursive = TRUE, showWarnings = FALSE)

  write.csv(summary_df, summary_out_csv, row.names = FALSE)
  saveRDS(me_dge_cpg, dge_out_rds)

  # Return the filtered DGEList
  return(me_dge_cpg)
}
