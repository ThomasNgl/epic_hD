# =============================================================================
# salmon_to_se.R
#
# Purpose : Import Salmon quantification output into a SummarizedExperiment
#           using tximeta, attach gene-level metadata (symbol, biotype), and
#           optionally save the result as an RDS file.
#
# Function : salmon_to_se(quant_dir, sample_meta, output_path, ...)
#
# Usage example:
#
#   source("src/R/salmon_to_se.R")
#
#   meta <- read.csv("metadata.csv", row.names = 1)
#
#   se <- salmon_to_se(
#       quant_dir   = "/mnt/thomas/phd_proj/results/sg_oocyte/05_Salmon_quant",
#       sample_meta = meta,
#       output_path = "/mnt/thomas/phd_proj/results/sg_oocyte/06_Salmon_se/se_gene.rds"
#   )
#
# Dependencies : tximeta, SummarizedExperiment (Bioconductor)
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD
# =============================================================================

suppressPackageStartupMessages({
    library(tximeta)
    library(SummarizedExperiment)
})

# =============================================================================
# salmon_to_se()
#
# Parameters:
#
#   quant_dir     (character)
#     Path to the Salmon quantification directory. Expected structure:
#       quant_dir/
#         <sample_id>_quant/
#           quant.sf
#     The sample subdirectory name is constructed as paste0(sample_id, "_quant").
#     sample_id values are taken from rownames(sample_meta).
#
#   sample_meta   (data.frame)
#     Sample metadata table. rownames must be the sample IDs matching the
#     Salmon output subdirectory names (without the '_quant' suffix).
#     All columns are carried into colData of the output SE.
#
#   output_path   (character or NULL)
#     If provided, the SummarizedExperiment is saved as an RDS to this path.
#     The directory is created if it does not exist.
#     If NULL, the SE is returned without saving.
#
#   level         (character, default "gene")
#     Summarisation level. One of:
#       "gene"       — summarise transcript counts to gene level via
#                      tximeta::summarizeToGene(). rownames are versioned
#                      Ensembl gene IDs (e.g. ENSMUSG00000051951.6).
#       "transcript" — keep transcript-level counts. rownames are versioned
#                      Ensembl transcript IDs.
#
#   gene_map_path (character)
#     Path to the gene map RDS produced by 13_build_gene_map.R.
#     Used to attach gene_symbol and gene_biotype to rowData.
#     Only used when level = "gene".
#     Default: "/mnt/auxiliary/gencode_primary/annotation/gencode.vM38.gene_map.rds"
#
#   bioc_cache    (character)
#     Path to the BiocFileCache directory used by tximeta to find the
#     linkedTxome JSON registered by 12_build_salmon_digest.R.
#     Default: "/mnt/auxiliary/biocache"
#
# Returns:
#   A SummarizedExperiment object with:
#     assays   : counts, abundance (TPM), length
#     colData  : all columns from sample_meta
#     rowData  : gene_id_base, gene_symbol, gene_biotype  (if level="gene")
# =============================================================================

salmon_to_se <- function(
    quant_dir,
    sample_meta,
    output_path   = NULL,
    level         = "gene",
    gene_map_path = "/mnt/auxiliary/gencode_primary/annotation/gencode.vM38.gene_map.rds",
    bioc_cache    = "/mnt/auxiliary/biocache"
) {

    # ── Argument validation ──────────────────────────────────────────────────

    stopifnot(
        "quant_dir must be a single character string" =
            is.character(quant_dir) && length(quant_dir) == 1,
        "sample_meta must be a data.frame" =
            is.data.frame(sample_meta),
        "sample_meta must have rownames (sample IDs)" =
            !is.null(rownames(sample_meta)),
        "level must be 'gene' or 'transcript'" =
            level %in% c("gene", "transcript")
    )

    # ── Locate quant.sf files ────────────────────────────────────────────────

    sample_ids <- rownames(sample_meta)
    quant_files <- file.path(
        quant_dir,
        paste0(sample_ids, "_quant"),
        "quant.sf"
    )
    names(quant_files) <- sample_ids

    missing <- quant_files[!file.exists(quant_files)]
    if (length(missing) > 0) {
        warning(
            length(missing), " quant.sf file(s) not found and will be excluded:\n",
            paste("  ", missing, collapse = "\n")
        )
        keep        <- file.exists(quant_files)
        quant_files <- quant_files[keep]
        sample_meta <- sample_meta[keep, , drop = FALSE]
    }

    if (length(quant_files) == 0) {
        stop("No quant.sf files found under: ", quant_dir)
    }

    message(sprintf("[salmon_to_se] Importing %d sample(s)", length(quant_files)))

    # ── Build tximeta coldata ────────────────────────────────────────────────
    # tximeta requires 'files' and 'names' as the first two columns

    coldata        <- sample_meta
    coldata$files  <- quant_files
    coldata$names  <- sample_ids
    coldata        <- coldata[, c("files", "names",
                                  setdiff(colnames(coldata), c("files", "names"))),
                               drop = FALSE]

    # ── Set BiocFileCache so tximeta finds the linkedTxome ───────────────────

    setTximetaBFC(bioc_cache)

    # ── Import at transcript level ───────────────────────────────────────────
    # useHub = FALSE: use local linkedTxome registered by 12_build_salmon_digest.R
    # instead of downloading from AnnotationHub

    se_tx <- tximeta(coldata, useHub = FALSE)
    message("[salmon_to_se] Transcript-level import complete")

    # ── Summarise to gene level if requested ─────────────────────────────────

    if (level == "gene") {
        se <- summarizeToGene(se_tx)
        message(sprintf(
            "[salmon_to_se] Summarised to %d genes across %d samples",
            nrow(se), ncol(se)
        ))
        rm(se_tx)

        # Attach gene_symbol and gene_biotype from pre-built gene map
        if (!file.exists(gene_map_path)) {
            warning(
                "gene_map_path not found: ", gene_map_path,
                "\nSkipping symbol/biotype annotation. ",
                "Run 13_build_gene_map.R first."
            )
        } else {
            gene_map <- readRDS(gene_map_path)

            # Match on versioned gene ID (rownames of SE match gene_map$gene_id)
            idx <- match(rownames(se), gene_map$gene_id)

            n_unmatched <- sum(is.na(idx))
            if (n_unmatched > 0) {
                message(sprintf(
                    "[salmon_to_se] %d gene(s) not matched in gene map (NA)",
                    n_unmatched
                ))
            }

            rowData(se)$gene_id_base  <- gene_map$gene_id_base[idx]
            rowData(se)$gene_symbol   <- gene_map$gene_symbol[idx]
            rowData(se)$gene_biotype  <- gene_map$gene_biotype[idx]

            message("[salmon_to_se] gene_symbol and gene_biotype attached to rowData")
        }

    } else {
        # transcript level — return as-is
        se <- se_tx
        message(sprintf(
            "[salmon_to_se] Returning transcript-level SE: %d transcripts x %d samples",
            nrow(se), ncol(se)
        ))
    }

    # ── Verify colData integrity ─────────────────────────────────────────────

    if (!identical(colnames(se), sample_ids[file.exists(quant_files)])) {
        warning("[salmon_to_se] Sample order in SE may not match sample_meta row order.")
    }

    # ── Save if output_path provided ─────────────────────────────────────────

    if (!is.null(output_path)) {
        dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(se, file = output_path)
        message("[salmon_to_se] SE saved to: ", output_path)
    }

    invisible(se)
}
