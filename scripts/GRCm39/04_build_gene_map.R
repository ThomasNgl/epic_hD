# =============================================================================
# 04_build_gene_map.R
#
# Purpose : Parse the GENCODE M38 GTF once and save a gene-level metadata
#           table as an RDS file. This avoids re-parsing the large GTF on
#           every analysis run.
#
#           The gene map is used by salmon_to_se() to attach gene symbols
#           and biotypes to the SummarizedExperiment rowData after Salmon
#           quantification import.
#
# Output  : /mnt/auxiliary/gencode_primary/annotation/gencode.vM38.gene_map.rds
#
#   A data.frame with columns:
#     gene_id      — Ensembl ID with version (e.g. ENSMUSG00000051951.6)
#                    kept with version to match tximeta output exactly
#     gene_id_base — version suffix stripped (ENSMUSG00000051951)
#     gene_symbol  — MGI approved name (gene_name in GTF)
#     gene_biotype — transcript class (gene_type in GTF)
#
# Usage   : Rscript 04_build_gene_map.R
#
# Dependencies : rtracklayer (Bioconductor)
#
# Input   : /mnt/auxiliary/gencode_primary/annotation/
#               gencode.vM38.primary_assembly.annotation.gtf
#           (produced by 00_download_gencode_m38.sh)
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

suppressPackageStartupMessages(library(rtracklayer))

log_msg <- function(...) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

# ─────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────

GTF     <- "/mnt/auxiliary/gencode_primary/annotation/gencode.vM38.primary_assembly.annotation.gtf"
OUT_RDS <- "/mnt/auxiliary/gencode_primary/annotation/gencode.vM38.gene_map.rds"

# ─────────────────────────────────────────────
# 1. INPUT CHECK
# ─────────────────────────────────────────────

if (!file.exists(GTF)) {
    stop("GTF not found: ", GTF, "\nRun 00_download_gencode_m38.sh first.")
}

if (file.exists(OUT_RDS)) {
    log_msg("Gene map already exists, skipping: ", OUT_RDS)
    quit(status = 0)
}

# ─────────────────────────────────────────────
# 2. PARSE GTF — GENE FEATURES ONLY
#
# rtracklayer::import() reads the full GTF into a GRanges object.
# We keep only 'gene' features (one entry per gene) to avoid
# duplicates that would arise from transcript/exon features.
#
# We retain the versioned gene_id (e.g. ENSMUSG00000051951.6) because
# tximeta/summarizeToGene returns versioned IDs in rownames by default,
# and match() requires the IDs to be identical.
# The stripped version (gene_id_base) is provided for convenience when
# joining to the annotation TSV (gencode.vM38.genes.tsv) which uses
# unversioned IDs.
# ─────────────────────────────────────────────

log_msg("Parsing GTF: ", GTF)
log_msg("This takes 30–60 seconds...")

gtf_gr  <- import(GTF)
gene_gr <- gtf_gr[gtf_gr$type == "gene"]
rm(gtf_gr)
gc()

log_msg("Gene features found: ", length(gene_gr))

gene_map <- data.frame(
    gene_id      = gene_gr$gene_id,
    gene_id_base = sub("\\.[0-9]+$", "", gene_gr$gene_id),
    gene_symbol  = ifelse(is.null(gene_gr$gene_name) | is.na(gene_gr$gene_name),
                          NA_character_, gene_gr$gene_name),
    gene_biotype = ifelse(is.null(gene_gr$gene_type) | is.na(gene_gr$gene_type),
                          NA_character_, gene_gr$gene_type),
    stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────
# 3. SANITY CHECK AND SAVE
# ─────────────────────────────────────────────

n_na_sym <- sum(is.na(gene_map$gene_symbol))
n_na_bio <- sum(is.na(gene_map$gene_biotype))

log_msg("Genes in map     : ", nrow(gene_map))
if (n_na_sym > 0) log_msg("WARNING: ", n_na_sym, " gene(s) with no gene_symbol (NA)")
if (n_na_bio > 0) log_msg("WARNING: ", n_na_bio, " gene(s) with no gene_biotype (NA)")

log_msg("Saving gene map to: ", OUT_RDS)
saveRDS(gene_map, file = OUT_RDS)
log_msg("Done.")
