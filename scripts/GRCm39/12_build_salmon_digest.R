# =============================================================================
# 12_build_salmon_digest.R
#
# Purpose : Register the GENCODE M38 Salmon index with tximeta by creating
#           a linked transcriptome (linkedTxome). This is required so that
#           tximeta::tximeta() can automatically attach transcript and gene
#           metadata when importing Salmon quantification output.
#
#           Run this script ONCE after building the Salmon index with
#           11_build_salmon_index.sh. The result is a JSON file stored in
#           the BiocFileCache that links the index to its source GTF and FASTA.
#
# Usage   : Rscript 12_build_salmon_digest.R
#        or: source("12_build_salmon_digest.R") from an R session
#
# Output  : A JSON file written to the BiocFileCache at /mnt/auxiliary/biocache
#           After this runs, tximeta() will recognise any quant.sf file
#           produced with mm39_vM38_salmon_index without further configuration.
#
# Dependencies : tximeta (Bioconductor)
#
# Inputs  : /mnt/auxiliary/salmon/mm39_vM38_salmon_index/   (index directory)
#           /mnt/auxiliary/gencode_primary/transcriptome/gencode.vM38.transcripts.fa
#           /mnt/auxiliary/gencode_primary/annotation/gencode.vM38.primary_assembly.annotation.gtf
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

suppressPackageStartupMessages(library(tximeta))

log_msg <- function(...) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

# ─────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────

INDEX_DIR  <- "/mnt/auxiliary/salmon/mm39_vM38_salmon_index"
FASTA      <- "/mnt/auxiliary/gencode_primary/transcriptome/gencode.vM38.transcripts.fa"
GTF        <- "/mnt/auxiliary/gencode_primary/annotation/gencode.vM38.primary_assembly.annotation.gtf"
BIOC_CACHE <- "/mnt/auxiliary/biocache"

# ─────────────────────────────────────────────
# 1. INPUT CHECK
# ─────────────────────────────────────────────

for (path in c(INDEX_DIR, FASTA, GTF)) {
    if (!file.exists(path)) {
        stop(sprintf("Expected input not found: %s\n  Run the upstream scripts first.", path))
    }
}

# ─────────────────────────────────────────────
# 2. SET BiocFileCache DIRECTORY
#
# tximeta stores TxDb objects and linkedTxome JSON files here.
# Must be set before calling makeLinkedTxome.
# ─────────────────────────────────────────────

log_msg("Setting BiocFileCache directory: ", BIOC_CACHE)
dir.create(BIOC_CACHE, showWarnings = FALSE, recursive = TRUE)
setTximetaBFC(BIOC_CACHE)

# ─────────────────────────────────────────────
# 3. REGISTER THE LINKED TRANSCRIPTOME
#
# makeLinkedTxome writes a JSON file that associates:
#   - the Salmon index (identified by its hash)
#   - the source FASTA and GTF used to build it
#   - release metadata (source, organism, release, genome)
#
# After this, tximeta() called on any quant.sf from this index will
# automatically load the correct transcript annotations without
# manual specification of the GTF.
# ─────────────────────────────────────────────

log_msg("Registering linkedTxome for GENCODE M38...")
log_msg("  Index  : ", INDEX_DIR)
log_msg("  FASTA  : ", FASTA)
log_msg("  GTF    : ", GTF)

makeLinkedTxome(
    indexDir  = INDEX_DIR,
    source    = "GENCODE",
    organism  = "Mus musculus",
    release   = "M38",
    genome    = "GRCm39",
    fasta     = FASTA,
    gtf       = GTF,
    write     = TRUE    # write JSON to BiocFileCache for permanent registration
)

log_msg("Done. linkedTxome registered in: ", BIOC_CACHE)
log_msg("tximeta() will now automatically recognise quant.sf files")
log_msg("produced with: ", INDEX_DIR)
