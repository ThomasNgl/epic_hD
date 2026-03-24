#!/usr/bin/env bash
# =============================================================================
# 02_build_bed12.sh
#
# Purpose : Convert the GENCODE M38 GTF to BED12 format.
#           Saved alongside the GTF in the annotation folder.
#
# Output  : /mnt/auxiliary/gencode_primary/annotation/gencode.vM38.mm39.bed12
#
#   BED12 encodes the full exon structure of every transcript in 12 columns:
#     col 1   chr
#     col 2   txStart    (0-based)
#     col 3   txEnd      (1-based exclusive)
#     col 4   name       (transcript ID, e.g. ENSMUST00000193812.2)
#     col 5   score      (0)
#     col 6   strand
#     col 7   thickStart (CDS start; same as txStart for non-coding)
#     col 8   thickEnd   (CDS end;   same as txEnd   for non-coding)
#     col 9   itemRgb    (0)
#     col 10  exonCount
#     col 11  exonSizes  (comma-separated)
#     col 12  exonStarts (comma-separated, relative to txStart)
#
#   Tools that require BED12:
#     RSeQC      : infer_experiment.py, read_distribution.py,
#                  junction_saturation.py — classify reads as exonic/intronic
#     deeptools  : computeMatrix with --metagene
#     IGV / UCSC : renders transcript structure directly from BED12
#
#   Conversion pipeline: GTF → genePred → BED12
#     gtfToGenePred and genePredToBed are UCSC utility tools available
#     via bioconda (ucsc-gtftogenepred, ucsc-genepredtobed).
#
# Usage   : bash 02_build_bed12.sh
#
# Dependencies : gtfToGenePred, genePredToBed, sort
#
# Input   : /mnt/auxiliary/gencode_primary/annotation/
#               gencode.vM38.primary_assembly.annotation.gtf
#           (produced by 00_download_gencode_m38.sh)
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────

ANNOT_DIR="/mnt/auxiliary/gencode_primary/annotation"
GTF="${ANNOT_DIR}/gencode.vM38.primary_assembly.annotation.gtf"
OUT="${ANNOT_DIR}/gencode.vM38.mm39.bed12"
LOG="${ANNOT_DIR}/02_build_bed12.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "======================================================"
log "02_build_bed12.sh — epic_hD Aim 1"
log "GTF : ${GTF}"
log "Out : ${OUT}"
log "======================================================"

# ─────────────────────────────────────────────
# 1. DEPENDENCY AND INPUT CHECK
# ─────────────────────────────────────────────

for bin in gtfToGenePred genePredToBed sort; do
    command -v "${bin}" &>/dev/null || {
        log "ERROR: '${bin}' not found."
        log "       Install via: conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed"
        exit 1
    }
done

[[ -f "${GTF}" ]] || {
    log "ERROR: GTF not found: ${GTF}"
    log "       Run 00_download_gencode_m38.sh first."
    exit 1
}

# ─────────────────────────────────────────────
# 2. BUILD BED12
# ─────────────────────────────────────────────

if [[ -f "${OUT}" ]]; then
    log "Output already exists, skipping: ${OUT}"
    exit 0
fi

log "Converting GTF → genePred → BED12..."
log "This typically takes 1–2 minutes."

# -genePredExt : write extended genePred with extra annotation fields,
#                required by some RSeQC modules
# -allErrors   : report all conversion warnings instead of stopping at first;
#                GENCODE GTFs can produce harmless warnings for non-standard
#                features (e.g. readthrough transcripts) that should not abort
#                the run
# Pipe directly through genePredToBed to avoid a large temporary file,
# then sort by chromosome and start position (required by RSeQC and deeptools)
gtfToGenePred \
    -genePredExt \
    -allErrors \
    "${GTF}" stdout \
    2>"${LOG}.warnings" \
| genePredToBed stdin stdout \
| sort -k1,1 -k2,2n \
> "${OUT}"

N_TX=$(wc -l < "${OUT}")
log "Done. ${N_TX} transcripts written to ${OUT}"

# Report conversion warnings if any
N_WARN=$(wc -l < "${LOG}.warnings")
if [[ "${N_WARN}" -gt 0 ]]; then
    log "  ${N_WARN} conversion warning(s) — see ${LOG}.warnings"
else
    log "  No conversion warnings."
    rm -f "${LOG}.warnings"
fi

log "======================================================"
log "Full log : ${LOG}"
log "======================================================"
