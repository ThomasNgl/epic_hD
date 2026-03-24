#!/usr/bin/env bash
# =============================================================================
# 10_prepare_salmon_inputs.sh
#
# Purpose : Build the two files Salmon needs before index construction:
#
#   decoys.txt  — one chromosome/scaffold name per line, extracted from the
#                 genome FASTA headers. Tells Salmon which sequences are
#                 genomic decoys (not transcribed) so that reads mapping to
#                 pseudogenes or unannotated genomic copies of genes are
#                 discarded rather than counted as transcript alignments.
#
#   gentrome.fa — the transcript FASTA concatenated with the genome FASTA
#                 (transcripts first, genome second). Salmon requires a single
#                 input file when building a decoy-aware index. The decoys.txt
#                 tells it where the transcriptome ends and the genome begins.
#
# Output  : /mnt/auxiliary/gencode_primary/salmon_prep/
#               decoys.txt
#               gentrome.fa
#
# Usage   : bash 10_prepare_salmon_inputs.sh
#
# Dependencies : grep, cut, sed, cat
#
# Inputs  : /mnt/auxiliary/gencode_primary/
#               genome/GRCm39.primary_assembly.genome.fa
#               transcriptome/gencode.vM38.transcripts.fa
#           (produced by 00_download_gencode_m38.sh)
#
# Notes   : gentrome.fa is ~2.8 Gb. It is only needed to build the Salmon
#           index and can be deleted afterwards if disk space is tight.
#           decoys.txt is tiny (~100 lines) and should be kept permanently.
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────

BASE_DIR="/mnt/auxiliary/gencode_primary"
GENOME_FA="${BASE_DIR}/genome/GRCm39.primary_assembly.genome.fa"
TX_FA="${BASE_DIR}/transcriptome/gencode.vM38.transcripts.fa"

SALMON_PREP_DIR="${BASE_DIR}/salmon_prep"
DECOYS="${SALMON_PREP_DIR}/decoys.txt"
GENTROME="${SALMON_PREP_DIR}/gentrome.fa"

LOG="${SALMON_PREP_DIR}/10_prepare_salmon_inputs.log"

mkdir -p "${SALMON_PREP_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "======================================================"
log "10_prepare_salmon_inputs.sh — epic_hD Aim 1"
log "Output dir : ${SALMON_PREP_DIR}"
log "======================================================"

# ─────────────────────────────────────────────
# 1. INPUT CHECK
# ─────────────────────────────────────────────

for f in "${GENOME_FA}" "${TX_FA}"; do
    [[ -f "${f}" ]] || {
        log "ERROR: Expected input not found: ${f}"
        log "       Run 00_download_gencode_m38.sh first."
        exit 1
    }
done

# ─────────────────────────────────────────────
# 2. DECOYS
#
# Extract sequence names from genome FASTA headers.
# Each header looks like: >chr1 1 dna:chromosome chromosome:GRCm39:1:1:195154279:1 REF
# We take only the first token after '>' to get the clean name (chr1, chrX, …)
# and strip the leading '>'.
# ─────────────────────────────────────────────

if [[ -f "${DECOYS}" ]]; then
    log "decoys.txt already exists, skipping."
else
    log "Building decoys.txt..."
    grep '^>' "${GENOME_FA}" \
        | cut -d ' ' -f1 \
        | sed 's/^>//' \
        > "${DECOYS}"
    log "decoys.txt: $(wc -l < "${DECOYS}") sequences"
fi

# ─────────────────────────────────────────────
# 3. GENTROME
#
# Concatenate transcripts then genome into one file.
# Order is mandatory: transcripts must come first.
# Salmon uses decoys.txt to identify where the genomic
# sequences begin inside this combined file.
# ─────────────────────────────────────────────

if [[ -f "${GENTROME}" ]]; then
    log "gentrome.fa already exists, skipping."
else
    log "Building gentrome.fa (transcripts + genome)..."
    log "This will take a few minutes — the genome is ~2.7 Gb."
    cat "${TX_FA}" "${GENOME_FA}" > "${GENTROME}"
    log "gentrome.fa: $(du -sh "${GENTROME}" | cut -f1)"
fi

# ─────────────────────────────────────────────
# 4. SUMMARY
# ─────────────────────────────────────────────

log "======================================================"
log "Salmon prep complete."
log ""
log "  Decoys   : ${DECOYS}"
log "  Gentrome : ${GENTROME}"
log ""
log "  Next step: build the Salmon index with"
log "    salmon index --transcripts ${GENTROME} \\"
log "                 --decoys ${DECOYS} \\"
log "                 --index <output_dir> \\"
log "                 --gencode --threads <N>"
log ""
log "  Note: gentrome.fa can be deleted after the index is built."
log "  Full log : ${LOG}"
log "======================================================"
