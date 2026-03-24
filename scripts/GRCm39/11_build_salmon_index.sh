#!/usr/bin/env bash
# =============================================================================
# 11_build_salmon_index.sh
#
# Purpose : Build the Salmon decoy-aware index from the gentrome and decoys
#           prepared by 10_prepare_salmon_inputs.sh.
#
# Output  : /mnt/auxiliary/salmon/mm39_vM38_salmon_index/
#
# Usage   : bash 11_build_salmon_index.sh [--threads N]
#
#   --threads N   Number of threads (default: 8)
#   -h, --help    Show this help message and exit
#
# Dependencies : salmon
#
# Inputs  : /mnt/auxiliary/gencode_primary/salmon_prep/
#               gentrome.fa    (produced by 10_prepare_salmon_inputs.sh)
#               decoys.txt     (produced by 10_prepare_salmon_inputs.sh)
#
# Notes on key flags:
#   --gencode   Strips GENCODE-specific version suffixes from transcript IDs
#               (e.g. ENSMUST00000001234.5 → ENSMUST00000001234) so that IDs
#               match between the index, the metadata files, and tximeta/tximport.
#               Always use this flag when indexing GENCODE transcriptomes.
#
#   --decoys    Enables the selective alignment strategy. Reads that map
#               preferentially to the genome (pseudogenes, unannotated copies)
#               are discarded rather than mis-assigned to transcripts.
#               Required when using a gentrome as input.
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# 0. DEFAULT PARAMETERS
# ─────────────────────────────────────────────

THREADS=8

usage() {
    sed -n '2,/^set -euo pipefail$/p' "$0" | sed '$d;s/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads) THREADS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; echo "Use --help." >&2; exit 1 ;;
    esac
done

if ! [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer." >&2; exit 1
fi

# ─────────────────────────────────────────────
# 1. PATHS
# ─────────────────────────────────────────────

SALMON_PREP_DIR="/mnt/auxiliary/gencode_primary/salmon_prep"
GENTROME="${SALMON_PREP_DIR}/gentrome.fa"
DECOYS="${SALMON_PREP_DIR}/decoys.txt"

INDEX_DIR="/mnt/auxiliary/salmon/mm39_vM38_salmon_index"
LOG="${INDEX_DIR}/11_build_salmon_index.log"

mkdir -p "${INDEX_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "======================================================"
log "11_build_salmon_index.sh — epic_hD Aim 1"
log "Gentrome : ${GENTROME}"
log "Decoys   : ${DECOYS}"
log "Index    : ${INDEX_DIR}"
log "Threads  : ${THREADS}"
log "======================================================"

# ─────────────────────────────────────────────
# 2. DEPENDENCY AND INPUT CHECK
# ─────────────────────────────────────────────

command -v salmon &>/dev/null || {
    log "ERROR: salmon not found. Activate your conda environment."
    exit 1
}

for f in "${GENTROME}" "${DECOYS}"; do
    [[ -f "${f}" ]] || {
        log "ERROR: Expected input not found: ${f}"
        log "       Run 10_prepare_salmon_inputs.sh first."
        exit 1
    }
done

# ─────────────────────────────────────────────
# 3. BUILD INDEX
# ─────────────────────────────────────────────

# Check if index already exists (Salmon writes an info.json on completion)
if [[ -f "${INDEX_DIR}/info.json" ]]; then
    log "Salmon index already exists, skipping: ${INDEX_DIR}"
    exit 0
fi

log "Building Salmon index..."
log "Typical runtime: 10–20 min | RAM: ~8 GB"

salmon index \
    --threads "${THREADS}" \
    --transcripts "${GENTROME}" \
    --decoys "${DECOYS}" \
    --index "${INDEX_DIR}" \
    --gencode

log "======================================================"
log "Salmon index complete."
log ""
log "  Index : ${INDEX_DIR}"
log ""
log "  Note: gentrome.fa (~3 Gb) can now be deleted if disk space is needed:"
log "    rm ${GENTROME}"
log ""
log "  Full log : ${LOG}"
log "======================================================"
