#!/usr/bin/env bash
# =============================================================================
# 20_build_star_index.sh
#
# Purpose : Build the STAR genome index for splice-aware alignment of
#           long RNA-seq reads (mRNA, lncRNA) from sperm, oocyte, and
#           zygote public datasets.
#
# Output  : /mnt/auxiliary/star/mm39_vM38_star_index_oh{OVERHANG}/
#
#   The output directory name encodes the overhang value so that indices
#   built for different read lengths coexist without conflict.
#
# Usage   : bash 20_build_star_index.sh [--threads N] [--overhang N]
#
#   --threads N    Number of threads (default: 8)
#   --overhang N   sjdbOverhang = read_length - 1 (default: 100)
#                  Set to read_length - 1 for a specific dataset.
#                  100 is a safe general-purpose value for reads 75–150 bp.
#                  Examples: 74 for 75 bp reads, 124 for 125 bp, 149 for 150 bp.
#   -h, --help     Show this help message and exit
#
# Dependencies : STAR
#
# Inputs  : /mnt/auxiliary/gencode_primary/
#               genome/GRCm39.primary_assembly.genome.fa
#               annotation/gencode.vM38.primary_assembly.annotation.gtf
#           (produced by 00_download_gencode_m38.sh)
#
# Notes on key parameters:
#   --sjdbOverhang         Controls the length of the junction donor/acceptor
#                          sequence stored in the index. Optimal = read_length - 1.
#                          100 works well across 75–150 bp; only rebuild if all
#                          your datasets have a single known read length.
#
#   --genomeSAindexNbases  14 is standard for mammalian genomes (~2.7 Gb).
#                          Formula: min(14, log2(GenomeLength)/2 - 1).
#
#   RAM requirement        ~32 GB during index build. Run as a cluster job
#                          or on a node with sufficient memory.
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# 0. DEFAULT PARAMETERS
# ─────────────────────────────────────────────

THREADS=8
OVERHANG=100

usage() {
    sed -n '2,/^set -euo pipefail$/p' "$0" | sed '$d;s/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)  THREADS="$2";  shift 2 ;;
        --overhang) OVERHANG="$2"; shift 2 ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; echo "Use --help." >&2; exit 1 ;;
    esac
done

if ! [[ "${THREADS}"  =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer." >&2; exit 1
fi
if ! [[ "${OVERHANG}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --overhang must be a positive integer." >&2; exit 1
fi

# ─────────────────────────────────────────────
# 1. PATHS
# ─────────────────────────────────────────────

GENOME_FA="/mnt/auxiliary/gencode_primary/genome/GRCm39.primary_assembly.genome.fa"
ANNOT_GTF="/mnt/auxiliary/gencode_primary/annotation/gencode.vM38.primary_assembly.annotation.gtf"

INDEX_DIR="/mnt/auxiliary/star/mm39_vM38_star_index_oh${OVERHANG}"
LOG="${INDEX_DIR}/20_build_star_index.log"

mkdir -p "${INDEX_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "======================================================"
log "20_build_star_index.sh — epic_hD Aim 1"
log "Genome   : ${GENOME_FA}"
log "GTF      : ${ANNOT_GTF}"
log "Index    : ${INDEX_DIR}"
log "Overhang : ${OVERHANG}  (read length assumed: $((OVERHANG + 1)) bp)"
log "Threads  : ${THREADS}"
log "======================================================"

# ─────────────────────────────────────────────
# 2. DEPENDENCY AND INPUT CHECK
# ─────────────────────────────────────────────

command -v STAR &>/dev/null || {
    log "ERROR: STAR not found. Activate your conda environment."
    exit 1
}

for f in "${GENOME_FA}" "${ANNOT_GTF}"; do
    [[ -f "${f}" ]] || {
        log "ERROR: Expected input not found: ${f}"
        log "       Run 00_download_gencode_m38.sh first."
        exit 1
    }
done

# ─────────────────────────────────────────────
# 3. BUILD INDEX
# ─────────────────────────────────────────────

# STAR writes genomeParameters.txt on successful completion
if [[ -f "${INDEX_DIR}/genomeParameters.txt" ]]; then
    log "STAR index already exists, skipping: ${INDEX_DIR}"
    exit 0
fi

log "Building STAR index..."
log "Typical runtime: 30–60 min | RAM: ~32 GB"

STAR \
    --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${INDEX_DIR}" \
    --genomeFastaFiles "${GENOME_FA}" \
    --sjdbGTFfile "${ANNOT_GTF}" \
    --sjdbOverhang "${OVERHANG}" \
    --genomeSAindexNbases 14 \
    --outTmpDir "${INDEX_DIR}/_tmp_star_$(date +%s)"

log "======================================================"
log "STAR index complete."
log ""
log "  Index    : ${INDEX_DIR}"
log "  Overhang : ${OVERHANG}"
log ""
log "  To align reads with this index:"
log "    STAR --runThreadN ${THREADS} \\"
log "         --genomeDir ${INDEX_DIR} \\"
log "         --readFilesIn <R1.fastq.gz> [<R2.fastq.gz>] \\"
log "         --readFilesCommand zcat \\"
log "         --outSAMtype BAM SortedByCoordinate \\"
log "         --outFileNamePrefix <prefix>"
log ""
log "  Full log : ${LOG}"
log "======================================================"
