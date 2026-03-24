#!/usr/bin/env bash
# =============================================================================
# 00_download_gencode_m38.sh
#
# Purpose : Download, verify (MD5), decompress, and prepare the GENCODE M38
#           mouse reference files used by this workflow.
#
# Output directory layout under /mnt/auxiliary/gencode_primary/ :
#
#   gencode_primary/
#   ├── genome/
#   │   ├── GRCm39.primary_assembly.genome.fa.gz   (raw download, kept)
#   │   ├── GRCm39.primary_assembly.genome.fa      (decompressed)
#   │   ├── GRCm39.primary_assembly.genome.fa.fai  (samtools faidx)
#   │   └── GRCm39.primary_assembly.chrom.sizes    (derived from .fai)
#   │
#   ├── annotation/
#   │   ├── gencode.vM38.primary_assembly.annotation.gtf.gz
#   │   ├── gencode.vM38.primary_assembly.annotation.gtf
#   │   ├── gencode.vM38.tRNAs.gtf.gz
#   │   └── gencode.vM38.tRNAs.gtf
#   │
#   ├── transcriptome/
#   │   ├── gencode.vM38.transcripts.fa.gz
#   │   └── gencode.vM38.transcripts.fa
#   │
#   └── metadata/
#       ├── gencode.vM38.metadata.MGI.gz
#       ├── gencode.vM38.metadata.EntrezGene.gz
#       └── gencode.vM38.metadata.RefSeq.gz
#
# Usage   : bash 00_download_gencode_m38.sh [OPTIONS]
#
#   --threads N      Number of threads for pigz decompression (default: 8)
#   -h, --help       Show this help message and exit
#
# Dependencies (must be in PATH / active conda environment):
#   wget, md5sum, samtools
# Optional:
#   pigz             Used for parallel decompression if available
#
# Notes:
#   - This script does NOT build STAR, Salmon, or bowtie1 indexes.
#   - It downloads the comprehensive primary-assembly GTF, not the "basic"
#     annotation subset.
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# 0. DEFAULT PARAMETERS
# ─────────────────────────────────────────────

BASE_DIR="/mnt/auxiliary/gencode_primary"
THREADS=8

usage() {
    sed -n '1,/^set -euo pipefail$/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown argument: $1" >&2
            echo "Use --help to see valid options." >&2
            exit 1
            ;;
    esac
done

if ! [[ "${THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer." >&2
    exit 1
fi

# ─────────────────────────────────────────────
# 1. DIRECTORY STRUCTURE
# ─────────────────────────────────────────────

GENOME_DIR="${BASE_DIR}/genome"
ANNOT_DIR="${BASE_DIR}/annotation"
TX_DIR="${BASE_DIR}/transcriptome"
META_DIR="${BASE_DIR}/metadata"

mkdir -p "${GENOME_DIR}" "${ANNOT_DIR}" "${TX_DIR}" "${META_DIR}"

# ─────────────────────────────────────────────
# 2. GENCODE M38 FTP BASE URL
# ─────────────────────────────────────────────

FTP="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38"

# ─────────────────────────────────────────────
# 3. LOGGING HELPER
# ─────────────────────────────────────────────

LOG="${BASE_DIR}/00_download_gencode_m38.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

authoritative_md5_file="${BASE_DIR}/MD5SUMS"

log "======================================================"
log "GENCODE M38 reference setup — epic_hD Aim 1"
log "Base directory : ${BASE_DIR}"
log "Threads        : ${THREADS}"
log "======================================================"

# ─────────────────────────────────────────────
# 4. DEPENDENCY CHECK
# ─────────────────────────────────────────────

log "Checking dependencies..."
REQUIRED=(wget md5sum samtools awk grep cut)

for bin in "${REQUIRED[@]}"; do
    if ! command -v "${bin}" &>/dev/null; then
        log "ERROR: '${bin}' not found in PATH. Activate your conda environment first."
        exit 1
    fi
done
log "All required tools found."

# ─────────────────────────────────────────────
# 5. HELPER FUNCTIONS
# ─────────────────────────────────────────────

get_md5_from_manifest() {
    local filename="$1"
    local md5

    md5=$(awk -v f="$filename" '$2 == f {print $1; exit}' "${authoritative_md5_file}")

    if [[ -z "${md5}" ]]; then
        log "ERROR: Could not find MD5 entry for ${filename} in ${authoritative_md5_file}"
        exit 1
    fi

    printf '%s\n' "${md5}"
}

download_and_verify() {
    local url="$1"
    local dest="$2"
    local expected_md5="$3"

    if [[ -f "${dest}" ]]; then
        log "Already present, skipping download: $(basename "${dest}")"
    else
        log "Downloading: $(basename "${dest}")"
        wget --quiet --show-progress -O "${dest}" "${url}"
    fi

    log "Verifying MD5: $(basename "${dest}")"
    local actual_md5
    actual_md5=$(md5sum "${dest}" | awk '{print $1}')

    if [[ "${actual_md5}" != "${expected_md5}" ]]; then
        log "ERROR: MD5 mismatch for $(basename "${dest}")"
        log "  Expected : ${expected_md5}"
        log "  Actual   : ${actual_md5}"
        log "  Check    : ${FTP}/MD5SUMS"
        exit 1
    fi

    log "MD5 OK: $(basename "${dest}")"
}

decompress_if_needed() {
    local src="$1"
    local dest="$2"

    if [[ -f "${dest}" ]]; then
        log "Already decompressed, skipping: $(basename "${dest}")"
    else
        log "Decompressing: $(basename "${src}")"
        if command -v pigz &>/dev/null; then
            pigz -d -k -c -p "${THREADS}" "${src}" > "${dest}"
        else
            gunzip -c "${src}" > "${dest}"
        fi
        log "Done: $(basename "${dest}")"
    fi
}

# ─────────────────────────────────────────────
# 6. DOWNLOAD PHASE
# ─────────────────────────────────────────────

log "--- Phase 1: Download ---"

wget --quiet -O "${authoritative_md5_file}" "${FTP}/MD5SUMS"
log "MD5SUMS file saved to ${authoritative_md5_file}"

# Genome FASTA: primary assembly (chromosomes + scaffolds; excludes alt loci / patches)
GENOME_GZ="${GENOME_DIR}/GRCm39.primary_assembly.genome.fa.gz"
GENOME_MD5=$(get_md5_from_manifest "GRCm39.primary_assembly.genome.fa.gz")
download_and_verify \
    "${FTP}/GRCm39.primary_assembly.genome.fa.gz" \
    "${GENOME_GZ}" \
    "${GENOME_MD5}"

# Comprehensive annotation on primary assembly regions.
# This is not the "basic" annotation subset.
ANNOT_GZ="${ANNOT_DIR}/gencode.vM38.primary_assembly.annotation.gtf.gz"
ANNOT_MD5=$(get_md5_from_manifest "gencode.vM38.primary_assembly.annotation.gtf.gz")
download_and_verify \
    "${FTP}/gencode.vM38.primary_assembly.annotation.gtf.gz" \
    "${ANNOT_GZ}" \
    "${ANNOT_MD5}"

# tRNA annotation generated by Ensembl using tRNAscan-SE.
# Distributed separately from the main GTF.
TRNA_GZ="${ANNOT_DIR}/gencode.vM38.tRNAs.gtf.gz"
TRNA_MD5=$(get_md5_from_manifest "gencode.vM38.tRNAs.gtf.gz")
download_and_verify \
    "${FTP}/gencode.vM38.tRNAs.gtf.gz" \
    "${TRNA_GZ}" \
    "${TRNA_MD5}"

# Transcript FASTA: all transcript sequences distributed for release M38.
TX_GZ="${TX_DIR}/gencode.vM38.transcripts.fa.gz"
TX_MD5=$(get_md5_from_manifest "gencode.vM38.transcripts.fa.gz")
download_and_verify \
    "${FTP}/gencode.vM38.transcripts.fa.gz" \
    "${TX_GZ}" \
    "${TX_MD5}"

# Metadata tables
META_MGI="${META_DIR}/gencode.vM38.metadata.MGI.gz"
MGI_MD5=$(get_md5_from_manifest "gencode.vM38.metadata.MGI.gz")
download_and_verify \
    "${FTP}/gencode.vM38.metadata.MGI.gz" \
    "${META_MGI}" \
    "${MGI_MD5}"

META_ENTREZ="${META_DIR}/gencode.vM38.metadata.EntrezGene.gz"
ENTREZ_MD5=$(get_md5_from_manifest "gencode.vM38.metadata.EntrezGene.gz")
download_and_verify \
    "${FTP}/gencode.vM38.metadata.EntrezGene.gz" \
    "${META_ENTREZ}" \
    "${ENTREZ_MD5}"

META_REFSEQ="${META_DIR}/gencode.vM38.metadata.RefSeq.gz"
REFSEQ_MD5=$(get_md5_from_manifest "gencode.vM38.metadata.RefSeq.gz")
download_and_verify \
    "${FTP}/gencode.vM38.metadata.RefSeq.gz" \
    "${META_REFSEQ}" \
    "${REFSEQ_MD5}"

# ─────────────────────────────────────────────
# 7. DECOMPRESS PHASE
# ─────────────────────────────────────────────

log "--- Phase 2: Decompress ---"

GENOME_FA="${GENOME_DIR}/GRCm39.primary_assembly.genome.fa"
ANNOT_GTF="${ANNOT_DIR}/gencode.vM38.primary_assembly.annotation.gtf"
TRNA_GTF="${ANNOT_DIR}/gencode.vM38.tRNAs.gtf"
TX_FA="${TX_DIR}/gencode.vM38.transcripts.fa"

decompress_if_needed "${GENOME_GZ}" "${GENOME_FA}"
decompress_if_needed "${ANNOT_GZ}"  "${ANNOT_GTF}"
decompress_if_needed "${TRNA_GZ}"   "${TRNA_GTF}"
decompress_if_needed "${TX_GZ}"     "${TX_FA}"

# ─────────────────────────────────────────────
# 8. GENOME INDEXING (samtools faidx)
# ─────────────────────────────────────────────

log "--- Phase 3: samtools faidx ---"

GENOME_FAI="${GENOME_FA}.fai"
if [[ -f "${GENOME_FAI}" ]]; then
    log "FAI index already exists, skipping."
else
    log "Building samtools faidx index..."
    samtools faidx "${GENOME_FA}"
    log "FAI index built: ${GENOME_FAI}"
fi

CHROM_SIZES="${GENOME_DIR}/GRCm39.primary_assembly.chrom.sizes"
if [[ -f "${CHROM_SIZES}" ]]; then
    log "chrom.sizes already exists, skipping."
else
    cut -f1,2 "${GENOME_FAI}" > "${CHROM_SIZES}"
    log "chrom.sizes written: ${CHROM_SIZES}"
fi

# ─────────────────────────────────────────────
# 9. FINAL SUMMARY
# ─────────────────────────────────────────────

log "======================================================"
log "Setup complete. Summary of key files:"
log ""
log "  Genome FASTA  : ${GENOME_FA}"
log "  Genome FAI    : ${GENOME_FAI}"
log "  Chrom sizes   : ${CHROM_SIZES}"
log "  Main GTF      : ${ANNOT_GTF}"
log "  tRNA GTF      : ${TRNA_GTF}"
log "  Transcripts   : ${TX_FA}"
log ""
log "  Metadata (compressed, gunzip on demand):"
log "    MGI symbols : ${META_MGI}"
log "    Entrez IDs  : ${META_ENTREZ}"
log "    RefSeq IDs  : ${META_REFSEQ}"
log ""
log "  Full log      : ${LOG}"
log "======================================================"
