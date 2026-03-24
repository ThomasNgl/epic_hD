#!/usr/bin/env bash
# =============================================================================
# 01_build_annotation.sh
#
# Purpose : Extract a flat gene annotation table from the GENCODE M38 GTF.
#           One row per gene. Saved alongside the GTF in the annotation folder.
#
# Output  : /mnt/auxiliary/gencode_primary/annotation/gencode.vM38.genes.tsv
#
#   Columns:
#     gene_id.          — stable Ensembl ID, version suffix stripped
#                         (ENSMUSG00000051951.6 → ENSMUSG00000051951)
#     gene_symbol       — MGI approved gene name (e.g. Xkr4)
#     gene_biotype      — transcript class (e.g. protein_coding, lncRNA,
#                         miRNA, snRNA, pseudogene, …)
#                         Useful to filter: keep protein_coding + lncRNA for
#                         long RNA-seq; keep miRNA + misc_RNA for small RNA
#     chr               — chromosome (UCSC-style: chr1 … chrM)
#     start             — 0-based start (BED convention)
#     end               — 1-based exclusive end
#     strand            — + or -
#
#   From this one file you can:
#     gene ID   → symbol    : awk '$1=="ENSMUSG..." {print $2}'
#     symbol    → gene IDs  : awk '$2=="Dnmt3l"  {print $1}'  (may return >1 row)
#     gene ID   → region    : awk '$1=="ENSMUSG..." {print $4,$5,$6}'
#     biotype filter        : awk '$3=="lncRNA"'
#
# Usage   : bash 01_build_annotation.sh
#
# Dependencies : awk (any POSIX awk)
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
OUT="${ANNOT_DIR}/gencode.vM38.genes.tsv"
LOG="${ANNOT_DIR}/01_build_annotation.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "======================================================"
log "01_build_annotation.sh — epic_hD Aim 1"
log "GTF : ${GTF}"
log "Out : ${OUT}"
log "======================================================"

# ─────────────────────────────────────────────
# 1. DEPENDENCY AND INPUT CHECK
# ─────────────────────────────────────────────

command -v awk &>/dev/null || { log "ERROR: awk not found."; exit 1; }

[[ -f "${GTF}" ]] || {
    log "ERROR: GTF not found: ${GTF}"
    log "       Run 00_download_gencode_m38.sh first."
    exit 1
}

# ─────────────────────────────────────────────
# 2. BUILD THE TABLE
# ─────────────────────────────────────────────

if [[ -f "${OUT}" ]]; then
    log "Output already exists, skipping: ${OUT}"
    exit 0
fi

log "Parsing GTF (gene features only)..."

awk -F'\t' '
BEGIN {
    OFS = "\t"
    print "gene_id", "gene_symbol", "gene_biotype", "chr", "start", "end", "strand"
}

# Skip comment lines
/^#/ { next }

# Process only gene-level features (one row per gene, no per-transcript duplicates)
$3 == "gene" {

    chr    = $1
    start  = $4 - 1    # convert GTF 1-based to BED 0-based
    end    = $5
    strand = $7

    gid     = ""
    sym     = "NA"      # gene_name in GTF = MGI-approved gene symbol
    biotype = "NA"      # gene_type in GTF = transcript class

    # Parse the semicolon-delimited attributes in field 9
    # GTF attribute format: key "value"; key "value";
    n = split($9, attrs, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+/, "", attrs[i])   # strip leading whitespace

        if (attrs[i] ~ /^gene_id /) {
            gsub(/gene_id "|"/, "", attrs[i])
            gid = attrs[i]
            gsub(/\.[0-9]+$/, "", gid)  # strip version suffix
        }
        if (attrs[i] ~ /^gene_name /) {
            gsub(/gene_name "|"/, "", attrs[i])
            sym = attrs[i]
        }
        if (attrs[i] ~ /^gene_type /) {
            gsub(/gene_type "|"/, "", attrs[i])
            biotype = attrs[i]
        }
    }

    # Require only gene_id — sym and biotype fall back to "NA" if absent
    if (gid != "")
        print gid, sym, biotype, chr, start, end, strand
}
' "${GTF}" > "${OUT}"

N_GENES=$(( $(wc -l < "${OUT}") - 1 ))
N_NA_SYM=$(awk 'NR>1 && $2=="NA"' "${OUT}" | wc -l)
N_NA_BIO=$(awk 'NR>1 && $3=="NA"' "${OUT}" | wc -l)
log "Done. ${N_GENES} genes written to ${OUT}"
[[ "${N_NA_SYM}" -gt 0 ]] && log "  WARNING: ${N_NA_SYM} gene(s) with no gene_name  → gene_symbol = NA"
[[ "${N_NA_BIO}" -gt 0 ]] && log "  WARNING: ${N_NA_BIO} gene(s) with no gene_type  → gene_biotype = NA"

# ─────────────────────────────────────────────
# 3. QUICK BIOTYPE SUMMARY
# ─────────────────────────────────────────────

log ""
log "Biotype breakdown:"
awk 'NR > 1 {print $3}' "${OUT}" \
    | sort | uniq -c | sort -rn \
    | awk '{printf "  %7s  %s\n", $1, $2}' \
    | tee -a "${LOG}"

log "======================================================"
log "Full log : ${LOG}"
log "======================================================"
