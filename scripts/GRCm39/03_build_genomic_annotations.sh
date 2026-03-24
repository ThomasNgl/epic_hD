#!/usr/bin/env bash
# =============================================================================
# 03_build_genomic_annotations.sh
#
# Purpose : Build all genomic annotation BED files needed for region-level
#           analysis: CpG context, gene structure, RepeatMasker classes,
#           and optionally ENCODE cCREs.
#
#           These BED files are the reference layer used by:
#             - annotate_regions.sh  (annotate any query BED per chromosome)
#             - region-level DMR aggregation in R (sgCpG2region_explore)
#             - any bedtools intersection against gene structure or CpG context
#
# Output  : /mnt/auxiliary/gencode_primary/genomic_annotations/
#
#   Gene structure (from GTF):
#     mm39_anno.genes.primary.bed         col4 = Ensembl gene ID (no version suffix)
#     mm39_anno.tss.primary.bed           col4 = gene_symbol:transcript_id (1 bp per transcript)
#     mm39_anno.promoters_2kb.primary.bed TSS ± 2kb, strand-aware, merged
#     mm39_anno.exons.primary.bed         all exons flattened and merged
#     mm39_anno.introns.primary.bed       gene_body minus exons
#     mm39_anno.intergenic.primary.bed    complement of all gene bodies
#
#   CpG context (from UCSC cpgIslandExt table):
#     mm39_anno.cpg_islands.primary.bed   CpG islands (UCSC definition)
#     mm39_anno.cpg_shores.primary.bed    ± 2kb from island edge, minus island
#     mm39_anno.cpg_shelves.primary.bed   2–4kb from island edge, minus shore
#     mm39_anno.cpg_opensea.primary.bed   genome minus islands+shores+shelves
#
#   RepeatMasker (from UCSC rmsk table):
#     mm39_anno.repeatmaskerDNA.primary.bed
#     mm39_anno.repeatmaskerLINE.primary.bed
#     mm39_anno.repeatmaskerLTR.primary.bed
#     mm39_anno.repeatmaskerSINE.primary.bed
#     mm39_anno.repeatmaskerSatellite.primary.bed
#     mm39_anno.repeatmaskerSimple_repeat.primary.bed
#     mm39_anno.repeatmaskerLow_complexity.primary.bed
#     mm39_anno.repeatmaskerNoRepeat.primary.bed
#
#   ENCODE cCREs (optional, requires liftOver):
#     mm39_anno.cCRE_ELS.enhancers.primary.bed   dELS + pELS enhancers
#     mm39_anno.cCRE_PLS.promoters.primary.bed   PLS promoter-like signatures
#     mm39_anno.cCRE_CTCF.primary.bed            CTCF-bound elements
#
#   Auxiliary:
#     mm39_anno.chrom.sizes                      all sequences
#     mm39_anno.chrom.sizes.primary              primary chromosomes only
#
# Usage   : bash 03_build_genomic_annotations.sh
#
# Dependencies : samtools, bedtools, awk, curl, wget
# Optional     : liftOver (ucsc-liftover) — for ENCODE cCRE section only
#
# Inputs  : /mnt/auxiliary/gencode_primary/
#               genome/GRCm39.primary_assembly.genome.fa   (+ .fai)
#               annotation/gencode.vM38.primary_assembly.annotation.gtf
#           (produced by 00_download_gencode_m38.sh)
#
# Author  : Thomas Negrello, negrello@hifo.uzh.ch
# Project : epic_hD — Aim 1 reference setup
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────

GENOME_FA="/mnt/auxiliary/gencode_primary/genome/GRCm39.primary_assembly.genome.fa"
GENOME_FAI="${GENOME_FA}.fai"
GTF="/mnt/auxiliary/gencode_primary/annotation/gencode.vM38.primary_assembly.annotation.gtf"

ANNOT_DIR="/mnt/auxiliary/gencode_primary/genomic_annotations"
OUT="${ANNOT_DIR}/mm39_anno"
LOG="${ANNOT_DIR}/03_build_genomic_annotations.log"

mkdir -p "${ANNOT_DIR}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "======================================================"
log "03_build_genomic_annotations.sh — epic_hD Aim 1"
log "GTF    : ${GTF}"
log "Output : ${ANNOT_DIR}"
log "======================================================"

# ─────────────────────────────────────────────
# 1. DEPENDENCY AND INPUT CHECK
# ─────────────────────────────────────────────

log "Checking dependencies..."
for bin in samtools bedtools awk curl; do
    command -v "${bin}" &>/dev/null || {
        log "ERROR: '${bin}' not found. Activate your conda environment."
        exit 1
    }
done

for f in "${GENOME_FA}" "${GENOME_FAI}" "${GTF}"; do
    [[ -f "${f}" ]] || {
        log "ERROR: Expected input not found: ${f}"
        log "       Run 00_download_gencode_m38.sh first."
        exit 1
    }
done
log "All inputs found."

# ─────────────────────────────────────────────
# 2. CHROMOSOME SIZES
#
# Two files:
#   chrom.sizes         — all sequences in the FASTA (chromosomes + scaffolds)
#   chrom.sizes.primary — filtered to canonical chromosomes only
#                         (chr1-19, chrX, chrY, chrM)
#
# The primary filter uses the regex ^chr([1-9]|1[0-9]|X|Y|M)$ which matches
# exactly chr1..chr19, chrX, chrY, chrM. The unfiltered file is needed as
# an intermediate by some bedtools steps; only the primary file is used for
# final outputs.
# ─────────────────────────────────────────────

log "--- Section 1: Chromosome sizes ---"

CHROMS="${OUT}.chrom.sizes"
CHROMS_PRI="${OUT}.chrom.sizes.primary"

if [[ -f "${CHROMS_PRI}" ]]; then
    log "chrom.sizes already exists, skipping."
else
    # FAI already built by 00_download_gencode_m38.sh — just extract cols 1-2
    cut -f1,2 "${GENOME_FAI}" > "${CHROMS}"
    awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${CHROMS}" > "${CHROMS_PRI}"
    log "chrom.sizes.primary: $(wc -l < "${CHROMS_PRI}") sequences"
fi

# ─────────────────────────────────────────────
# 3. GENE STRUCTURE FROM GTF
#
# All four sections parse the GTF with awk using -F'\t' to correctly
# split on tabs only (GTF attribute strings contain spaces which would
# shift field numbers with default FS).
#
# Key design decision vs original script:
#   The original used gene_name (symbol) as col4 in genes.bed.
#   Here we use the Ensembl gene ID (version stripped) as col4 so that
#   bedtools intersection output can be joined directly to the gene
#   annotation TSV (gencode.vM38.genes.tsv) by gene_id. The symbol is
#   in the TSV — no need to duplicate it in the BED.
# ─────────────────────────────────────────────

log "--- Section 2: Gene structure BED files ---"

# ── 3a. Genes ────────────────────────────────────────────────────────────────
# One row per gene. Coordinates are the full gene span (GTF start-1 to end).
#
# Columns:
#   col1 = chr
#   col2 = start (0-based)
#   col3 = end
#   col4 = gene_id    (Ensembl ID, version suffix stripped)
#   col5 = 0          (score placeholder, required by BED format)
#   col6 = strand
#   col7 = gene_symbol (MGI approved name, e.g. Dnmt3l)
#   col8 = gene_biotype (e.g. protein_coding, lncRNA)
#
# col4 = gene_id allows direct join to gencode.vM38.genes.tsv by key.
# col7 and col8 are carried so bedtools intersection output is self-contained
# without requiring a separate lookup.

if [[ -f "${OUT}.genes.primary.bed" ]]; then
    log "genes.primary.bed already exists, skipping."
else
    log "Building genes BED..."
    awk -F'\t' 'BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="gene" {
        gid="NA"; sym="NA"; bio="NA"; str=$7
        n=split($9,a,";")
        for(i=1;i<=n;i++){
            gsub(/^[ \t]+/,"",a[i])
            if(a[i] ~ /^gene_id /)   {gsub(/gene_id "|"/,"",a[i]);   gid=a[i]; gsub(/\.[0-9]+$/,"",gid)}
            if(a[i] ~ /^gene_name /) {gsub(/gene_name "|"/,"",a[i]); sym=a[i]}
            if(a[i] ~ /^gene_type /) {gsub(/gene_type "|"/,"",a[i]); bio=a[i]}
        }
        if(gid!="NA") print $1,$4-1,$5,gid,0,str,sym,bio
    }' "${GTF}" \
    | bedtools sort -faidx "${CHROMS}" -i - > "${OUT}.genes.bed"

    awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${OUT}.genes.bed" \
    > "${OUT}.genes.primary.bed"
    log "genes.primary.bed: $(wc -l < "${OUT}.genes.primary.bed") genes"
fi

# ── 3b. TSS (1 bp per transcript) ────────────────────────────────────────────
# TSS position depends on strand:
#   + strand: TSS = transcript start (GTF col4 - 1, 0-based)
#   - strand: TSS = transcript end   (GTF col5 - 1, 0-based)
#
# Columns:
#   col1 = chr
#   col2 = tss_start (0-based, 1bp)
#   col3 = tss_end
#   col4 = transcript_id  (ENSMUST..., version kept for uniqueness)
#   col5 = 0              (score placeholder)
#   col6 = strand
#   col7 = gene_symbol    (MGI name, e.g. Dnmt3l)
#   col8 = gene_id        (ENSMUSG..., version stripped, for joining to TSV)
#
# Storing transcript_id in col4 (not a concatenated string) makes each row
# uniquely identifiable. gene_symbol in col7 is what annotate_regions.sh
# will display as the nearest TSS label. gene_id in col8 allows joining back
# to the full annotation table.

if [[ -f "${OUT}.tss.primary.bed" ]]; then
    log "tss.primary.bed already exists, skipping."
else
    log "Building TSS BED..."
    awk -F'\t' 'BEGIN{OFS="\t"}
    /^#/ {next}
    $3=="transcript" {
        gid="NA"; sym="NA"; tid="NA"; str=$7; ts=$4-1; te=$5
        n=split($9,a,";")
        for(i=1;i<=n;i++){
            gsub(/^[ \t]+/,"",a[i])
            if(a[i] ~ /^gene_id /)       {gsub(/gene_id "|"/,"",a[i]);       gid=a[i]; gsub(/\.[0-9]+$/,"",gid)}
            if(a[i] ~ /^gene_name /)     {gsub(/gene_name "|"/,"",a[i]);     sym=a[i]}
            if(a[i] ~ /^transcript_id /) {gsub(/transcript_id "|"/,"",a[i]); tid=a[i]}
        }
        tss=(str=="+")?ts:(te-1)
        if(tid!="NA") print $1,tss,tss+1,tid,0,str,sym,gid
    }' "${GTF}" \
    | bedtools sort -faidx "${CHROMS}" -i - > "${OUT}.tss.bed"

    awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${OUT}.tss.bed" \
    > "${OUT}.tss.primary.bed"
    log "tss.primary.bed: $(wc -l < "${OUT}.tss.primary.bed") transcripts"
fi

# ── 3c. Promoters (TSS ± 2kb, strand-aware, merged) ──────────────────────────
# bedtools slop -s extends strand-aware: -l = upstream, -r = downstream.
# With -s, -l 2000 -r 2000 means 2kb upstream of TSS and 2kb downstream,
# respecting strand direction.
# After merging overlapping promoter windows:
#   col4 = distinct transcript_ids (comma-separated)
#   col5 = distinct strands

if [[ -f "${OUT}.promoters_2kb.primary.bed" ]]; then
    log "promoters_2kb.primary.bed already exists, skipping."
else
    log "Building promoters BED (TSS ± 2kb)..."
    bedtools slop -s -l 2000 -r 2000 \
        -i "${OUT}.tss.primary.bed" \
        -g "${CHROMS_PRI}" \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    | bedtools merge -i - -c 4,6 -o distinct,distinct \
    > "${OUT}.promoters_2kb.primary.bed"
    log "promoters_2kb.primary.bed: $(wc -l < "${OUT}.promoters_2kb.primary.bed") regions"
fi

# ── 3d. Exons (flattened and merged) ─────────────────────────────────────────
# All exon features across all transcripts are merged into a non-redundant
# flat exon set. This is used for gene structure classification:
# a region overlapping this file is exonic.

if [[ -f "${OUT}.exons.primary.bed" ]]; then
    log "exons.primary.bed already exists, skipping."
else
    log "Building exons BED..."
    awk -F'\t' 'BEGIN{OFS="\t"} /^#/{next} $3=="exon"{print $1,$4-1,$5}' "${GTF}" \
    | bedtools sort -faidx "${CHROMS}" -i - \
    | bedtools merge -i - \
    > "${OUT}.exons.bed"

    awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${OUT}.exons.bed" \
    > "${OUT}.exons.primary.bed"
    log "exons.primary.bed: $(wc -l < "${OUT}.exons.primary.bed") merged exon regions"
fi

# ── 3e. Introns (gene_body minus exons) ──────────────────────────────────────
# bedtools subtract removes exon intervals from gene body intervals.
# What remains is intronic sequence. Note: this is not per-transcript —
# it is the union of all intronic regions across all transcripts of a gene.

if [[ -f "${OUT}.introns.primary.bed" ]]; then
    log "introns.primary.bed already exists, skipping."
else
    log "Building introns BED (genes minus exons)..."
    bedtools subtract \
        -a "${OUT}.genes.primary.bed" \
        -b "${OUT}.exons.primary.bed" \
    > "${OUT}.introns.primary.bed"
    log "introns.primary.bed: $(wc -l < "${OUT}.introns.primary.bed") intronic regions"
fi

# ── 3f. Intergenic (complement of all gene bodies) ────────────────────────────
# bedtools complement returns all intervals NOT covered by the input.
# Input must be sorted; output covers the entire genome minus all gene spans.

if [[ -f "${OUT}.intergenic.primary.bed" ]]; then
    log "intergenic.primary.bed already exists, skipping."
else
    log "Building intergenic BED..."
    bedtools sort -faidx "${CHROMS_PRI}" -i "${OUT}.genes.primary.bed" \
    | bedtools complement -i - -g "${CHROMS_PRI}" \
    > "${OUT}.intergenic.primary.bed"
    log "intergenic.primary.bed: $(wc -l < "${OUT}.intergenic.primary.bed") intergenic regions"
fi

# ─────────────────────────────────────────────
# 4. CpG CONTEXT (from UCSC cpgIslandExt table)
#
# The four CpG context classes are mutually exclusive and cover the
# entire genome:
#
#   island   = CpG island as defined by UCSC (Gardner-Altman criteria:
#              length > 200bp, CpG obs/exp > 0.6, GC content > 50%)
#   shore    = ± 2kb from island boundary, minus the island itself
#   shelf    = 2–4kb from island boundary, minus shore
#   opensea  = everything else (background)
#
# Construction logic:
#   1. Download raw island table, filter to primary chr, sort → islands BED
#   2. Expand islands by 2kb → islands+2k; subtract islands → shores
#   3. Expand islands by 4kb → islands+4k; subtract islands+2k → shelves
#   4. Genome complement of (islands ∪ shores ∪ shelves) → opensea
#
# Two intermediate files (islands_plus2k, islands_plus4k) are kept because
# they are needed for the shelf calculation. They are deleted at the end.
# ─────────────────────────────────────────────

log "--- Section 3: CpG context BED files ---"

if [[ -f "${OUT}.cpg_opensea.primary.bed" ]]; then
    log "CpG context BEDs already exist, skipping."
else
    log "Downloading CpG island table from UCSC..."
    # UCSC cpgIslandExt table columns (0-indexed):
    # 0=bin 1=chrom 2=chromStart 3=chromEnd 4=name 5=length
    # 6=cpgNum 7=gcNum 8=perCpg 9=perGc 10=obsExp
    # We keep: chrom(1), start(2), end(3), name(4)
    curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cpgIslandExt.txt.gz" \
    | gunzip -c \
    | awk 'BEGIN{OFS="\t"} $2 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/ {print $2,$3,$4,$5}' \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.cpg_islands.primary.bed"
    log "cpg_islands.primary.bed: $(wc -l < "${OUT}.cpg_islands.primary.bed") islands"

    # Shores: expand 2kb each side, subtract islands
    bedtools slop -g "${CHROMS_PRI}" -b 2000 \
        -i "${OUT}.cpg_islands.primary.bed" \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.cpg_islands_plus2k.primary.bed"

    bedtools subtract \
        -a "${OUT}.cpg_islands_plus2k.primary.bed" \
        -b "${OUT}.cpg_islands.primary.bed" \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.cpg_shores.primary.bed"
    log "cpg_shores.primary.bed: $(wc -l < "${OUT}.cpg_shores.primary.bed") shore regions"

    # Shelves: expand 4kb, subtract islands+shores (= islands_plus2k)
    bedtools slop -g "${CHROMS_PRI}" -b 4000 \
        -i "${OUT}.cpg_islands.primary.bed" \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.cpg_islands_plus4k.primary.bed"

    bedtools subtract \
        -a "${OUT}.cpg_islands_plus4k.primary.bed" \
        -b "${OUT}.cpg_islands_plus2k.primary.bed" \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.cpg_shelves.primary.bed"
    log "cpg_shelves.primary.bed: $(wc -l < "${OUT}.cpg_shelves.primary.bed") shelf regions"

    # Opensea: genome minus all three classes
    # Build the union of island+shore+shelf, merge overlaps, then complement
    bedtools subtract \
        -a <(awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "${CHROMS_PRI}") \
        -b <(cat \
                "${OUT}.cpg_islands.primary.bed" \
                "${OUT}.cpg_shores.primary.bed" \
                "${OUT}.cpg_shelves.primary.bed" \
             | bedtools sort -faidx "${CHROMS_PRI}" -i - \
             | bedtools merge -i -) \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.cpg_opensea.primary.bed"
    log "cpg_opensea.primary.bed: $(wc -l < "${OUT}.cpg_opensea.primary.bed") opensea regions"

    # Remove intermediates
    rm -f "${OUT}.cpg_islands_plus2k.primary.bed" \
          "${OUT}.cpg_islands_plus4k.primary.bed"
fi

# ─────────────────────────────────────────────
# 5. REPEATMASKER (from UCSC rmsk table)
#
# The UCSC rmsk table contains all repeat annotations from RepeatMasker.
# We download it, filter to primary chromosomes, split by repeat class,
# and build a NoRepeat file (genome minus all repeats).
#
# UCSC rmsk table columns (relevant ones, 0-indexed):
#   5=genoName(chr)  6=genoStart  7=genoEnd
#   10=repName  11=repClass  12=repFamily
# After filtering: col1=chr, col2=start, col3=end, col4=repName, col5=repClass, col6=repFamily
#
# Classes split: DNA, LINE, LTR, SINE, Satellite, Simple_repeat, Low_complexity
# Note the original used $5==c to match repClass (column 5 in the output BED).
# ─────────────────────────────────────────────

log "--- Section 4: RepeatMasker BED files ---"

if [[ -f "${OUT}.repeatmaskerNoRepeat.primary.bed" ]]; then
    log "RepeatMasker BEDs already exist, skipping."
else
    log "Downloading RepeatMasker table from UCSC (this is large, ~10 min)..."
    curl -s "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz" \
    | gunzip -c \
    | awk -v OFS='\t' '$6 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/ {print $6,$7,$8,$11,$12,$13}' \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.repeatmasker.all.primary.bed"
    log "repeatmasker.all: $(wc -l < "${OUT}.repeatmasker.all.primary.bed") elements"

    # Split by class (col5 = repClass in the output BED)
    for cls in DNA LINE LTR SINE Satellite Simple_repeat Low_complexity; do
        awk -v c="${cls}" 'BEGIN{OFS="\t"} $5==c {print $1,$2,$3,$4,$5,$6}' \
            "${OUT}.repeatmasker.all.primary.bed" \
        > "${OUT}.repeatmasker${cls}.primary.bed"
        log "repeatmasker${cls}: $(wc -l < "${OUT}.repeatmasker${cls}.primary.bed") elements"
    done

    # NoRepeat: genome minus all repeat-annotated intervals
    bedtools subtract \
        -a <(awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "${CHROMS_PRI}") \
        -b "${OUT}.repeatmasker.all.primary.bed" \
    | bedtools sort -faidx "${CHROMS_PRI}" -i - \
    > "${OUT}.repeatmaskerNoRepeat.primary.bed"
    log "repeatmaskerNoRepeat: $(wc -l < "${OUT}.repeatmaskerNoRepeat.primary.bed") non-repeat regions"
fi

# ─────────────────────────────────────────────
# 6. ENCODE cCREs (optional — requires liftOver)
#
# ENCODE candidate cis-regulatory elements (cCREs) for mouse are only
# available mapped to mm10. We download the mm10 BED, liftOver to mm39
# using the UCSC chain file, then split by cCRE class:
#   ELS  (dELS + pELS) = enhancer-like signatures
#   PLS                = promoter-like signatures
#   CTCF               = CTCF-bound elements
#
# If liftOver is not in PATH, this section is skipped with a note.
# Install with: conda install -c bioconda ucsc-liftover
# ─────────────────────────────────────────────

log "--- Section 5: ENCODE cCREs (optional) ---"

if ! command -v liftOver &>/dev/null; then
    log "NOTE: liftOver not found — skipping cCRE section."
    log "      Install with: conda install -c bioconda ucsc-liftover"
else
    if [[ -f "${OUT}.cCRE_ELS.enhancers.primary.bed" ]]; then
        log "cCRE BEDs already exist, skipping."
    else
        CCRE_MM10_GZ="${ANNOT_DIR}/mm10_cCREs.bed.gz"
        CCRE_MM10="${ANNOT_DIR}/mm10_cCREs.bed"
        CHAIN="${ANNOT_DIR}/mm10ToMm39.over.chain.gz"
        CCRE_MM39_RAW="${ANNOT_DIR}/mm39_cCREs.lifted.bed"
        CCRE_UNMAP="${ANNOT_DIR}/mm39_cCREs.unmapped.bed"
        CCRE_MM39_ALL="${ANNOT_DIR}/mm39_cCREs.all.bed"

        log "Downloading mm10 cCREs from ENCODE..."
        curl -L -o "${CCRE_MM10_GZ}" \
            "https://www.encodeproject.org/files/ENCFF167FJQ/@@download/ENCFF167FJQ.bed.gz"
        gunzip -c "${CCRE_MM10_GZ}" > "${CCRE_MM10}"

        log "Downloading mm10→mm39 liftOver chain..."
        curl -L -o "${CHAIN}" \
            "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"

        log "Running liftOver mm10 → mm39..."
        liftOver -bedPlus=9 "${CCRE_MM10}" "${CHAIN}" "${CCRE_MM39_RAW}" "${CCRE_UNMAP}"

        # Filter to primary chromosomes and sort
        # After liftOver the cCRE label is in column 10
        awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${CCRE_MM39_RAW}" \
        | bedtools sort -faidx "${CHROMS_PRI}" -i - \
        > "${CCRE_MM39_ALL}"

        # ELS enhancers (distal + proximal)
        awk 'BEGIN{OFS="\t"} ($10=="dELS" || $10=="pELS") {print $1,$2,$3}' "${CCRE_MM39_ALL}" \
        | bedtools sort -faidx "${CHROMS_PRI}" -i - \
        | bedtools merge -i - \
        > "${OUT}.cCRE_ELS.enhancers.primary.bed"

        # PLS promoter-like
        awk 'BEGIN{OFS="\t"} $10=="PLS" {print $1,$2,$3}' "${CCRE_MM39_ALL}" \
        | bedtools sort -faidx "${CHROMS_PRI}" -i - \
        | bedtools merge -i - \
        > "${OUT}.cCRE_PLS.promoters.primary.bed"

        # CTCF-bound
        awk 'BEGIN{OFS="\t"} $10 ~ /CTCF/ {print $1,$2,$3}' "${CCRE_MM39_ALL}" \
        | bedtools sort -faidx "${CHROMS_PRI}" -i - \
        | bedtools merge -i - \
        > "${OUT}.cCRE_CTCF.primary.bed"

        log "cCRE_ELS.enhancers : $(wc -l < "${OUT}.cCRE_ELS.enhancers.primary.bed") regions"
        log "cCRE_PLS.promoters : $(wc -l < "${OUT}.cCRE_PLS.promoters.primary.bed") regions"
        log "cCRE_CTCF          : $(wc -l < "${OUT}.cCRE_CTCF.primary.bed") regions"

        # Remove download intermediates (chain and raw liftover kept for reproducibility)
        rm -f "${CCRE_MM10}" "${CCRE_MM10_GZ}"
    fi
fi

# ─────────────────────────────────────────────
# 7. FINAL SORT PASS
#
# Ensures all output BEDs satisfy bedtools -sorted requirements.
# This is a safety pass — individual sections already sort their outputs,
# but re-sorting guarantees consistency if any section was skipped or
# partially completed from a previous interrupted run.
# ─────────────────────────────────────────────

log "--- Section 6: Final sort pass ---"

shopt -s nullglob
for f in "${ANNOT_DIR}"/mm39_anno.*.primary.bed; do
    [[ -e "${f}" ]] || continue
    tmp="${f}.sorted.tmp"
    if bedtools sort -faidx "${CHROMS_PRI}" -i "${f}" > "${tmp}" 2>/dev/null; then
        mv "${tmp}" "${f}"
    else
        log "WARNING: Could not sort $(basename "${f}") — skipping"
        rm -f "${tmp}"
    fi
done
shopt -u nullglob
log "Sort pass complete."

# ─────────────────────────────────────────────
# 8. SANITY CHECK — LINE COUNTS
# ─────────────────────────────────────────────

log "--- Summary: line counts ---"
for f in \
    genes.primary.bed tss.primary.bed promoters_2kb.primary.bed \
    exons.primary.bed introns.primary.bed intergenic.primary.bed \
    cpg_islands.primary.bed cpg_shores.primary.bed \
    cpg_shelves.primary.bed cpg_opensea.primary.bed \
    repeatmaskerDNA.primary.bed repeatmaskerLINE.primary.bed \
    repeatmaskerLTR.primary.bed repeatmaskerSINE.primary.bed \
    repeatmaskerSatellite.primary.bed repeatmaskerSimple_repeat.primary.bed \
    repeatmaskerLow_complexity.primary.bed repeatmaskerNoRepeat.primary.bed
do
    target="${OUT}.${f}"
    [[ -f "${target}" ]] \
        && log "  $(printf '%8d' $(wc -l < "${target}"))  ${f}" \
        || log "  MISSING: ${f}"
done

log "======================================================"
log "Done. All annotation BEDs written to: ${ANNOT_DIR}"
log "Full log : ${LOG}"
log "======================================================"
