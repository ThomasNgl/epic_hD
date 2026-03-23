#!/usr/bin/env bash
# Annotate regions (mm39): primary category + CpG context + RepeatMasker + nearest TSS (name/id/distance)
# Writes a single final TSV **with header**. Intermediates are deleted unless KEEP_TMP=1.

# Usage:
#   bash af_annotate_regions.sh \
#     /mnt/thomas/matriline/results/GRCm39/annotations \
#     /path/to/regions.bed \
#     chr1 \
#     /path/to/outdir

set -euo pipefail
# Set default for temporary file cleanup.
# User can override this by exporting KEEP_TMP=1 before running.
: ${KEEP_TMP:=0}

# Display help message and exit if arguments are incorrect
if [ "$#" -ne 4 ]; then
  echo "--- Genomic Region Annotator (mm39) ---" >&2
  echo "Purpose: Annotates genomic regions with primary category, CpG context," >&2
  echo "         RepeatMasker class, and nearest TSS features." >&2
  echo " " >&2
  echo "Usage: $0 <ANNOT_DIR> <REGIONS_BED> <CHR> <OUT_DIR>" >&2
  echo " " >&2
  echo "Arguments:" >&2
  echo "  <ANNOT_DIR>  : Directory containing required annotation BED files." >&2
  echo "  <REGIONS_BED>: Input BED file of regions to annotate (e.g., peaks.bed)." >&2
  echo "  <CHR>        : Target chromosome (e.g., 'chr1', 'chrX')." >&2
  echo "  <OUT_DIR>    : Output directory for final TSV and summaries." >&2
  echo " " >&2
  echo "Optional: Set KEEP_TMP=1 environment variable to retain intermediate files." >&2
  exit 1
fi

ANNOT_DIR="$(readlink -f "$1")"
REGS="$(readlink -f "$2")"
CHR="$3"
OUT_DIR="$(readlink -f "$4")"
mkdir -p "$OUT_DIR"

# References
CHROMS="${ANNOT_DIR}/mm39_anno.chrom.sizes.primary"
GENE="${ANNOT_DIR}/mm39_anno.genes.primary.bed"
TSS="${ANNOT_DIR}/mm39_anno.tss.primary.bed"
PROM="${ANNOT_DIR}/mm39_anno.promoters_2kb.primary.bed"

# Optional enhancers (ENCODE cCRE ELS)
ENH="${ANNOT_DIR}/mm39_anno.cCRE_ELS.enhancers.primary.bed"
HAS_ENH=0; [ -s "$ENH" ] && HAS_ENH=1

# Gene structure
EXON="${ANNOT_DIR}/mm39_anno.exons.primary.bed"
INTRON="${ANNOT_DIR}/mm39_anno.introns.primary.bed"

# CpG context
CG_ISL="${ANNOT_DIR}/mm39_anno.cpg_islands.primary.bed"
CG_SHO="${ANNOT_DIR}/mm39_anno.cpg_shores.primary.bed"
CG_SHE="${ANNOT_DIR}/mm39_anno.cpg_shelves.primary.bed"
CG_OPEN="${ANNOT_DIR}/mm39_anno.cpg_opensea.primary.bed"

# RepeatMasker
RM_DNA="${ANNOT_DIR}/mm39_anno.repeatmaskerDNA.primary.bed"
RM_LINE="${ANNOT_DIR}/mm39_anno.repeatmaskerLINE.primary.bed"
RM_LTR="${ANNOT_DIR}/mm39_anno.repeatmaskerLTR.primary.bed"
RM_SINE="${ANNOT_DIR}/mm39_anno.repeatmaskerSINE.primary.bed"
RM_SAT="${ANNOT_DIR}/mm39_anno.repeatmaskerSatellite.primary.bed"
RM_SIM="${ANNOT_DIR}/mm39_anno.repeatmaskerSimple_repeat.primary.bed"
RM_LOW="${ANNOT_DIR}/mm39_anno.repeatmaskerLow_complexity.primary.bed"
RM_NOREP="${ANNOT_DIR}/mm39_anno.repeatmaskerNoRepeat.primary.bed"

BASE="$(basename "${REGS%.*}")"
RCHR="${OUT_DIR}/${BASE}.processed"

# Regions → same chr only, then sort by genome order
awk -v C="$CHR" 'BEGIN{OFS="\t"} $1==C{print $1,$2,$3}' "$REGS" > "${RCHR}.bed"
bedtools sort -faidx "$CHROMS" -i "${RCHR}.bed" > "${RCHR}.sorted.bed"

echo "== debug overlaps =="
for B in "$PROM" "$ENH" "$EXON" "$INTRON" "$GENE" "$CG_ISL" "$CG_SHO" "$CG_SHE" "$CG_OPEN"; do
  [ -e "$B" ] || continue
  n=$(bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$B" | wc -l)
  echo "  $(basename "$B"): $n"
done

########################
# 1) PRIMARY CATEGORY  #
########################
# Priority: promoter > enhancer > exon > intron > gene_body > intergenic
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$PROM" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"promoter"}' > "${RCHR}.cat.prom"

if [ $HAS_ENH -eq 1 ]; then
  bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$ENH" \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"enhancer"}' > "${RCHR}.cat.enh"
else
  : > "${RCHR}.cat.enh"
fi

bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$EXON" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"exon"}' > "${RCHR}.cat.exon"

bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$INTRON" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"intron"}' > "${RCHR}.cat.intron"

bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$GENE" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"gene_body"}' > "${RCHR}.cat.gene"

cat "${RCHR}.cat.prom" "${RCHR}.cat.enh" "${RCHR}.cat.exon" "${RCHR}.cat.intron" "${RCHR}.cat.gene" \
| awk 'BEGIN{OFS="\t"}{k=$1 FS $2 FS $3; if(!(k in seen)){print; seen[k]=1}}' > "${RCHR}.cat.priority"

# Primary category table (default to intergenic when no hits)
if [ ! -s "${RCHR}.cat.priority" ]; then
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"intergenic"}' "${RCHR}.sorted.bed" > "${RCHR}.primary_category.tsv"
else
  awk 'BEGIN{OFS="\t"} NR==FNR{ann[$1 FS $2 FS $3]=$4; next}
       {k=$1 FS $2 FS $3; print $1,$2,$3, ((k in ann)?ann[k]:"intergenic")}' \
      "${RCHR}.cat.priority" "${RCHR}.sorted.bed" > "${RCHR}.primary_category.tsv"
fi

########################
# 2) CpG CONTEXT       #
########################
# island > shore > shelf > opensea (mutually exclusive by design)
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$CG_ISL" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"CpG_island"}' > "${RCHR}.cpg.isl"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$CG_SHO" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"CpG_shore"}' > "${RCHR}.cpg.sho"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$CG_SHE" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"CpG_shelf"}' > "${RCHR}.cpg.she"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$CG_OPEN" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"CpG_opensea"}' > "${RCHR}.cpg.open"

cat "${RCHR}.cpg.isl" "${RCHR}.cpg.sho" "${RCHR}.cpg.she" "${RCHR}.cpg.open" \
| awk 'BEGIN{OFS="\t"}{k=$1 FS $2 FS $3; if(!(k in seen)){print; seen[k]=1}}' > "${RCHR}.cpg.priority"

awk 'BEGIN{OFS="\t"} NR==FNR{ctx[$1 FS $2 FS $3]=$4; next}
     {k=$1 FS $2 FS $3; print $1,$2,$3, ((k in ctx)?ctx[k]:"CpG_opensea")}' \
    "${RCHR}.cpg.priority" "${RCHR}.sorted.bed" > "${RCHR}.cpg_context.tsv"

########################
# 3) REPEAT CLASS      #
########################
# First matching class wins; fallback to NoRepeat (if in non-repetitive set) else "None"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_DNA"  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"DNA"}'             > "${RCHR}.rep.dna"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_LINE" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"LINE"}'            > "${RCHR}.rep.line"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_LTR"  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"LTR"}'             > "${RCHR}.rep.ltr"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_SINE" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"SINE"}'            > "${RCHR}.rep.sine"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_SAT"  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"Satellite"}'       > "${RCHR}.rep.sat"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_SIM"  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"Simple_repeat"}'   > "${RCHR}.rep.sim"
bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_LOW"  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"Low_complexity"}'  > "${RCHR}.rep.low"

cat "${RCHR}.rep.dna" "${RCHR}.rep.line" "${RCHR}.rep.ltr" "${RCHR}.rep.sine" "${RCHR}.rep.sat" "${RCHR}.rep.sim" "${RCHR}.rep.low" \
| awk 'BEGIN{OFS="\t"}{k=$1 FS $2 FS $3; if(!(k in seen)){print; seen[k]=1}}' > "${RCHR}.rep.priority"

bedtools intersect -sorted -g "$CHROMS" -u -a "${RCHR}.sorted.bed" -b "$RM_NOREP" \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"NoRepeat"}' > "${RCHR}.rep.norep"

awk 'BEGIN{OFS="\t"} NR==FNR{cls[$1 FS $2 FS $3]=$4; next}
     {k=$1 FS $2 FS $3; if(k in cls){print $1,$2,$3,cls[k]} else {print $1,$2,$3,"None"}}' \
    "${RCHR}.rep.priority" "${RCHR}.sorted.bed" > "${RCHR}.repeat_class.tmp"

awk 'BEGIN{OFS="\t"} NR==FNR{nr[$1 FS $2 FS $3]=1; next}
     {k=$1 FS $2 FS $3; if($4=="None" && (k in nr)) $4="NoRepeat"; print}' \
    "${RCHR}.rep.norep" "${RCHR}.repeat_class.tmp" > "${RCHR}.repeat_class.tsv"

########################################
# 4) NEAREST TSS (name, id, distance)  #
########################################
bedtools closest -sorted -g "$CHROMS" -a "${RCHR}.sorted.bed" -b "$TSS" -D a -t first -wb \
| awk 'BEGIN{OFS="\t"}{
         name=$7; dist=$NF; tid=name;
         n=split(name,parts,":"); if(n>=2){tid=parts[n]}
         print $1,$2,$3,name,tid,dist
       }' > "${RCHR}.nearest_tss_name_id_dist.tsv"

##############################
# 5) FINAL MERGE + SUMMARY   #
##############################
# Final table **with header** (single file)
{
  echo -e "chr\tstart\tend\tprimary_category\tcpg_context\trepeat_class\tnearest_TSS_name\tnearest_TSS_id\tnearest_TSS_distance"
  awk 'BEGIN{OFS="\t"} NR==FNR{pc[$1 FS $2 FS $3]=$4; next}
       {k=$1 FS $2 FS $3; print $1,$2,$3, ((k in pc)?pc[k]:"intergenic")}' \
     "${RCHR}.primary_category.tsv" "${RCHR}.sorted.bed" \
  | awk 'BEGIN{OFS="\t"} NR==FNR{cg[$1 FS $2 FS $3]=$4; next}
         {k=$1 FS $2 FS $3; print $0, ((k in cg)?cg[k]:"CpG_opensea")}' \
        "${RCHR}.cpg_context.tsv" - \
  | awk 'BEGIN{OFS="\t"} NR==FNR{rp[$1 FS $2 FS $3]=$4; next}
         {k=$1 FS $2 FS $3; print $0, ((k in rp)?rp[k]:"None")}' \
        "${RCHR}.repeat_class.tsv" - \
  | awk 'BEGIN{OFS="\t"} NR==FNR{n[$1 FS $2 FS $3]=$4; i[$1 FS $2 FS $3]=$5; d[$1 FS $2 FS $3]=$6; next}
         {k=$1 FS $2 FS $3; print $0, ((k in n)?n[k]:"NA"), ((k in i)?i[k]:"NA"), ((k in d)?d[k]:"NA")}' \
        "${RCHR}.nearest_tss_name_id_dist.tsv" -
} > "${RCHR}.annotated.tsv"

# Summaries (kept)
cut -f4 "${RCHR}.annotated.tsv" | tail -n +2 | sort | uniq -c > "${RCHR}.primary_category_counts.txt"
cut -f5 "${RCHR}.annotated.tsv" | tail -n +2 | sort | uniq -c > "${RCHR}.cpg_context_counts.txt"
cut -f6 "${RCHR}.annotated.tsv" | tail -n +2 | sort | uniq -c > "${RCHR}.repeat_class_counts.txt"

# Cleanup (delete intermediates unless KEEP_TMP=1)
if [ "${KEEP_TMP:-0}" -ne 1 ]; then
  rm -f "${RCHR}.bed" "${RCHR}.sorted.bed" \
        "${RCHR}.cat."* "${RCHR}.cpg."* "${RCHR}.rep."* \
        "${RCHR}.cat.priority" "${RCHR}.cpg.priority" "${RCHR}.repeat_class.tmp" \
        "${RCHR}.primary_category.tsv" "${RCHR}.cpg_context.tsv" "${RCHR}.repeat_class.tsv" \
        "${RCHR}.nearest_tss_name_id_dist.tsv"
fi

echo "Wrote:"
echo "  ${RCHR}.annotated.tsv"
echo "  ${RCHR}.primary_category_counts.txt"
echo "  ${RCHR}.cpg_context_counts.txt"
echo "  ${RCHR}.repeat_class_counts.txt"

# Add friendly symbol column from nearest_TSS_name (token before ':')
awk 'BEGIN{OFS="\t"} NR==1{print $0,"nearest_TSS_symbol"; next} {sym=$7; split($7,a,":"); if(a[1]!="") sym=a[1]; print $0,sym}' \
  "${RCHR}.annotated.tsv" > "${RCHR}.annotated.with_symbol.tsv" \
&& mv "${RCHR}.annotated.with_symbol.tsv" "${RCHR}.annotated.tsv"
