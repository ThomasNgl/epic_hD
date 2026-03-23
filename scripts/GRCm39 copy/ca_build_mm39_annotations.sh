#!/usr/bin/env bash
# Build mm39 annotations: genes/TSS/promoters + CpG (island/shore/shelf/opensea) +
# exons/introns/intergenic + RepeatMasker + cCREs (mm10→mm39).

# TO RUN:
#   bash ca_build_mm39_annotations.sh /mnt/thomas/matriline/results/GRCm39/gencode_primary

set -euo pipefail

# --- I/O & required inputs ---
GP_DIR="${1:?Usage: $0 /path/to/gencode_primary}"
GP_DIR="$(readlink -f "$GP_DIR")"
BASE_DIR="$(dirname "$GP_DIR")"
ANNOT_DIR="${BASE_DIR}/annotations"
GENOME="${GP_DIR}/GRCm39.primary_assembly.genome.fa"
GTF="${GP_DIR}/gencode.vM37.primary_assembly.annotation.gtf"
OUT="${ANNOT_DIR}/mm39_anno"

mkdir -p "$ANNOT_DIR"

# --- Dependencies sanity check ---
for bin in samtools bedtools curl awk; do
  command -v "$bin" >/dev/null || { echo "ERROR: $bin not found in PATH"; exit 1; }
done

# --- Chrom sizes (all + primary) ---
samtools faidx "$GENOME"
cut -f1,2 "${GENOME}.fai" > "${OUT}.chrom.sizes"
awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${OUT}.chrom.sizes" > "${OUT}.chrom.sizes.primary"

# --- Genes (BED) ---
awk 'BEGIN{OFS="\t"} $3=="gene"{s=$4-1;e=$5;gid="";gname="";str=$7;
     if(match($0,/gene_id "([^"]+)"/,a))gid=a[1];
     if(match($0,/gene_name "([^"]+)"/,b))gname=b[1];
     nm=(gname!=""?gname:gid);print $1,s,e,nm,0,str}' "$GTF" \
| bedtools sort -faidx "${OUT}.chrom.sizes" -i - > "${OUT}.genes.bed"
awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${OUT}.genes.bed" > "${OUT}.genes.primary.bed"

# --- TSS (1 bp per transcript) ---
awk 'BEGIN{OFS="\t"} $3=="transcript"{ts=$4-1;te=$5;str=$7;gid="";gname="";tid="";
     if(match($0,/gene_id "([^"]+)"/,a))gid=a[1];
     if(match($0,/gene_name "([^"]+)"/,b))gname=b[1];
     if(match($0,/transcript_id "([^"]+)"/,c))tid=c[1];
     nm=(gname!=""?gname:gid)":"tid;
     tss=(str=="+")?ts:(te-1);print $1,tss,tss+1,nm,0,str}' "$GTF" \
| bedtools sort -faidx "${OUT}.chrom.sizes" -i - > "${OUT}.tss.bed"
awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${OUT}.tss.bed" > "${OUT}.tss.primary.bed"

# --- Promoters (TSS ±2 kb; strand-aware; merged) ---
bedtools slop -s -l 2000 -r 2000 -i "${OUT}.tss.primary.bed" -g "${OUT}.chrom.sizes.primary" \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - \
| bedtools merge -i - -c 4,6 -o distinct,distinct > "${OUT}.promoters_2kb.primary.bed"

# --- CpG islands (UCSC), shores (±2kb minus islands), shelves (2–4kb), opensea (genome minus all) ---
curl -s https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cpgIslandExt.txt.gz | gunzip -c \
| awk 'BEGIN{OFS="\t"} $2 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/{print $2,$3,$4,$5}' \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "${OUT}.cpg_islands.primary.bed"

bedtools slop -g "${OUT}.chrom.sizes.primary" -b 2000 -i "${OUT}.cpg_islands.primary.bed" \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "${OUT}.cpg_islands_plus2k.primary.bed"

bedtools subtract -a "${OUT}.cpg_islands_plus2k.primary.bed" -b "${OUT}.cpg_islands.primary.bed" \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "${OUT}.cpg_shores.primary.bed"

bedtools slop -g "${OUT}.chrom.sizes.primary" -b 4000 -i "${OUT}.cpg_islands.primary.bed" \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "${OUT}.cpg_islands_plus4k.primary.bed"

bedtools subtract -a "${OUT}.cpg_islands_plus4k.primary.bed" -b "${OUT}.cpg_islands_plus2k.primary.bed" \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "${OUT}.cpg_shelves.primary.bed"

bedtools subtract -a <(awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "${OUT}.chrom.sizes.primary") \
                  -b <(cat "${OUT}.cpg_islands.primary.bed" "${OUT}.cpg_shores.primary.bed" "${OUT}.cpg_shelves.primary.bed" \
                      | bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - | bedtools merge -i -) \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "${OUT}.cpg_opensea.primary.bed"

# --- Exons (flattened), introns (genes−exons), intergenic (complement of genes) ---
awk 'BEGIN{OFS="\t"} $3=="exon"{print $1,$4-1,$5}' "$GTF" \
| bedtools sort -faidx "${OUT}.chrom.sizes" -i - | bedtools merge -i - > "${OUT}.exons.bed"
awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "${OUT}.exons.bed" > "${OUT}.exons.primary.bed"
bedtools subtract -a "${OUT}.genes.primary.bed" -b "${OUT}.exons.primary.bed" > "${OUT}.introns.primary.bed"
bedtools complement -i <(bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i "${OUT}.genes.primary.bed") -g "${OUT}.chrom.sizes.primary" \
> "${OUT}.intergenic.primary.bed"

# --- RepeatMasker (UCSC rmsk): all + split by class + NoRepeat ---
curl -s https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz | gunzip -c \
| awk -vOFS='\t' '$6 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/ {print $6,$7,$8,$11,$12,$13}' \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - \
> "${OUT}.repeatmasker.all.primary.bed"

for cls in DNA LINE LTR SINE Satellite Simple_repeat Low_complexity; do
  awk -v c="$cls" 'BEGIN{OFS="\t"} $5==c {print $1,$2,$3,$4,$5,$6}' "${OUT}.repeatmasker.all.primary.bed" \
  > "${OUT}.repeatmasker${cls}.primary.bed"
done

bedtools subtract -a <(awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "${OUT}.chrom.sizes.primary") -b "${OUT}.repeatmasker.all.primary.bed" \
| bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "${OUT}.repeatmaskerNoRepeat.primary.bed"

# --- ENCODE cCREs: download mm10, liftOver to mm39, keep primaries, split (dELS/pELS) / PLS / CTCF ---
if command -v liftOver >/dev/null 2>&1; then
  CCRE_MM10_GZ="${ANNOT_DIR}/mm10_cCREs.bed.gz"
  CCRE_MM10="${ANNOT_DIR}/mm10_cCREs.bed"
  CHAIN="${ANNOT_DIR}/mm10ToMm39.over.chain.gz"
  CCRE_MM39_RAW="${ANNOT_DIR}/mm39_cCREs.lifted.bed"
  CCRE_UNMAP="${ANNOT_DIR}/mm39_cCREs.unmapped.bed"
  CCRE_MM39_ALL="${ANNOT_DIR}/mm39_cCREs.all.bed"

  curl -L -o "$CCRE_MM10_GZ" "https://www.encodeproject.org/files/ENCFF167FJQ/@@download/ENCFF167FJQ.bed.gz"
  gunzip -c "$CCRE_MM10_GZ" > "$CCRE_MM10"
  curl -L -o "$CHAIN" "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz"
  liftOver -bedPlus=9 "$CCRE_MM10" "$CHAIN" "$CCRE_MM39_RAW" "$CCRE_UNMAP"

  awk '$1 ~ /^chr([1-9]|1[0-9]|X|Y|M)$/' "$CCRE_MM39_RAW" \
  | bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - > "$CCRE_MM39_ALL"

  # NOTE: ENCODE cCRE label is in column 10 after liftOver (e.g. dELS, pELS, PLS, CA-CTCF, ...)
  awk 'BEGIN{OFS="\t"} ($10=="dELS" || $10=="pELS") {print $1,$2,$3}' "$CCRE_MM39_ALL" \
  | bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - | bedtools merge -i - \
  > "${OUT}.cCRE_ELS.enhancers.primary.bed"

  awk 'BEGIN{OFS="\t"} $10=="PLS" {print $1,$2,$3}' "$CCRE_MM39_ALL" \
  | bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - | bedtools merge -i - \
  > "${OUT}.cCRE_PLS.promoters.primary.bed"

  awk 'BEGIN{OFS="\t"} $10 ~ /CTCF/ {print $1,$2,$3}' "$CCRE_MM39_ALL" \
  | bedtools sort -faidx "${OUT}.chrom.sizes.primary" -i - | bedtools merge -i - \
  > "${OUT}.cCRE_CTCF.primary.bed"

  echo -e "$(wc -l < "$CCRE_MM39_ALL")\t$CCRE_MM39_ALL"
  echo -e "$(wc -l < "${OUT}.cCRE_ELS.enhancers.primary.bed")\t${OUT}.cCRE_ELS.enhancers.primary.bed"
  echo -e "$(wc -l < "${OUT}.cCRE_PLS.promoters.primary.bed")\t${OUT}.cCRE_PLS.promoters.primary.bed"
  echo -e "$(wc -l < "${OUT}.cCRE_CTCF.primary.bed")\t${OUT}.cCRE_CTCF.primary.bed"
else
  echo "NOTE: liftOver not found; skipping cCRE liftOver (install: conda install -c bioconda ucsc-liftover)"
fi

# Resort all reference BEDs to satisfy -sorted operations (fail hard if sort fails)
CHROMS="${OUT}.chrom.sizes.primary"
shopt -s nullglob

for f in \
  ${ANNOT_DIR}/mm39_anno.promoters_2kb.primary.bed \
  ${ANNOT_DIR}/mm39_anno.genes.primary.bed \
  ${ANNOT_DIR}/mm39_anno.exons.primary.bed \
  ${ANNOT_DIR}/mm39_anno.introns.primary.bed \
  ${ANNOT_DIR}/mm39_anno.intergenic.primary.bed \
  ${ANNOT_DIR}/mm39_anno.cpg_islands.primary.bed \
  ${ANNOT_DIR}/mm39_anno.cpg_shores.primary.bed \
  ${ANNOT_DIR}/mm39_anno.cpg_shelves.primary.bed \
  ${ANNOT_DIR}/mm39_anno.cpg_opensea.primary.bed \
  ${ANNOT_DIR}/mm39_anno.repeatmasker*.primary.bed \
  ${ANNOT_DIR}/mm39_anno.cCRE_*.primary.bed
do
  [ -e "$f" ] || continue
  echo "sorting $(basename "$f")"
  tmp="${f}.sorted.tmp"
  if bedtools sort -faidx "$CHROMS" -i "$f" > "$tmp"; then
    mv "$tmp" "$f"
  else
    echo "ERROR: bedtools sort failed for $f" >&2
    rm -f "$tmp"
    exit 1
  fi
done

shopt -u nullglob

# --- Minimal sanity checks (print line counts) ---
for f in genes.primary.bed tss.primary.bed promoters_2kb.primary.bed \
         cpg_islands.primary.bed cpg_shores.primary.bed cpg_shelves.primary.bed cpg_opensea.primary.bed \
         exons.primary.bed introns.primary.bed intergenic.primary.bed \
         repeatmaskerDNA.primary.bed repeatmaskerLINE.primary.bed repeatmaskerLTR.primary.bed \
         repeatmaskerSINE.primary.bed repeatmaskerSatellite.primary.bed repeatmaskerSimple_repeat.primary.bed \
         repeatmaskerLow_complexity.primary.bed repeatmaskerNoRepeat.primary.bed; do
  echo -e "$(wc -l < ${OUT}.${f})\t${OUT}.${f}"
done
