set -euo pipefail  # fail on error, unset, and pipe

# --- I/O & required inputs ---
GP_DIR="${1:?Usage: $0 /path/to/gencode_primary}"; GP_DIR="$(readlink -f "$GP_DIR")"  # resolve input dir
BASE_DIR="$(dirname "$GP_DIR")"                 # parent dir
ANNOT_DIR="${BASE_DIR}/annotations"             # output dir (sibling to gencode_primary)
GENOME="${GP_DIR}/GRCm39.primary_assembly.genome.fa"                 # mm39 FASTA (UCSC chr names)
GTF="${GP_DIR}/gencode.vM37.primary_assembly.annotation.gtf"         # GENCODE GTF
OUT="${ANNOT_DIR}/mm39_anno"                    # output prefix

# What chromosome names are in the FASTA index?
head -n 5 "${GENOME}.fai" | cut -f1

# What chromosome names are in the GTF?
awk 'BEGIN{FS="\t"} $0 !~ /^#/ {print $1; exit}' "$GTF"

# Are you filtering everything away?
wc -l "${OUT}.chrom.sizes" "${OUT}.chrom.sizes.primary" \
      "${OUT}.genes.bed" "${OUT}.genes.primary.bed" \
      "${OUT}.tss.bed"   "${OUT}.tss.primary.bed"
