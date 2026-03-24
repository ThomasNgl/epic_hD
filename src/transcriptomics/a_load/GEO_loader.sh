# ============================================================
# PART 2 — GET THE SRR ACCESSION LIST
# ============================================================
# A GSE accession (GEO) does not work directly with prefetch.
# You need the SRR accession numbers (SRA), which are linked
# to the GSE but live in a separate database.

# OPTION A — Web (recommended for first-time exploration)
# 1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81216
# 2. click "SRA Run Selector" link
# 3. Filter samples as needed (organism, library type, etc.)
# 4. Click "Accession List" → downloads SRR_Acc_List.txt
#    (one SRR accession per line)

# OPTION B — Command line using Entrez Direct tools (esearch/efetch)
# Requires: sudo apt-get install ncbi-entrez-direct
# Use the BioProject accession (PRJNA...) found on the GEO page
esearch -db sra -query PRJNA320896 | \
    efetch -format runinfo | \       # fetches full metadata as CSV
    cut -d',' -f1 | \                # extracts first column (SRR accessions)
    grep SRR \                       # keeps only lines starting with SRR
    > SRR_Acc_List.txt               # saves to file


# ============================================================
# PART 3 — DOWNLOAD
# ============================================================

# OPTION A — Download a single run
prefetch SRR_accession_number

# OPTION B — Download all runs from a list (recommended for a full dataset)
# prefetch will resume automatically if interrupted — safe to re-run
prefetch --option-file SRR_Acc_List.txt


# ============================================================
# PART 4 — CONVERT .sra → .fastq
# ============================================================
# Run fasterq-dump from the SAME directory where you ran prefetch
# fasterq-dump is faster than fastq-dump (multithreaded)
# Note: fasterq-dump has no --gzip option — compress manually after

# --- Understanding split options ---
#
# --split-files  → always produces exactly 2 files (R1, R2)
#                  unmatched reads are DISCARDED
#                  output: SRR_1.fastq  SRR_2.fastq
#
# --split-3      → DEFAULT in fasterq-dump
#                  produces 2 files for paired-end (R1, R2)
#                  + optional 3rd file for unmatched/singleton reads
#                  output: SRR_1.fastq  SRR_2.fastq  SRR.fastq (if singletons exist)
#
# For single-end data (e.g. small RNA-seq like GSE81216):
#   both options produce a single file — no practical difference

# Single-end (e.g. small RNA-seq) — produces one file per run
fasterq-dump SRR_accession_number
gzip *.fastq

# Paired-end — produces two (or three) files per run
fasterq-dump --split-3 SRR_accession_number
gzip *.fastq

# Loop over all runs in your list
while read SRR; do
    fasterq-dump --split-3 "$SRR"
done < SRR_Acc_List.txt
gzip *.fastq

