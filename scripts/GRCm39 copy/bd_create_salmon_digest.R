###Code to unambiguously assign transcripts to genes (suggested by O3)
###Just in case no UMI

rm(list = ls())

setwd("/mnt/thomas/matriline/epic_msus")

# this is to create a linkedTxome for tximeta to recognize your index + GTF

# this needs to be run only once for the transcriptome version you’re using.


library(tximeta)

#tximeta needs a BiocFileCache directory to access and save TxDb objects.
#set it now:
setTximetaBFC("/mnt/auxiliary/biocache")
makeLinkedTxome(
  indexDir = "/mnt/auxiliary/salmon/mm39_salmon_index",  # path to your salmon index
  source = "GENCODE",
  organism = "Mus musculus",
  release = "M37",
  genome = "GRCm39",
  fasta = "/mnt/auxiliary/salmon/mm39_raw_gencode_primary/gencode.vM37.transcripts.fa",
  gtf = "/mnt/auxiliary/salmon/mm39_raw_gencode_primary/gencode.vM37.annotation.gtf",
  write = TRUE  # Creates a JSON file linking metadata
)
