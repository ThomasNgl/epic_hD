library(AnnotationHub)
library(ensembldb)

# ===== GENE ID CONVERSION =====

convertGeneNames <- function(ids, gene_map, strip_version = TRUE) {
  
  # Accept either a file path or a loaded gene_map object
  if (is.character(gene_map) && length(gene_map) == 1 && file.exists(gene_map)) {
    gene_map <- readRDS(gene_map)
  } else if (!is.data.frame(gene_map)) {
    stop("gene_map must be either a path to an RDS file or a data.frame")
  }
  
  # Optionally strip version suffixes
  ids_lookup <- if (strip_version) sub("\\.\\d+$", "", ids) else ids
  map_ids    <- if (strip_version) sub("\\.\\d+$", "", gene_map$gene_id) else gene_map$gene_id
  
  # Match Ensembl IDs to symbols
  symbols <- gene_map$gene_symbol[ match(ids_lookup, map_ids) ]
  
  # Report unmatched
  n_na <- sum(is.na(symbols))
  if (n_na > 0) message(n_na, " IDs had no symbol match — original ID kept as fallback")
  
  # Return data.frame with both
  data.frame(
    gene_id     = ids,
    gene_symbol = ifelse(is.na(symbols), ids, symbols),
    stringsAsFactors = FALSE
  )
}

# ===== FETCH EnsDb =====

# Load EnsDb for Mus musculus vM37 
ah <- AnnotationHub(
    cache = "/mnt/auxiliary/.cache/R/AnnotationHub",
    ask = FALSE
)

query_results <- query(ah, c("EnsDb", "Mus musculus", "110"))

# Extract AH ID for GRCm39
ah_id <- rownames(mcols(query_results))[mcols(query_results)$genome == "GRCm39"]

# Safety check — should be exactly one
stopifnot(length(ah_id) == 1)

# Load it
edb <- ah[[ah_id]] # Ensembl 110 / GENCODE vM37

# Get extra annotation
gene_info <- genes(edb, columns = c(
  "gene_id",
  "gene_name", 
  "gene_biotype",
  "seq_name",        # chromosome
  "gene_seq_start",
  "gene_seq_end"
))
gene_info_df <- as.data.frame(gene_info)

# Strip versions from your SE rownames for matching
gene_ids_noversion <- sub("\\.\\d+$", "", rownames(se_gene))

# Join by gene_id
idx <- match(gene_ids_noversion, gene_info_df$gene_id)

# Add columns to rowData
rowData(se_gene)$gene_biotype  <- gene_info_df$gene_biotype[idx]
rowData(se_gene)$chromosome    <- gene_info_df$seq_name[idx]
rowData(se_gene)$gene_start    <- gene_info_df$gene_seq_start[idx]
rowData(se_gene)$gene_end      <- gene_info_df$gene_seq_end[idx]
  cache = "/mnt/auxiliary/.cache/R/AnnotationHub"

