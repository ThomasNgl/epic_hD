# ── Dependencies ───────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(scales)

# ── Shared colour palette ──────────────────────────────────────────────────
BIOTYPE_COLORS <- c(
  protein_coding = "#3266ad",
  lncRNA         = "#1d9e75",
  pseudogene     = "#888780",
  small_ncRNA    = "#ef9f27",
  IG_TR_gene     = "#d85a30",
  TEC            = "#b4b2a9",
  rRNA_tRNA      = "#e24b4a",
  unknown        = "#d4537e",
  other          = "#534ab7"
)

simplify_biotype <- function(bt) {
  dplyr::case_when(
    bt == "protein_coding"                                         ~ "protein_coding",
    bt == "lncRNA"                                                 ~ "lncRNA",
    bt %in% c("processed_pseudogene", "unprocessed_pseudogene",
               "transcribed_processed_pseudogene",
               "transcribed_unprocessed_pseudogene",
               "transcribed_unitary_pseudogene",
               "unitary_pseudogene",
               "translated_unprocessed_pseudogene",
               "IG_V_pseudogene", "TR_V_pseudogene",
               "TR_J_pseudogene", "IG_D_pseudogene")              ~ "pseudogene",
    bt %in% c("miRNA", "snRNA", "snoRNA",
               "scaRNA", "misc_RNA", "ribozyme")                  ~ "small_ncRNA",
    bt %in% c("IG_V_gene", "IG_D_gene", "IG_J_gene", "IG_C_gene",
               "TR_V_gene", "TR_J_gene", "TR_C_gene", "TR_D_gene") ~ "IG_TR_gene",
    bt == "TEC"                                                    ~ "TEC",
    bt %in% c("rRNA", "Mt_tRNA", "Mt_rRNA")                      ~ "rRNA_tRNA",
    is.na(bt)                                                      ~ "unknown",
    TRUE                                                           ~ "other"
  )
}

# ── Main function ──────────────────────────────────────────────────────────
plot_biotype <- function(se, gene_map, save_dir, label = NULL) {

  # -- Load se if path
  if (is.character(se) && length(se) == 1 && file.exists(se)) {
    se <- readRDS(se)
  } else if (!is(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object or a path to an RDS file")
  }

  # -- Load gene_map if path
  if (is.character(gene_map) && length(gene_map) == 1 && file.exists(gene_map)) {
    gene_map <- readRDS(gene_map)
  } else if (!is.data.frame(gene_map)) {
    stop("gene_map must be a data.frame or a path to an RDS file")
  }

  # -- Create save directory
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  # -- Attach biotypes
  ids     <- sub("\\.\\d+$", "", rownames(se))
  map_ids <- sub("\\.\\d+$", "", gene_map$gene_id)
  bt      <- gene_map$gene_biotype[match(ids, map_ids)]
  n_na    <- sum(is.na(bt))
  if (n_na > 0) message(n_na, " genes with no biotype match — labelled 'unknown'")
  rowData(se)$gene_biotype <- bt

  # -- Count biotypes
  n_features  <- nrow(se)
  plot_label  <- if (!is.null(label)) label else paste0("n=", scales::comma(n_features))
  bt_simple   <- simplify_biotype(rowData(se)$gene_biotype)
  bt_df       <- as.data.frame(table(biotype_simple = bt_simple, useNA = "ifany"))
  bt_df$pct   <- round(100 * bt_df$Freq / sum(bt_df$Freq), 1)
  bt_df$biotype_simple <- factor(bt_df$biotype_simple, levels = names(BIOTYPE_COLORS))

  # -- Plot counts
  p_count <- ggplot(bt_df, aes(x = biotype_simple, y = Freq, fill = biotype_simple)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = scales::comma(Freq)),
              vjust = -0.3, size = 2.6, color = "grey30") +
    scale_fill_manual(values = BIOTYPE_COLORS, na.value = "#d4537e") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                       labels = scales::comma) +
    labs(title    = paste0("Biotype distribution — counts (", plot_label, ")"),
         x = NULL, y = "Number of genes") +
    theme_minimal(base_size = 12) +
    theme(legend.position    = "none",
          axis.text.x        = element_text(angle = 35, hjust = 1),
          panel.grid.major.x = element_blank(),
          plot.title         = element_text(face = "bold"))

  # -- Plot percentages
  p_pct <- ggplot(bt_df, aes(x = biotype_simple, y = pct, fill = biotype_simple)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = paste0(pct, "%")),
              vjust = -0.3, size = 2.6, color = "grey30") +
    scale_fill_manual(values = BIOTYPE_COLORS, na.value = "#d4537e") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                       labels = function(x) paste0(x, "%")) +
    labs(title    = paste0("Biotype distribution — percentage (", plot_label, ")"),
         x = NULL, y = "% of genes") +
    theme_minimal(base_size = 12) +
    theme(legend.position    = "none",
          axis.text.x        = element_text(angle = 35, hjust = 1),
          panel.grid.major.x = element_blank(),
          plot.title         = element_text(face = "bold"))

  # -- Save
  ggsave(file.path(save_dir, "biotype_count.pdf"),  p_count, width = 10, height = 6)
  ggsave(file.path(save_dir, "biotype_pct.pdf"),    p_pct,   width = 10, height = 6)
  write.csv(bt_df, file.path(save_dir, "biotype_summary.csv"), row.names = FALSE)

  message("Saved to: ", save_dir)
  invisible(list(count = p_count, pct = p_pct, table = bt_df, se = se))
}