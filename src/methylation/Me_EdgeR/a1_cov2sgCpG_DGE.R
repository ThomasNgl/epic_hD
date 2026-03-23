library(data.table)
library(dplyr)
library(readr)
library(edgeR)
library(tibble)

attach_pheno_to_samples <- function(me_dge, pheno) {
  stopifnot("Sample" %in% names(pheno))
  
  samp_raw <- rownames(me_dge$samples)
  if (is.null(samp_raw)) stop("me_dge$samples has no rownames; see Option B.")
  
  key <- tibble(
    .idx = seq_along(samp_raw),
    .sample_raw = samp_raw,
    Sample = sub("-(Me|Un)$", "", samp_raw)   # strip suffix only at end
  )
  
  # join pheno onto the samples table
  samp2 <- as_tibble(me_dge$samples, rownames = ".sample_raw") %>%
    left_join(key %>% select(.sample_raw, Sample), by = ".sample_raw") %>%
    left_join(pheno, by = "Sample") %>%
    column_to_rownames(".sample_raw")
  
  me_dge$samples <- samp2
  me_dge
}

# --- Function Definition ---

#' Reads Bismark coverage files and creates an edgeR DGEList object.
#'
#' This function performs the initial data loading and filtering steps for
#' differential methylation analysis using edgeR, converting Bismark
#' coverage files into a DGEList object.
#'
#' @param meth_extracted_dir Character. Directory containing the Bismark .cov.gz files.
#' @param pheno_dir Character. Path to the CSV file containing sample metadata (phenotype).
#' @param output_dir Character. Directory where the final DGEList object will be saved.
#' @param output_name Character. The filename for the saved DGEList object (e.g., "me_dge_initial.rds").
#' @return The edgeR DGEList object (invisibly).
#' @export
cov.gz2me_dge <- function(meth_extracted_dir, pheno_dir, output_dir, output_name, missing_samples =NULL
) {
  
  # 1. Load and Prepare Phenotype Data
  # ----------------------------------
  pheno <- read_csv(pheno_dir, show_col_types = FALSE)
  
  # Define samples to be excluded
  if (!is.null(missing_samples) && length(missing_samples) > 0) {
    pheno <- pheno %>%
      filter(!Sample %in% missing_samples)
  }
  # 2. Create the Targets Data Frame
  # --------------------------------
  # Create a table of sample names and corresponding file paths
  files <- list.files(meth_extracted_dir,
                      pattern = "\\.bismark\\.cov\\.gz$",
                      recursive = TRUE,
                      full.names = TRUE)

  targets <- tibble(Sample = pheno$Sample) %>%
    rowwise() %>%
    mutate(File = {
      hit <- files[grepl(paste0("/", Sample, "/"), files)]        # matches folder name
      if (length(hit) != 1) hit <- files[grepl(Sample, files)]   # fallback: anywhere in path
      if (length(hit) == 0) stop("No .bismark.cov.gz found for Sample: ", Sample)
      if (length(hit) > 1)  stop("Multiple .bismark.cov.gz found for Sample: ", Sample,
                                "\n", paste(hit, collapse = "\n"))
      hit
    }) %>%
    ungroup() %>%
    left_join(pheno, by = "Sample")
  
  
  
  # 3. Read Bismark Data into DGEList
  # ---------------------------------
  # readBismark2DGE converts the Bismark coverage format into an edgeR DGEList.
  # The counts are stored in $counts (methylated) and $unmethylated (unmethylated)
  me_dge <- readBismark2DGE(targets$File,
                            sample.names = targets$Sample,
                            readr = TRUE,
                            verbose = FALSE)
  
  # 4. Filter and Sort by Chromosome
  # --------------------------------
  # Define valid chromosomes (mm39 primary assembly)
  ChrNames <- paste0("chr", c(1:19, "X", "Y", "M"))
  
  # Filter to keep only CpG sites on primary chromosomes
  me_dge <- me_dge[me_dge$genes$Chr %in% ChrNames, ]
  
  # Sort the CpG sites in chromosomal order for better data structure
  me_dge$genes$Chr <- factor(me_dge$genes$Chr, levels = ChrNames)
  o <- order(me_dge$genes$Chr, me_dge$genes$Locus)
  me_dge <- me_dge[o, ]
  
  message(paste("Initial DGEList created with", nrow(me_dge), "CpG sites."))
  
  # attach pheno data
  me_dge <- attach_pheno_to_samples(me_dge, pheno)
  # 5. Save the DGEList Object
  # --------------------------
  saveRDS(me_dge, file = file.path(output_dir, output_name))
  message(paste("Saved DGEList to:", file.path(output_dir, output_name)))
  
  return(invisible(me_dge))
}

# --- Example Usage (Commented Out) ---

# # Define your paths
# cpg_dir <- "/mnt/thomas/matriline/results/sg_oocyte_MSUS37_DK/WGBS/04_edgeR/01_meth_extracted/"
# pheno_path <- "/mnt/thomas/matriline/results/sg_oocyte_MSUS37_DK/WGBS/metadata/metadata_oocyte_scWGBS.csv"
# out_dir <- "/mnt/thomas/matriline/results/sg_oocyte_MSUS37_DK/WGBS/04_edgeR/03_DGE"
# out_file <- "me_dge_initial.rds"
# 
# # Make sure output directory exists
# if (!dir.exists(out_dir)) {
#   dir.create(out_dir, recursive = TRUE)
# }
# 
# # Run the function
# # me_dge <- create_me_dge(cpg_dir, pheno_path, out_dir, out_file)

