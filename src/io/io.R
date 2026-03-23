# Load required libraries
library(SummarizedExperiment)
library(openxlsx)
library(edgeR)
library(sva)
library(limma)

load_se_list <- function(base_path, inventory_path, sheet_name, suffix = NULL) {
     # Load inventory
     inventory <- read.xlsx(inventory_path, sheet = sheet_name)
     
     # Filter rows where 'se_saved' == "X"
     filtered_inventory <- subset(inventory, se_saved == "X")
     
     # Create a named list to store SummarizedExperiment objects
     se_list <- list()
     
     # Determine pattern based on suffix
     pattern <- if (is.null(suffix)) {
          "^se_.*\\.rds$"
     } else {
          paste0("^se_", suffix, ".*\\.rds$")
     }
     
     # Load each SE object with the inventory
     for (i in seq_len(nrow(filtered_inventory))) {
          # Build full directory path
          dir_path <- file.path(base_path, filtered_inventory$se_save_path[i])
          
          # List RDS files starting with "se_" but NOT "se_filtered"
          rds_files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
          rds_files <- rds_files[!grepl("se_filtered", rds_files)]
          
          # Safety check: make sure exactly one match
          if (length(rds_files) != 1) {
               warning(paste("Unexpected number of matching RDS files in", dir_path))
               next
          }
          
          # Load the SummarizedExperiment object
          se_object <- readRDS(rds_files[1])
          
          # Use the filename (without .rds) as list key
          name <- tools::file_path_sans_ext(basename(rds_files[1]))
          
          # Add to the list
          se_list[[name]] <- se_object
     }
     
     return(se_list)
}
