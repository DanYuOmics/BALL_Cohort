#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_subclone_evolution_figures.R
# Description: Import normalized PyClone-VI output matrix, 
#              de-multiplex patient barcodes and gene symbols, integrate 
#              clinical metadata, and structure data for Figure 3 visualization.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(dplyr)
  library(readr)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/B_ALL_105_Clinical_Data_Final.xlsx")
pyclone_file  <- Sys.getenv("PYCLONE_TXT_OUT", unset = "./article/Fig2_Upgraded_Clonal_Architecture_TXT.tsv")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Load Clinical Metadata
# ------------------------------------------------------------------------------
message("Step 1: Loading clinical metadata...")

if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical metadata file at: ", clinical_file)
}

clinical_data <- read_xlsx(clinical_file) %>%
  mutate(Event_Label = if_else(Event == 1, "Relapsed/Deceased", "Event-free"))

message(sprintf("Loaded clinical metadata for %d patients.", nrow(clinical_data)))

# ------------------------------------------------------------------------------
# 3. Re-import Normalized Text Matrix and De-multiplex Metadata
# ------------------------------------------------------------------------------
message("Step 2: Ingesting PyClone-VI output matrix and de-multiplexing mutation IDs...")

if (!file.exists(pyclone_file)) {
  stop("⚠️ Error: Missing PyClone output file at: ", pyclone_file)
}

pyclone_results <- read_tsv(pyclone_file, show_col_types = FALSE)

# De-multiplex barcodes and gene symbols from composite mutation_id
fig2_ready_df <- pyclone_results %>%
  mutate(
    # De-multiplex patient barcodes from the mutation_id string prefix
    real_sample_id = sapply(strsplit(as.character(mutation_id), "_"), `[`, 1),
    # Re-extract standard gene symbols for targeted recurrent driver visualization
    hugo_symbol    = sapply(strsplit(as.character(mutation_id), "_"), `[`, 2)
  ) %>%
  dplyr::select(mutation_id, sample_id = real_sample_id, hugo_symbol, cluster_id, cellular_fraction) %>%
  # Filter to ensure alignment with clinical cohort
  filter(sample_id %in% clinical_data$Patient_ID) %>%
  left_join(
    clinical_data %>% dplyr::select(Patient_ID, Event_Label),
    by = c("sample_id" = "Patient_ID")
  )

# ------------------------------------------------------------------------------
# 4. Preview and Save Standardized Multi-Omics Dataframe
# ------------------------------------------------------------------------------
message("Step 3: Previewing final standardized multi-omics dataset...")
print(head(fig2_ready_df))

output_rds <- file.path(output_dir, "fig3_subclone_architecture_ready.rds")
saveRDS(fig2_ready_df, output_rds)

message("🎉 Successfully parsed and saved subclonal architecture dataframe to: ", output_rds)
