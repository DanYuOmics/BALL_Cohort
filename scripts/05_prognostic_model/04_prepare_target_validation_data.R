#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_prepare_target_validation_data.R
# Description: Preprocess and clean the independent external validation cohort 
#              (TARGET-ALL Phase 2):
#              1. Filter B-ALL clinical records and resolve duplicate sample IDs.
#              2. Quantify subclonal mutation counts per patient from PyClone outputs.
#              3. Standardize clinical parameters (NCI Risk Score, MRD Status, Subtype).
#              4. Serialize the ready-to-use matrix for survival validation.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# 1. Environment Configurations & Input Verification
# ------------------------------------------------------------------------------
message("Step 1: Setting up environment and loading raw TARGET datasets...")

target_snv_file  <- Sys.getenv("TARGET_SNV", unset = "./TARGET_ALL_Phase2_Anchored_SNV.csv")
pyclone_out_file <- Sys.getenv("TARGET_PYCLONE", unset = "./TARGET_Pyclone_Output.tsv")
target_clin_file <- Sys.getenv("TARGET_CLIN", unset = "./df_clinical_base.csv")
output_dir       <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(target_snv_file) || !file.exists(pyclone_out_file) || !file.exists(target_clin_file)) {
  stop("⚠️ Error: Missing one or more raw TARGET input files. Please verify file paths.")
}

target_snv_raw   <- fread(target_snv_file)
pyclone_real_out <- fread(pyclone_out_file)
df_clinical_base <- fread(target_clin_file)

# ------------------------------------------------------------------------------
# 2. Filter Clinical Records and Standardize Patient Identifiers
# ------------------------------------------------------------------------------
message("Step 2: Resolving unique MATCH_IDs and filtering sample types...")

df_clinical_filtered <- as_tibble(df_clinical_base) %>%
  mutate(sample_type = substr(SAMPLE_ID, 18, 19)) %>%
  arrange(MATCH_ID, sample_type) %>%
  group_by(MATCH_ID) %>%
  filter(row_number() == 1) %>% # Deduplicate and retain primary tumor barcode
  ungroup() %>%
  distinct(MATCH_ID, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 3. Extract Subclonal Mutation Burden from PyClone Outputs
# ------------------------------------------------------------------------------
message("Step 3: Calculating patient-level subclonal mutation counts from PyClone results...")

target_subclonal_counts <- pyclone_real_out %>%
  mutate(MATCH_ID = substr(mutation_id, 1, 16)) %>%
  group_by(MATCH_ID) %>%
  summarise(Subclonal_Mut_Count = sum(cluster_id >= 1), .groups = "drop")

snv_for_pyclone <- target_snv_raw %>%
  mutate(MATCH_ID = substr(Tumor_Sample_Barcode, 1, 16)) %>%
  semi_join(target_subclonal_counts, by = "MATCH_ID")

# ------------------------------------------------------------------------------
# 4. Harmonize Clinical Variables & Survival Outcomes
# ------------------------------------------------------------------------------
message("Step 4: Harmonizing survival outcomes, NCI score, MRD binary status, and molecular subtype...")

target_matrix_clean <- df_clinical_filtered %>%
  inner_join(target_subclonal_counts, by = "MATCH_ID") %>%
  inner_join(snv_for_pyclone %>% distinct(MATCH_ID), by = "MATCH_ID") %>%
  mutate(
    Surv_Months            = as.numeric(DAYS_TO_EVENT) / 30.4375,
    Event_Fixed            = if_else(is.na(FIRST_EVENT) | FIRST_EVENT == "None", 0, 1),
    NCI_Score              = if_else(AGE >= 10 | WBC >= 50, 1, 0),
    MRD_Binary             = if_else(is.na(MRD_PERCENT_DAY_29) | MRD_PERCENT_DAY_29 <= 0.01, 0, 1),
    Molecular_Subtype_Fact = as.factor(replace_na(MOLECULAR_SUBTYPE, "Not_Determined"))
  ) %>%
  filter(!is.na(Surv_Months) & Surv_Months > 0, !is.na(Event_Fixed))

message(sprintf("✅ TARGET Validation Matrix successfully assembled!"))
message(sprintf("   - Total Samples (N) = %d", nrow(target_matrix_clean)))
message(sprintf("   - Total Events      = %d", sum(target_matrix_clean$Event_Fixed)))

# ------------------------------------------------------------------------------
# 5. Export Clean Validation Matrix
# ------------------------------------------------------------------------------
output_csv <- file.path(output_dir, "TARGET_Validation_Cohort_Clean.csv")
write.csv(target_matrix_clean, output_csv, row.names = FALSE)

message(sprintf("🎉 Clean validation cohort table serialized to: %s", output_csv))
