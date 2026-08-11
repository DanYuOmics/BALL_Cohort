#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_subclone_diversity_survival.R
# Description: Calculate Shannon and Simpson Diversity Indices,
#              fit Cox proportional hazards models for Event-Free Survival (EFS),
#              determine optimal cutoff points via survminer, plot Kaplan-Meier
#              curves, and export summary statistics to Excel.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(survival)
  library(survminer)
  library(openxlsx)
  library(broom)
})

# ------------------------------------------------------------------------------
# 1. Environment Setup & Directory Configurations
# ------------------------------------------------------------------------------
message("Step 1: Setting up environment and loading clinical metadata...")

clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
pyclone_txt   <- Sys.getenv("PYCLONE_TXT_OUT", unset = "./article/Fig2_Upgraded_Clonal_Architecture_TXT.tsv")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load clinical metadata
if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical metadata file at: ", clinical_file)
}

clinical_data_raw_updated <- read_excel(clinical_file) %>%
  mutate(
    event_binary = case_when(
      Event %in% c(1, "1", "Relapse/Died", "RELAPSED") ~ 1,
      Event %in% c(0, "0", "Event-free", "NON-RELAPSED") ~ 0,
      TRUE ~ NA_real_
    )
  )

message(sprintf("Loaded clinical metadata for %d patients.", nrow(clinical_data_raw_updated)))

# ------------------------------------------------------------------------------
# 2. Load PyClone Outputs & Calculate Diversity Metrics
# ------------------------------------------------------------------------------
message("Step 2: Loading PyClone outputs and computing Shannon & Simpson indices...")

if (!file.exists(pyclone_txt)) {
  stop("⚠️ Error: Missing PyClone output file at: ", pyclone_txt)
}

pyclone_results <- read.table(pyclone_txt, header = TRUE, sep = "\t")

fig2_ready_df <- pyclone_results %>%
  mutate(
    sample_id = sapply(strsplit(as.character(mutation_id), "_"), `[`, 1)
  )

# Calculate Shannon Diversity Index per patient
diversity_df <- fig2_ready_df %>%
  group_by(sample_id) %>%
  summarise(shannon = -sum(cellular_fraction * log(cellular_fraction + 1e-10)), .groups = "drop") %>%
  left_join(
    select(clinical_data_raw_updated, Patient_ID, EFS_Month, event_binary),
    by = c("sample_id" = "Patient_ID")
  ) %>%
  filter(event_binary %in% c(0, 1)) %>%
  rename(event_binary_final = event_binary, EFS_time = EFS_Month) %>%
  mutate(EFS_time = as.numeric(EFS_time)) %>%
  filter(!is.na(EFS_time))

# Calculate Simpson Diversity Index per patient
simpson_df <- fig2_ready_df %>%
  group_by(sample_id) %>%
  summarise(simpson = 1 - sum(cellular_fraction^2), .groups = "drop") %>%
  left_join(
    select(clinical_data_raw_updated, Patient_ID, EFS_Month, event_binary),
    by = c("sample_id" = "Patient_ID")
  ) %>%
  filter(event_binary %in% c(0, 1)) %>%
  rename(event_binary_final = event_binary, EFS_time = EFS_Month) %>%
  mutate(EFS_time = as.numeric(EFS_time)) %>%
  filter(!is.na(EFS_time))

# ------------------------------------------------------------------------------
# 3. Fit Continuous Cox Proportional Hazards Models
# ------------------------------------------------------------------------------
message("Step 3: Fitting Cox proportional hazards regression models...")

cox_shannon <- coxph(Surv(EFS_time, event_binary_final) ~ shannon, data = diversity_df)
cox_simpson <- coxph(Surv(EFS_time, event_binary_final) ~ simpson, data = simpson_df)

message("--- Cox Proportional Hazards Model (Shannon Diversity) ---")
print(summary(cox_shannon))

message("--- Cox Proportional Hazards Model (Simpson Diversity) ---")
print(summary(cox_simpson))

# ------------------------------------------------------------------------------
# 4. Optimal Cutoff Determination & Kaplan-Meier Survival Analysis
# ------------------------------------------------------------------------------
message("Step 4: Determining optimal cutoff and plotting Kaplan-Meier survival curves...")

res.cut <- surv_cutpoint(
  diversity_df,
  time      = "EFS_time",
  event     = "event_binary_final",
  variables = "shannon"
)

message("--- Optimal Cutoff Summary ---")
print(summary(res.cut))

# Stratify cohort into High vs. Low Shannon Diversity groups
diversity_df$shannon_group <- ifelse(
  diversity_df$shannon > res.cut$cutpoint$cutpoint,
  "High", "Low"
)

fit <- survfit(Surv(EFS_time, event_binary_final) ~ shannon_group, data = diversity_df)

# Render KM Curve Plot
p_km <- ggsurvplot(
  fit,
  data          = diversity_df,
  pval          = TRUE,
  risk.table    = TRUE,
  palette       = c("darkgreen", "red"),
  legend.title  = "Shannon Diversity",
  legend.labs   = c("Low", "High"),
  xlab          = "Months",
  ylab          = "Event-Free Survival Probability",
  title         = "KM Curve by Shannon Diversity (Optimal Cutoff)"
)

# Export KM Plot PDF
km_pdf_path <- file.path(output_dir, "Figure3_Shannon_Diversity_KM_Curve.pdf")
pdf(km_pdf_path, width = 7, height = 7, onefile = FALSE)
print(p_km)
dev.off()

message("Saved KM Survival Plot to: ", km_pdf_path)

# ------------------------------------------------------------------------------
# 5. Export Analytical Results and Source Data to Excel
# ------------------------------------------------------------------------------
message("Step 5: Exporting source data and Cox summary to Excel...")

shannon_source <- diversity_df %>%
  select(sample_id, shannon, shannon_group, EFS_time, event_binary_final)

cox_shannon_summary <- broom::tidy(cox_shannon, conf.int = TRUE, exponentiate = TRUE) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

excel_output_path <- file.path(output_dir, "Shannon_Diversity_Survival.xlsx")

write.xlsx(
  list(
    Shannon_Source = shannon_source,
    Cox_Summary    = cox_shannon_summary
  ),
  file     = excel_output_path,
  rowNames = FALSE
)

message("🎉 Analytical workflow completed successfully!")
message("  - Excel Summary: ", excel_output_path)
