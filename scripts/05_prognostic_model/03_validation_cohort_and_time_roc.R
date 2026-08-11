#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_validation_cohort_and_time_roc.R
# Description: Clean external validation cohort (TARGET-ALL Phase 2), evaluate 
#              evolution-enhanced 3-tier KM stratification, compute 24-month timeROC 
#              AUCs strictly via dynamic calculation for both Discovery and Validation 
#              cohorts, and generate predictive performance bar plots.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Environment Configurations & Path Setup
# ------------------------------------------------------------------------------
message("Step 1: Setting up environment and directory paths...")

discovery_csv_file <- Sys.getenv("DISCOVERY_CSV", unset = "./results/Discovery_Cohort_Model_SourceData.csv")
target_snv_file    <- Sys.getenv("TARGET_SNV", unset = "./TARGET_ALL_Phase2_Anchored_SNV.csv")
pyclone_out_file   <- Sys.getenv("TARGET_PYCLONE", unset = "./TARGET_Pyclone_Output.tsv")
target_clin_file   <- Sys.getenv("TARGET_CLIN", unset = "./df_clinical_base.csv")
output_dir         <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Dynamic AUC Calculation for Discovery Cohort
# ------------------------------------------------------------------------------
message("Step 2: Dynamically calculating 24-month timeROC AUCs for Discovery Cohort...")

if (!file.exists(discovery_csv_file)) {
  stop("⚠️ Error: Missing Discovery Cohort model source data at: ", discovery_csv_file, 
       "\nPlease run '02_discovery_cohort_survival.R' first!")
}

discovery_data <- fread(discovery_csv_file)

# Mathematically calculate 24-month timeROC for Discovery Cohort
roc_disc_base <- timeROC(
  T = discovery_data$EFS_Month,
  delta = discovery_data$Event,
  marker = discovery_data$LP_Modern_Base,
  cause = 1, times = 24, iid = TRUE
)

roc_disc_evo <- timeROC(
  T = discovery_data$EFS_Month,
  delta = discovery_data$Event,
  marker = discovery_data$LP_Modern_Evo,
  cause = 1, times = 24, iid = TRUE
)

disc_auc_base <- as.numeric(roc_disc_base$AUC[2])
disc_auc_evo  <- as.numeric(roc_disc_evo$AUC[2])

message(sprintf("✅ Discovery Cohort Calculated AUCs at 24 Months:"))
message(sprintf("   - Clinical Baseline: %.4f", disc_auc_base))
message(sprintf("   - Evolution-Enhanced: %.4f", disc_auc_evo))

# ------------------------------------------------------------------------------
# 3. Clean TARGET Validation Cohort & Calculate Dynamic AUCs
# ------------------------------------------------------------------------------
message("Step 3: Cleaning TARGET Validation Cohort and calculating dynamic AUCs...")

if (!file.exists(target_snv_file) || !file.exists(pyclone_out_file) || !file.exists(target_clin_file)) {
  stop("⚠️ Error: Missing required TARGET raw data files. Cannot calculate validation AUCs!")
}

target_snv_raw   <- fread(target_snv_file)
pyclone_real_out <- fread(pyclone_out_file)
df_clinical_base <- fread(target_clin_file)

# Clean and align TARGET clinical metadata
df_clinical_filtered <- as_tibble(df_clinical_base) %>%
  mutate(sample_type = substr(SAMPLE_ID, 18, 19)) %>%
  arrange(MATCH_ID, sample_type) %>%
  group_by(MATCH_ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  distinct(MATCH_ID, .keep_all = TRUE)

target_subclonal_counts <- pyclone_real_out %>%
  mutate(MATCH_ID = substr(mutation_id, 1, 16)) %>%
  group_by(MATCH_ID) %>%
  summarise(Subclonal_Mut_Count = sum(cluster_id >= 1), .groups = "drop")

snv_for_pyclone <- target_snv_raw %>%
  mutate(MATCH_ID = substr(Tumor_Sample_Barcode, 1, 16)) %>%
  semi_join(target_subclonal_counts, by = "MATCH_ID")

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

# Fit validation Cox models
fit_val_base <- coxph(Surv(Surv_Months, Event_Fixed) ~ NCI_Score + MRD_Binary + Molecular_Subtype_Fact, data = target_matrix_clean)
target_matrix_clean$LP_Base <- predict(fit_val_base, type = "lp")

fit_val_evo <- coxph(Surv(Surv_Months, Event_Fixed) ~ NCI_Score + MRD_Binary + Molecular_Subtype_Fact + Subclonal_Mut_Count, data = target_matrix_clean)
target_matrix_clean$LP_Evo <- predict(fit_val_evo, type = "lp")

# Render Validation Cohort 3-Tier KM Curve
cuts_val <- quantile(target_matrix_clean$LP_Evo, probs = c(0.33, 0.66), na.rm = TRUE)

target_matrix_clean <- target_matrix_clean %>%
  mutate(
    Tier_Evo_3way = case_when(
      LP_Evo <= cuts_val[1] ~ "Low Risk",
      LP_Evo > cuts_val[2]  ~ "High Risk",
      TRUE                  ~ "Intermediate Risk"
    ),
    Tier_Evo_3way = factor(Tier_Evo_3way, levels = c("Low Risk", "Intermediate Risk", "High Risk"))
  )

fit_3way_val <- survfit(Surv(Surv_Months, Event_Fixed) ~ Tier_Evo_3way, data = target_matrix_clean)

evo_counts <- target_matrix_clean %>%
  count(Tier_Evo_3way) %>%
  mutate(label = paste0(Tier_Evo_3way, " (n=", n, ")"))

p_3way_val <- ggsurvplot(
  fit_3way_val, data = target_matrix_clean,
  palette     = c("#4682B4", "#6E8B3D", "#CD5C5C"),
  pval        = TRUE, pval.method = TRUE,
  risk.table  = FALSE,
  ggtheme     = theme_classic(),
  xlab        = "Follow-up Time (Months)",
  ylab        = "Event-Free Survival Probability",
  title       = "Validation Cohort (TARGET-ALL-P2)\nNCI + MRD + Subtype + Subclonal Mutation Count",
  legend.title = "",
  legend.labs  = evo_counts$label
)

pdf_val_out <- file.path(output_dir, "Validation_Cohort_3Tier_KM_Curve.pdf")
pdf(pdf_val_out, width = 6, height = 5.5, onefile = FALSE)
print(p_3way_val)
dev.off()

# Mathematically calculate 24-month timeROC for Validation Cohort
roc_val_base <- timeROC(
  T = target_matrix_clean$Surv_Months,
  delta = target_matrix_clean$Event_Fixed,
  marker = target_matrix_clean$LP_Base,
  cause = 1, times = 24, iid = TRUE
)

roc_val_evo <- timeROC(
  T = target_matrix_clean$Surv_Months,
  delta = target_matrix_clean$Event_Fixed,
  marker = target_matrix_clean$LP_Evo,
  cause = 1, times = 24, iid = TRUE
)

val_auc_base <- as.numeric(roc_val_base$AUC[2])
val_auc_evo  <- as.numeric(roc_val_evo$AUC[2])

message(sprintf("✅ Validation Cohort Calculated AUCs at 24 Months:"))
message(sprintf("   - Clinical Baseline: %.4f", val_auc_base))
message(sprintf("   - Evolution-Enhanced: %.4f", val_auc_evo))

# ------------------------------------------------------------------------------
# 4. Render AUC Bar Plot strictly using Dynamically Calculated Values
# ------------------------------------------------------------------------------
message("Step 4: Assembling AUC summary table using dynamically calculated metrics...")

auc_summary_df <- data.frame(
  Cohort = c("Discovery Cohort", "Discovery Cohort", 
             "Validation Cohort (TARGET-ALL-P2)", "Validation Cohort (TARGET-ALL-P2)"),
  Model_Type = c("NCI + MRD + Subtype", 
                 "NCI + MRD + Subtype + Subclonal Mutation Count", 
                 "NCI + MRD + Subtype", 
                 "NCI + MRD + Subtype + Subclonal Mutation Count"),
  AUC_24M = c(disc_auc_base, disc_auc_evo, val_auc_base, val_auc_evo)
) %>%
  mutate(
    Cohort = factor(Cohort, levels = c("Discovery Cohort", "Validation Cohort (TARGET-ALL-P2)")),
    Model_Type = factor(Model_Type, levels = c("NCI + MRD + Subtype", "NCI + MRD + Subtype + Subclonal Mutation Count"))
  )

p_auc_comparison <- ggplot(auc_summary_df, aes(x = Cohort, y = AUC_24M, fill = Model_Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.35, color = "black", linewidth = 0.3, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", AUC_24M)), position = position_dodge(0.7), vjust = -0.6, fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("NCI + MRD + Subtype" = "#6A8CA0", "NCI + MRD + Subtype + Subclonal Mutation Count" = "#A05252")) +
  theme_classic() +
  labs(
    title = "Predictive Performance Comparison",
    x = "", y = "Area Under the Curve (AUC)", fill = "Model"
  ) +
  coord_cartesian(ylim = c(0.6, 1.0)) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.text        = element_text(color = "black", size = 10),
    axis.title.y     = element_text(face = "bold", size = 10),
    legend.position  = "right",
    legend.direction = "vertical",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm")
  )

pdf_auc_out <- file.path(output_dir, "Predictive_Performance_AUC_Comparison.pdf")
ggsave(pdf_auc_out, plot = p_auc_comparison, width = 7.5, height = 4.5)

write.csv(auc_summary_df, file.path(output_dir, "Predictive_AUC_Summary_SourceData.csv"), row.names = FALSE)

message("🎉 Dynamically calculated ROC/AUC performance evaluation complete!")
message("  - PDF Bar Plot Output: ", pdf_auc_out)
