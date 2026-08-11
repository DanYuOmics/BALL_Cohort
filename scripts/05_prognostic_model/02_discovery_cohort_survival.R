#!/usr/bin/env Rscript
# ==============================================================================
# Script: 02_discovery_cohort_survival.R
# Description: Fit baseline clinical and evolution-enhanced multivariate Cox 
#              models in the Discovery Cohort, compute Linear Predictors (LP), 
#              stratify into 3-tier risk groups (tertiles), and plot KM curves.
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(readxl)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical metadata file at: ", clinical_file)
}

clinical_data_augmented <- read_excel(clinical_file)

# ------------------------------------------------------------------------------
# 2. Data Cleaning & Model Fitting (Modern Clinical Baseline vs Evolution-Enhanced)
# ------------------------------------------------------------------------------
message("Step 1: Filtering complete cases for Discovery Cohort risk modeling...")

our_modern_data <- clinical_data_augmented %>%
  filter(
    !is.na(EFS_Month), !is.na(Event),
    !is.na(NCI_Risk_Group), !is.na(Molecular_Subtype), !is.na(True_Molecular_MRD),
    !is.na(Subclonal_Mut_Count)
  ) %>%
  mutate(
    EFS_Month = as.numeric(EFS_Month),
    Event     = as.numeric(Event)
  )

# Model 1: Clinical & Molecular Baseline (NCI + Subtype + MRD)
fit_modern_base <- coxph(
  Surv(EFS_Month, Event) ~ NCI_Risk_Group + Molecular_Subtype + True_Molecular_MRD, 
  data = our_modern_data
)
our_modern_data$LP_Modern_Base <- predict(fit_modern_base, type = "lp")

# Model 2: Evolution-Enhanced Model (+ Subclonal Mutation Count)
fit_modern_evo <- coxph(
  Surv(EFS_Month, Event) ~ NCI_Risk_Group + Molecular_Subtype + True_Molecular_MRD + Subclonal_Mut_Count, 
  data = our_modern_data
)
our_modern_data$LP_Modern_Evo <- predict(fit_modern_evo, type = "lp")

message(sprintf("Discovery Cohort Concordance Index (Clinical Baseline) : %.4f", concordance(fit_modern_base)$concordance))
message(sprintf("Discovery Cohort Concordance Index (Evolution-Enhanced): %.4f", concordance(fit_modern_evo)$concordance))

# ------------------------------------------------------------------------------
# 3. 3-Tier Risk Stratification (Tertiles) & KM Curve Plotting
# ------------------------------------------------------------------------------
message("Step 2: Stratifying Discovery Cohort into 3-tier risk groups using tertiles...")

cuts <- quantile(our_modern_data$LP_Modern_Evo, probs = c(0.33, 0.66), na.rm = TRUE)

our_modern_data <- our_modern_data %>%
  mutate(
    Tier_Evo_3way = case_when(
      LP_Modern_Evo <= cuts[1] ~ "Low Risk",
      LP_Modern_Evo > cuts[2]  ~ "High Risk",
      TRUE                     ~ "Intermediate Risk"
    ),
    Tier_Evo_3way = factor(Tier_Evo_3way, levels = c("Low Risk", "Intermediate Risk", "High Risk"))
  )

fit_3way <- survfit(Surv(EFS_Month, Event) ~ Tier_Evo_3way, data = our_modern_data)

group_counts <- our_modern_data %>%
  count(Tier_Evo_3way) %>%
  mutate(label = paste0(Tier_Evo_3way, " (n=", n, ")"))

p_3way <- ggsurvplot(
  fit_3way, data = our_modern_data,
  palette     = c("#4682B4", "#6E8B3D", "#CD5C5C"), # Blue, Green, Red
  pval        = TRUE, 
  pval.method = TRUE,
  risk.table  = FALSE,
  ggtheme     = theme_classic(),
  xlab        = "Follow-up Time (Months)",
  ylab        = "Event-Free Survival Probability",
  title       = "Discovery Cohort: NCI + Subtype + MRD + Subclonal Mutation Count",
  legend.title = "",
  legend.labs  = group_counts$label
)

p_3way$plot <- p_3way$plot +
  theme(
    legend.text      = element_text(size = 10, face = "bold"),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 10, face = "bold"),
    axis.text        = element_text(size = 10, color = "black"),
    panel.border     = element_rect(color = "grey40", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

pdf_out <- file.path(output_dir, "Discovery_Cohort_3Tier_KM_Curve.pdf")
pdf(pdf_out, width = 6, height = 5.5, onefile = FALSE)
print(p_3way)
dev.off()

write.csv(our_modern_data, file.path(output_dir, "Discovery_Cohort_Model_SourceData.csv"), row.names = FALSE)

message("✅ Discovery Cohort 3-tier survival stratification complete: ", pdf_out)
