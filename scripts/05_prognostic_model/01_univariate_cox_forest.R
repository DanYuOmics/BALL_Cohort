#!/usr/bin/env Rscript
# ==============================================================================
# Script: 01_univariate_cox_forest.R
# Description: Run univariate Cox proportional hazards regressions 
#              for baseline clinical factors, locus-specific genomic alterations, 
#              and subclonal evolutionary metrics. Render publication-ready forest 
#              plots on a logarithmic scale with HR (95% CI) and P-value annotations.
# ==============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(openxlsx)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
message("Step 1: Setting up environment and loading augmented clinical metadata...")

clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical metadata file at: ", clinical_file)
}

clinical_data_augmented <- read_excel(clinical_file)

# ------------------------------------------------------------------------------
# 2. Data Standardization & Reference Level Locking
# ------------------------------------------------------------------------------
message("Step 2: Standardizing variables and locking factor reference levels...")

analysis_df <- clinical_data_augmented %>%
  dplyr::mutate(
    EFS_Month                = as.numeric(EFS_Month),
    event_binary             = as.numeric(Event),
    Subclonal_Mut_Count      = as.numeric(Subclonal_Mut_Count),
    Cumulative_Genomic_Score = as.numeric(Cumulative_Genomic_Score),
    
    # Lock categorical variables reference levels
    NCI_Risk_Group     = relevel(factor(NCI_Risk_Group), ref = "Low Risk (LR)"),
    True_Molecular_MRD = relevel(factor(True_Molecular_MRD), ref = "Negative (<0.01%)"),
    HSCT               = relevel(factor(HSCT), ref = "Yes"),
    
    `11q13.3_Amp`      = relevel(factor(`11q13.3_Amp`, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    `14p13.3_Del`      = relevel(factor(`14p13.3_Del`, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    `17p11.2_Amp`      = relevel(factor(`17p11.2_Amp`, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    `12q22.3_Amp`      = relevel(factor(`12q22.3_Amp`, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    CREBBP             = relevel(factor(CREBBP, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    EGR1               = relevel(factor(EGR1, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    IRF2BPL            = relevel(factor(IRF2BPL, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    PAX5               = relevel(factor(PAX5, levels = c(0, 1), labels = c("No", "Yes")), ref = "No"),
    NRAS               = relevel(factor(NRAS, levels = c(0, 1), labels = c("No", "Yes")), ref = "No")
  ) %>%
  dplyr::filter(!is.na(EFS_Month), !is.na(event_binary))

# ------------------------------------------------------------------------------
# 3. Univariate Cox Regression Loop (Core 8 Factors Configuration)
# ------------------------------------------------------------------------------
message("Step 3: Fitting univariate Cox proportional hazards models for target features...")

target_variables <- c("NCI_Risk_Group", "True_Molecular_MRD", "HSCT", 
                      "Cumulative_Genomic_Score", "CREBBP", "EGR1", "Subclonal_Mut_Count")

forest_raw_results <- data.frame()

for (var in target_variables) {
  formula_uni <- as.formula(paste("Surv(EFS_Month, event_binary) ~ `", var, "`", sep = ""))
  fit_uni <- tryCatch(coxph(formula_uni, data = analysis_df), error = function(e) NULL)
  
  if (!is.null(fit_uni)) {
    sum_uni     <- summary(fit_uni)
    coef_matrix <- sum_uni$coefficients
    ci_matrix   <- sum_uni$conf.int
    row_names   <- rownames(coef_matrix)
    
    for (i in 1:nrow(coef_matrix)) {
      clean_term <- gsub("`", "", row_names[i])
      tmp_res <- data.frame(
        Raw_Variable = var,
        Term         = clean_term,
        HR           = ci_matrix[i, "exp(coef)"],
        Lower_95     = ci_matrix[i, "lower .95"],
        Upper_95     = ci_matrix[i, "upper .95"],
        P_Value      = coef_matrix[i, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
      forest_raw_results <- rbind(forest_raw_results, tmp_res)
    }
  }
}

# ------------------------------------------------------------------------------
# 4. Clean Labels, Module Categorization & Factor Ordering
# ------------------------------------------------------------------------------
message("Step 4: Cleaning term labels and mapping module groups...")

forest_clean_data <- forest_raw_results %>%
  dplyr::mutate(
    Clean_Label = case_when(
      Term == "NCI_Risk_GroupIntermediate Risk (IR)" ~ "NCI Intermediate Risk (vs. LR)",
      Term == "NCI_Risk_GroupHigh Risk/Escalated"   ~ "NCI High Risk (vs. LR)",
      Term == "True_Molecular_MRDPositive (≥0.01%)" ~ "MRD (Positive vs. Neg)",
      Term == "HSCTYes"                               ~ "HSCT (Yes vs. No)",
      Term == "Cumulative_Genomic_Score"              ~ "Cumulative Genomic Score (per point)",
      Term == "CREBBPYes"                             ~ "CREBBP Mutation (Yes vs. No)",
      Term == "EGR1Yes"                               ~ "EGR1 Mutation (Yes vs. No)",
      Term == "Subclonal_Mut_Count"                   ~ "Subclonal Mutation Count (per mut)",
      TRUE ~ Term
    ),
    Block_Group = case_when(
      Raw_Variable %in% c("NCI_Risk_Group", "True_Molecular_MRD", "HSCT") ~ "Clinical information",
      Raw_Variable %in% c("Cumulative_Genomic_Score", "CREBBP", "EGR1") ~ "Genomic features",
      Raw_Variable %in% c("Subclonal_Mut_Count") ~ "Clonal Evolutionary Metrics"
    )
  ) %>%
  dplyr::mutate(Block_Group = factor(Block_Group, levels = c("Clonal Evolutionary Metrics", 
                                                             "Genomic features", 
                                                             "Clinical information"))) %>%
  dplyr::arrange(Block_Group, HR) %>%
  dplyr::mutate(
    y_index      = dplyr::row_number(),
    HR_CI_Text   = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", Lower_95), "-", sprintf("%.2f", Upper_95), ")"),
    P_Value_Text = ifelse(P_Value < 0.001, "< 0.001", sprintf("%.3f", P_Value))
  )

# ------------------------------------------------------------------------------
# 5. Render & Export Forest Plot Chart
# ------------------------------------------------------------------------------
message("Step 5: Rendering publication-ready ggplot2 forest plot...")

x_max_limit <- 66.28  
x_pos_hr    <- x_max_limit * 1.8
x_pos_p     <- x_max_limit * 20.0

forest_final_chart <- ggplot(forest_clean_data, aes(y = y_index)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.6) +
  geom_pointrange(
    aes(x = HR, xmin = Lower_95, xmax = Upper_95, color = Block_Group), 
    size = 0.8, fatten = 3
  ) +
  scale_color_manual(values = c("Clinical information"        = "#377eb8", 
                                "Genomic features"            = "#4daf4a", 
                                "Clonal Evolutionary Metrics" = "#e41a1c")) +
  scale_x_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10, 20, 66.28), 
    labels = c("0.1", "0.5", "1.0", "2.0", "5.0", "10.0", "20.0", "66.3")
  ) +
  scale_y_continuous(breaks = forest_clean_data$y_index, labels = forest_clean_data$Clean_Label) +
  labs(
    x = "Univariate Hazard Ratio (95% CI, Log Scale)", 
    y = NULL, 
    title = "Univariate Predictors of Event-Free Survival"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text.y     = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x     = element_text(size = 9, color = "black"),
    axis.title.x    = element_text(size = 10, face = "bold", margin = margin(t = 8)),
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 9, face = "bold"),
    plot.margin     = margin(r = 180, l = 10, t = 10, b = 10) 
  ) +
  coord_cartesian(clip = "off") +
  geom_text(aes(x = x_pos_hr, label = HR_CI_Text), hjust = 0, size = 3.4, color = "black") +
  geom_text(
    aes(x = x_pos_p, label = P_Value_Text, fontface = ifelse(P_Value < 0.05, "bold", "plain")), 
    hjust = 0, size = 3.4, color = ifelse(forest_clean_data$P_Value < 0.05, "darkred", "black")
  )

pdf_out   <- file.path(output_dir, "Univariate_Cox_Forest_Plot.pdf")
excel_out <- file.path(output_dir, "Univariate_Cox_Forest_Results.xlsx")

pdf(pdf_out, width = 10, height = 6, onefile = FALSE)
print(forest_final_chart)
dev.off()

write.xlsx(forest_clean_data, file = excel_out, sheetName = "Forest_Results", rowNames = FALSE)

message("🎉 Univariate Cox forest plot generation complete!")
message("  - PDF Output   : ", pdf_out)
message("  - Excel Matrix : ", excel_out)
