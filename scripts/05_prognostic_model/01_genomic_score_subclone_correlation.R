#!/usr/bin/env Rscript
# ==============================================================================
# Script: 01_genomic_score_subclone_correlation.R
# Description: Evaluate the correlation between Cumulative Genomic Score 
#              and Subclonal Mutation Count using Spearman rank correlation.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(readxl)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing augmented clinical metadata file at: ", clinical_file)
}

clinical_data_augmented <- read_excel(clinical_file) %>%
  mutate(
    EFS_Month            = as.numeric(EFS_Month),
    Event                = as.numeric(Event),
    Subclonal_Mut_Count  = as.numeric(Subclonal_Mut_Count)
  )

# ------------------------------------------------------------------------------
# 2. Spearman Correlation Analysis and Scatter Plot Generation
# ------------------------------------------------------------------------------
message("Step 1: Calculating Spearman correlation between genomic score and subclonal mutation count...")

p_corr <- ggplot(clinical_data_augmented, aes(x = Cumulative_Genomic_Score, y = Subclonal_Mut_Count)) +
  geom_jitter(color = "#4A6FA5", alpha = 0.6, width = 0.15, height = 0.2, size = 2.5) +
  geom_smooth(method = "lm", color = "#E07A5F", fill = "#F4A261", alpha = 0.15, linewidth = 1.2) +
  stat_cor(
    method = "spearman", label.x = 0.5, 
    label.y = max(clinical_data_augmented$Subclonal_Mut_Count, na.rm = TRUE) * 0.9,
    size = 3.6, fontface = "italic", color = "black"
  ) +
  theme_classic() +
  labs(
    title = "Association Between Cumulative Genomic Score \nand Subclonal Mutation Count",
    x = "Cumulative Genomic Score",
    y = "Subclonal Mutation Count"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.3)
  )

pdf_out <- file.path(output_dir, "Score_vs_SubclonalMut_Correlation.pdf")
ggsave(pdf_out, plot = p_corr, width = 5.2, height = 4.5)

write.csv(
  clinical_data_augmented %>% select(Patient_ID, Cumulative_Genomic_Score, Subclonal_Mut_Count),
  file.path(output_dir, "Score_vs_SubclonalMut_SourceData.csv"),
  row.names = FALSE
)

message("✅ Score vs Subclonal Mutation correlation analysis complete: ", pdf_out)
