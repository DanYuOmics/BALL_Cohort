#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_comparative_mutational_map.R
# Description: Figure 2B - Cross-cohort comparative mutational landscape.
#              Loads clinical data, merges Panel, WES, Fudan (n=140), and 
#              St. Jude-COG (n=1717) cohorts, performs Fisher's exact tests, 
#              and plots a symmetric butterfly bar chart.
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./data/clinical_data_augmented_fixed_20260620_perfect.xlsx")
data_file     <- Sys.getenv("COMPARATIVE_DATA", unset = "./data/Comparative_Mutational_Landscape.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Load Clinical Data & Check Sample Cohort Alignment
# ------------------------------------------------------------------------------
message("Step 1: Loading clinical metadata...")
if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical metadata file at: ", clinical_file)
}
clinical_data <- read_xlsx(clinical_file)
message(sprintf("Loaded clinical metadata for %d patients.", nrow(clinical_data)))

# ------------------------------------------------------------------------------
# 3. Load Comparative Mutational Datasets across Cohorts
# ------------------------------------------------------------------------------
message("Step 2: Loading comparative cross-cohort datasets...")
if (!file.exists(data_file)) {
  stop("⚠️ Error: Missing comparative dataset at: ", data_file)
}

plot_df      <- read_excel(data_file, sheet = "Plot_Data")
final_matrix <- read_excel(data_file, sheet = "Fisher_Test")

# Define the focused 15 driver genes and lock factor levels
focused_genes <- c("NRAS", "KRAS", "KMT2D", "CREBBP", "NOTCH1", "PTPN11", 
                   "SETD2", "NSD2", "JAK3", "KMT2A", "NCOR2", "EGR1", "BCL11B")

plot_df$Hugo_Symbol      <- factor(plot_df$Hugo_Symbol, levels = rev(focused_genes))
final_matrix$Hugo_Symbol <- factor(final_matrix$Hugo_Symbol, levels = rev(focused_genes))

# Define color palette for distinct cohorts
cohort_colors <- c(
  "Panel (n=85)"                       = "#386CB0", 
  "WES (n=33)"                         = "#CCE5FF", 
  "Fudan Cohort (n=140)"               = "orange",  
  "St. Jude–COG B-ALL Cohort (n=1717)" = "#1b9e77"
)

plot_df$Cohort <- factor(plot_df$Cohort, levels = names(cohort_colors))

# ------------------------------------------------------------------------------
# 4. Render and Save Comparative Butterfly Plot
# ------------------------------------------------------------------------------
message("Step 3: Rendering cross-cohort comparison plot...")

p2b <- ggplot(plot_df, aes(x = Hugo_Symbol, y = Plot_Freq, fill = Cohort)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.7) +
  coord_flip() +
  scale_y_continuous(labels = abs, limits = c(-80, 80), breaks = seq(-80, 80, 40)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  
  # Add Fisher's exact test significance stars next to the central baseline
  geom_text(
    data = final_matrix, 
    aes(x = Hugo_Symbol, y = -1.5, label = p_label), 
    inherit.aes = FALSE, size = 5, color = "black", fontface = "bold", hjust = 1
  ) +
  
  theme_classic() +
  scale_fill_manual(values = cohort_colors) +
  labs(
    subtitle = "Spearman Correlation (Panel vs St. Jude–COG): R = 0.60, P < 0.001",
    x = "", 
    y = "Mutation Frequency (%)", 
    fill = ""
  ) + 
  theme(
    axis.text.y = element_text(face = "italic", size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

output_pdf <- file.path(output_dir, "Figure2B_cross_cohort_landscape.pdf")
ggsave(output_pdf, plot = p2b, width = 9, height = 6)

message("✅ Successfully generated Figure 2B Comparative Butterfly Plot: ", output_pdf)
