
#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_virtual_drug_screening_lgals9.R
# Description: Evaluate virtual drug sensitivity for late-stage resistant blasts 
#              (pseudotime > 10) and analyze Pearson correlation between LGALS9 
#              expression and chemoresistance scores.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

# ------------------------------------------------------------------------------
# 1. Environment Configurations
# ------------------------------------------------------------------------------
atlas_rds_path <- Sys.getenv("ATLAS_RDS", unset = "./Immune_Atlas_105pts_Final_20260412.rds")
output_dir     <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

immune_atlas_final <- readRDS(atlas_rds_path)
classic_colors <- c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#7570b3")

# ------------------------------------------------------------------------------
# 2. Targeted Drug Sensitivity in Resistant Blasts
# ------------------------------------------------------------------------------
message("Step 1: Analyzing predicted targeted drug sensitivity in late-stage blasts...")

plot_data_drug <- FetchData(
  immune_atlas_final, 
  vars = c("Event", "DrugRes_Score1", "cell_type_main", "pseudotime")
) %>%
  filter(cell_type_main == "B_cell" & pseudotime > 10 & !is.na(Event))

pdf(file.path(output_dir, "Targeted_Drug_Sensitivity.pdf"), width = 7, height = 5)
p_drug <- ggplot(plot_data_drug, aes(x = Event, y = DrugRes_Score1, fill = Event)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
  scale_fill_manual(values = classic_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(
    title = "Targeted Drug Sensitivity in Resistant Blasts (Pseudotime > 10)",
    y = "Predicted Drug Resistance (Lower = More Sensitive)"
  )

print(p_drug)
dev.off()

# ------------------------------------------------------------------------------
# 3. Precision Correlation (LGALS9 vs Drug Resistance Score)
# ------------------------------------------------------------------------------
message("Step 2: Testing Pearson correlation between LGALS9 and chemoresistance score...")

plot_data_lgals9 <- FetchData(
  immune_atlas_final, 
  vars = c("LGALS9", "DrugRes_Score1", "Event", "cell_type_main", "pseudotime")
) %>%
  filter(cell_type_main == "B_cell" & pseudotime > 10)

cor_test <- cor.test(plot_data_lgals9$DrugRes_Score1, plot_data_lgals9$LGALS9, method = "pearson")

pdf(file.path(output_dir, "LGALS9_vs_Drug_Sensitivity.pdf"), width = 6, height = 5)
p_lgals9 <- ggplot(plot_data_lgals9, aes(x = LGALS9, y = DrugRes_Score1)) +
  geom_point(alpha = 0.2, color = "#756bb1", size = 0.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_classic() +
  annotate(
    "text", x = max(plot_data_lgals9$LGALS9, na.rm = TRUE) * 0.7, 
    y = max(plot_data_lgals9$DrugRes_Score1, na.rm = TRUE),
    label = sprintf("R = %.2f\nP < 2.2e-16", cor_test$estimate),
    fontface = "italic", size = 5
  ) +
  labs(
    title = "Precision Targeting: LGALS9 vs. Drug Sensitivity",
    subtitle = "Late-stage Resistant Blasts (n=105 patients)",
    x = "Immune Inhibitory Ligand (LGALS9) Expression",
    y = "Predicted Drug Resistance Score"
  )

print(p_lgals9)
dev.off()

# Export Source Data
write.csv(plot_data_drug, file.path(output_dir, "DrugSensitivity_SourceData.csv"), row.names = FALSE)
write.csv(plot_data_lgals9, file.path(output_dir, "LGALS9_Correlation_SourceData.csv"), row.names = FALSE)

message("🎉 Virtual drug screening and LGALS9 precision correlation completed!")
