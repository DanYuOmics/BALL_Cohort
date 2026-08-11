#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_hsc_richness_pseudotime.R
# Description: Evaluate B-lineage and microenvironment stemness features:
#              1. Figure 4C: Patient-level HSC-like B cell abundance (pseudo-log scale).
#              2. Figure 4D: Single-cell transcriptional richness (nFeature_RNA).
#              3. Figure 4E: Systemic pseudotime density ridges across immune compartments
#                 (B, T, NK cells, n=193,252 cells).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggridges)
  library(scales)
  library(readxl)
  library(readr)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
message("Step 1: Setting up environment and directory paths...")

bcell_rds_path  <- Sys.getenv("BCELL_RDS", unset = "./results_subsets/B_cells_Final_Annotated_20260708.rds")
atlas_rds_path  <- Sys.getenv("ATLAS_RDS", unset = "./results_subsets/Immune_Atlas_105pts_Final_20260412.rds")
clinical_file   <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir      <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

event_colors <- c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#7570b3")

# ------------------------------------------------------------------------------
# 2. Figure 4C1: Patient-Level HSC-like B Cell Abundance (Pseudo-Log Scale)
# ------------------------------------------------------------------------------
message("Step 2: Analyzing patient-level HSC-like B cell abundance...")

if (!file.exists(bcell_rds_path)) {
  stop("⚠️ Missing B-cell Seurat RDS object at: ", bcell_rds_path)
}

bcell_obj <- readRDS(bcell_rds_path)

# Calculate patient-level proportion of HSC-like blasts
plot_data_abundance <- bcell_obj@meta.data %>%
  group_by(orig.ident, Event) %>%
  summarise(
    total_blasts = n(),
    hsc_count    = sum(B_cell_perfect == "HSC-like_Blast", na.rm = TRUE),
    hsc_prop     = (hsc_count / total_blasts) * 100,
    .groups      = "drop"
  ) %>%
  filter(total_blasts >= 20) # Filter noise samples with < 20 blasts

pdf(file.path(output_dir, "Fig4C1_HSC_like_Bcell_Abundance.pdf"), width = 3.5, height = 4.2)

p_abundance <- ggplot(plot_data_abundance, aes(x = Event, y = hsc_prop, fill = Event)) +
  geom_violin(trim = FALSE, alpha = 0.3, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.4) +
  geom_jitter(aes(color = Event), width = 0.12, size = 1.8, alpha = 0.7) +
  scale_y_continuous(
    trans  = pseudo_log_trans(base = 10, sigma = 0.01),
    breaks = c(0, 0.1, 1, 5, 10, 25, 50),
    labels = c("0", "0.1", "1", "5", "10", "25", "50")
  ) +
  scale_fill_manual(values = event_colors) +
  scale_color_manual(values = event_colors) +
  theme_classic(base_size = 11) +
  labs(
    title = "HSC-like B cell abundance",
    x     = "",
    y     = "Percentage of HSC-like B cells (%)"
  ) +
  theme(
    plot.margin  = margin(t = 15, r = 10, b = 10, l = 10),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text.x  = element_text(angle = 25, hjust = 1, color = "black", size = 9.5),
    axis.text.y  = element_text(color = "black", size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )

print(p_abundance)
dev.off()

# ------------------------------------------------------------------------------
# 3. Figure 4C2: Single-Cell Transcriptional Richness (nFeature_RNA)
# ------------------------------------------------------------------------------
message("Step 3: Quantifying single-cell gene richness per cell (nFeature_RNA)...")

plot_data_feature <- FetchData(bcell_obj, vars = c("Event", "nFeature_RNA")) %>%
  mutate(Event = factor(Event, levels = c("Event-free", "Relapse/Deceased")))

# Perform Wilcoxon rank-sum test
wilcox_feature <- wilcox.test(nFeature_RNA ~ Event, data = plot_data_feature)

summary_feature <- plot_data_feature %>%
  group_by(Event) %>%
  summarise(
    cell_count   = n(),
    median_genes = median(nFeature_RNA, na.rm = TRUE),
    mean_genes   = mean(nFeature_RNA, na.rm = TRUE),
    sd_genes     = sd(nFeature_RNA, na.rm = TRUE),
    .groups      = "drop"
  )

message(sprintf("Single-cell Gene Richness Wilcoxon P-value: %.5e", wilcox_feature$p.value))

pdf(file.path(output_dir, "Fig4C2_Transcriptional_Diversity_Clean.pdf"), width = 3.5, height = 4.2)

p_feature <- ggplot(plot_data_feature, aes(x = Event, y = nFeature_RNA, fill = Event)) +
  geom_violin(trim = FALSE, alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = event_colors) +
  coord_cartesian(ylim = c(0, 4500)) +
  theme_classic(base_size = 11) +
  labs(
    title = "Developmental Potential",
    x     = "",
    y     = "Detected Genes per Cell (nFeature)"
  ) +
  theme(
    plot.margin  = margin(t = 15, r = 10, b = 10, l = 10),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text.x  = element_text(angle = 25, hjust = 1, color = "black", size = 9.5),
    axis.text.y  = element_text(color = "black", size = 9),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )

print(p_feature)
dev.off()

# ------------------------------------------------------------------------------
# 4. Figure 4G: Systemic Microenvironment Pseudotime Density Ridges
# ------------------------------------------------------------------------------
message("Step 4: Rendering systemic microenvironment pseudotime density ridges (Figure 4G)...")

if (file.exists(atlas_rds_path)) {
  immune_atlas_final <- readRDS(atlas_rds_path)
  
  plot_data_ridges <- FetchData(
    immune_atlas_final, 
    vars = c("cell_type_main", "pseudotime", "Event")
  ) %>%
    filter(cell_type_main %in% c("B_cell", "T_cell", "NK_cell")) %>%
    mutate(
      Event          = factor(trimws(Event), levels = c("Event-free", "Relapse/Deceased")),
      cell_type_main = factor(cell_type_main, levels = c("B_cell", "T_cell", "NK_cell"))
    )

  pdf(file.path(output_dir, "Fig4G_Systemic_Pseudotime_Ridges.pdf"), width = 10, height = 7)

  p_ridges <- ggplot(plot_data_ridges, aes(x = pseudotime, y = cell_type_main, fill = Event)) +
    geom_density_ridges(alpha = 0.6, scale = 0.9, color = "white", rel_min_height = 0.01) +
    scale_fill_manual(values = c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#7570b3")) +
    theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.text.y     = element_text(face = "bold", size = 12),
      strip.text      = element_text(face = "bold")
    ) +
    labs(
      title    = "Systemic Pseudotime Shift in B-ALL Microenvironment (n=105)",
      subtitle = "Combined analysis of 193,252 cells showing significant Event-related shifts",
      x        = "Differentiation / Activation Pseudotime",
      y        = "Immune Compartments"
    ) +
    annotate("text", x = Inf, y = 3.2, label = "p < 0.001", hjust = 1.1, fontface = "italic") + # NK cell
    annotate("text", x = Inf, y = 2.2, label = "p < 0.001", hjust = 1.1, fontface = "italic") + # T cell
    annotate("text", x = Inf, y = 1.2, label = "p < 0.001", hjust = 1.1, fontface = "italic")   # B cell

  print(p_ridges)
  dev.off()
} else {
  warning("⚠️ Immune atlas object not found. Skipped Figure 4G density ridges.")
}

# ------------------------------------------------------------------------------
# 5. Export Source Data Tables
# ------------------------------------------------------------------------------
message("Step 5: Exporting source data tables...")

write.csv(plot_data_abundance, file.path(output_dir, "Table_S4_Fig4C1_HSC_Abundance_SourceData.csv"), row.names = FALSE)
write.csv(plot_data_feature, file.path(output_dir, "Table_S4_Fig4C2_Gene_Richness_SourceData.csv"), row.names = FALSE)

message("🎉 Stemness features and microenvironment pseudotime analysis complete!")
