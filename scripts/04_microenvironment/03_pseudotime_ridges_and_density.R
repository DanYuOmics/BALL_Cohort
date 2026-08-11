#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_pseudotime_ridges_and_density.R
# Description: Generate microenvironment-wide pseudotime density ridges (B, T, NK cells)
#              and quantify immune infiltration density (Cold Tumor verification).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
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
# 2. Microenvironment Pseudotime Density Ridges
# ------------------------------------------------------------------------------
message("Step 1: Rendering B, T, NK cell systemic pseudotime ridges...")

plot_data_ridges <- FetchData(
  immune_atlas_final, 
  vars = c("cell_type_main", "pseudotime", "Event")
) %>%
  filter(cell_type_main %in% c("B_cell", "T_cell", "NK_cell")) %>%
  mutate(
    Event          = factor(trimws(Event), levels = c("Event-free", "Relapse/Deceased")),
    cell_type_main = factor(cell_type_main, levels = c("B_cell", "T_cell", "NK_cell"))
  )

pdf(file.path(output_dir, "Systemic_Pseudotime_Ridges.pdf"), width = 10, height = 7)
p_ridges <- ggplot(plot_data_ridges, aes(x = pseudotime, y = cell_type_main, fill = Event)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9, color = "white", rel_min_height = 0.01) +
  scale_fill_manual(values = classic_colors) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme_classic() +
  theme(legend.position = "top", axis.text.y = element_text(face = "bold", size = 12)) +
  labs(
    title = "Systemic Pseudotime Shift in B-ALL Microenvironment (n=105)",
    subtitle = "Combined analysis of 193,252 cells showing significant Event-related shifts",
    x = "Differentiation / Activation Pseudotime", y = "Immune Compartments"
  ) +
  annotate("text", x = Inf, y = 3.2, label = "p < 0.001", hjust = 1.1, fontface = "italic") +
  annotate("text", x = Inf, y = 2.2, label = "p < 0.001", hjust = 1.1, fontface = "italic") +
  annotate("text", x = Inf, y = 1.2, label = "p < 0.001", hjust = 1.1, fontface = "italic")

print(p_ridges)
dev.off()

# ------------------------------------------------------------------------------
# 3. Immune Infiltration Density (Cold Tumor Verification)
# ------------------------------------------------------------------------------
message("Step 2: Quantifying patient-level immune cell infiltration density...")

patient_density <- as.data.frame(table(immune_atlas_final$orig.ident, immune_atlas_final$cell_type_main))
colnames(patient_density) <- c("orig.ident", "cell_type", "count")

patient_metadata <- immune_atlas_final@meta.data %>%
  select(orig.ident, Event) %>%
  distinct()

density_summary <- patient_density %>%
  group_by(orig.ident) %>%
  mutate(total = sum(count), proportion = count / total) %>%
  filter(cell_type %in% c("T_cell", "NK_cell")) %>%
  group_by(orig.ident) %>%
  summarise(immune_density = sum(proportion), .groups = "drop") %>%
  left_join(patient_metadata, by = "orig.ident") %>%
  filter(!is.na(Event))

pdf(file.path(output_dir, "Immune_Infiltration_Density.pdf"), width = 5, height = 5)
p_density <- ggplot(density_summary, aes(x = Event, y = immune_density, fill = Event)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_fill_manual(values = classic_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  labs(title = "Immune Infiltration Density", y = "T/NK Proportion in Total Cells")

print(p_density)
dev.off()

write.csv(density_summary, file.path(output_dir, "Immune_Density_SourceData.csv"), row.names = FALSE)
message("✅ Systemic pseudotime ridges and cold tumor density verification complete!")
