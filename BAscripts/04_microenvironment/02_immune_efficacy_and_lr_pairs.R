#!/usr/bin/env Rscript
# ==============================================================================
# Script: 02_immune_efficacy_and_lr_pairs.R
# Description: Generate Immune Efficacy Mapping (Blast Resistance vs Immune 
#              Cytotoxicity) and plot Ligand-Receptor inhibitory pair heatmaps 
#              (LGALS9-HAVCR2, PVR-TIGIT, CD80/86-CTLA4).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(readxl)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
atlas_rds_path <- Sys.getenv("ATLAS_RDS", unset = "./Immune_Atlas_105pts_Final_20260412.rds")
clinical_file  <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir     <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

classic_colors <- c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#7570b3")

# Load integrated Immune Atlas object
if (!file.exists(atlas_rds_path)) stop("⚠️ Missing Immune Atlas object at: ", atlas_rds_path)
immune_atlas_final <- readRDS(atlas_rds_path)

# ------------------------------------------------------------------------------
# 2. Immune Efficacy Mapping (Sword vs Shield)
# ------------------------------------------------------------------------------
message("Step 1: Mapping B-cell Resistance (Shield) vs T/NK Cytotoxicity (Sword)...")

# Blast Resistance Score (The Shield)
df_shield <- immune_atlas_final@meta.data %>%
  filter(cell_type_main == "B_cell") %>%
  group_by(orig.ident, Event) %>%
  summarise(avg_resistance = mean(as.numeric(DrugRes_Score1), na.rm = TRUE), .groups = "drop")

# T/NK Cytotoxicity Score (The Sword)
df_sword <- immune_atlas_final@meta.data %>%
  filter(cell_type_main %in% c("T_cell", "NK_cell")) %>%
  group_by(orig.ident) %>%
  summarise(avg_cytotoxicity = mean(as.numeric(Cytotoxicity_Score1), na.rm = TRUE), .groups = "drop")

patient_summary_fixed <- inner_join(df_shield, df_sword, by = "orig.ident") %>%
  filter(!is.na(Event))

p_val_sword <- wilcox.test(avg_cytotoxicity ~ Event, data = patient_summary_fixed)$p.value
p_label     <- formatC(p_val_sword, format = "e", digits = 1)

pdf(file.path(output_dir, "Immune_Efficacy_Mapping.pdf"), width = 7, height = 6)
p_eff <- ggplot(patient_summary_fixed, aes(x = avg_resistance, y = avg_cytotoxicity, color = Event)) +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(aes(fill = Event), geom = "polygon", alpha = 0.1, linetype = 0) +
  scale_color_manual(values = classic_colors) +
  scale_fill_manual(values = classic_colors) +
  theme_classic() +
  annotate(
    "text", x = min(patient_summary_fixed$avg_resistance), 
    y = max(patient_summary_fixed$avg_cytotoxicity), 
    label = paste0("p = ", p_label), hjust = 0, fontface = "italic", size = 5
  ) +
  theme(aspect.ratio = 1, legend.position = "top", axis.title = element_text(face = "bold", size = 12)) +
  labs(
    title = "Immune Efficacy Mapping (n=105)",
    x = "B-cell Resistance (The Shield)", 
    y = "Immune Cytotoxicity (The Sword)"
  )

print(p_eff)
dev.off()

write.csv(patient_summary_fixed, file.path(output_dir, "Immune_Efficacy_SourceData.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Inhibitory Ligand-Receptor Interaction Heatmap
# ------------------------------------------------------------------------------
message("Step 2: Plotting Inhibitory Ligand-Receptor Interaction Heatmap...")

b_ligands        <- c("LGALS9", "PVR", "CD80", "CD86")
immune_receptors <- c("HAVCR2", "TIGIT", "CTLA4", "LAG3")
target_pairs     <- intersect(c(b_ligands, immune_receptors), rownames(immune_atlas_final))

avg_exp <- AverageExpression(
  immune_atlas_final, 
  features = target_pairs, 
  group.by = c("cell_type_main", "Event"),
  layer = "data"
)$RNA

mat <- t(scale(t(avg_exp)))
col_fun_nc <- colorRamp2(c(-1.5, 0, 1.5), c("#2166ac", "#f7f7f7", "#b2182b"))

pdf(file.path(output_dir, "Inhibitory_Axe_Heatmap.pdf"), width = 9, height = 6)
ht <- Heatmap(
  mat, name = "Z-score", col = col_fun_nc,
  column_title = "Inhibitory Interaction Landscape",
  rect_gp = gpar(col = "white", lwd = 1.2),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11, fontface = "italic"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows = TRUE, cluster_columns = FALSE,
  heatmap_legend_param = list(title = "Z-score", title_position = "topleft")
)

draw(ht, padding = unit(c(2, 10, 2, 2), "mm"))
dev.off()

write.csv(avg_exp, file.path(output_dir, "LR_Interaction_SourceData.csv"))
message("✅ Immune efficacy mapping and LR interaction heatmap complete!")
