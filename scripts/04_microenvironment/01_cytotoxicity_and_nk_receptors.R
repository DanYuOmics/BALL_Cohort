#!/usr/bin/env Rscript
# ==============================================================================
# Script: 01_cytotoxicity_and_nk_receptors.R
# Description: Quantify T/NK cell cytotoxicity scores (GZMB, PRF1, GNLY, GZMA, NKG7)
#              across clinical outcome groups, and plot NK cell activation/inhibition 
#              receptor balance heatmaps.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(readxl)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
t_rds_path     <- Sys.getenv("T_CELL_RDS", unset = "./results_subsets/seurat_t_cell.rds")
nk_rds_path    <- Sys.getenv("NK_CELL_RDS", unset = "./results_subsets/seurat_nk_pure.rds")
clinical_file  <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir     <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

event_colors <- c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#7570b3")

# Load clinical metadata
if (!file.exists(clinical_file)) stop("⚠️ Missing clinical metadata file at: ", clinical_file)
clinical_data <- read_xlsx(clinical_file)

# ------------------------------------------------------------------------------
# 2. Combined T/NK Cytotoxicity Scoring
# ------------------------------------------------------------------------------
message("Step 1: Quantifying T/NK combined cytotoxicity scores...")

seurat_t  <- readRDS(t_rds_path)
seurat_nk <- readRDS(nk_rds_path)

seurat_t$cell_type_main  <- "T_cell"
seurat_nk$cell_type_main <- "NK_cell"

seurat_immune <- merge(seurat_t, y = seurat_nk, add.cell.ids = c("T", "NK"))

kill_genes <- list(Cytotoxicity = c("GZMB", "PRF1", "GNLY", "GZMA", "NKG7"))
kill_genes_present <- list(Cytotoxicity = intersect(kill_genes$Cytotoxicity, rownames(seurat_immune)))

seurat_immune <- AddModuleScore(seurat_immune, features = kill_genes_present, name = "Cytotoxicity_Score")

plot_data_kill <- FetchData(seurat_immune, vars = c("Cytotoxicity_Score1", "Event", "cell_type_main")) %>%
  filter(!is.na(Event))

pdf(file.path(output_dir, "Combined_Immune_Cytotoxicity.pdf"), width = 6, height = 5)
p_kill <- ggviolin(
  plot_data_kill, x = "Event", y = "Cytotoxicity_Score1", fill = "Event",
  palette = event_colors, add = "boxplot", add.params = list(fill = "white")
) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.3, size = 5, fontface = "italic") +
  theme_classic() +
  labs(
    title = "Immune Surveillance Collapse", 
    subtitle = "Combined NK & T-cell Cytotoxicity Score",
    y = "Cytotoxicity Score (GZMB/PRF1/GNLY)", x = ""
  ) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_kill)
dev.off()

write.csv(plot_data_kill, file.path(output_dir, "Cytotoxicity_SourceData.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. NK Cell Receptor Balance Heatmap
# ------------------------------------------------------------------------------
message("Step 2: Plotting NK cell activation/inhibition receptor balance heatmap...")

nk_receptors_fixed <- intersect(
  c("KLRC1", "KLRB1", "KLRD1", "KLRG1", "KLRK1", "NCR1", "FCGR3A"),
  rownames(seurat_nk)
)

DefaultAssay(seurat_nk) <- "RNA"
seurat_nk <- ScaleData(seurat_nk, features = nk_receptors_fixed)

pdf(file.path(output_dir, "NK_Receptor_Heatmap.pdf"), width = 10, height = 6)
p_nk_ht <- DoHeatmap(
  seurat_nk, features = nk_receptors_fixed, group.by = "Event", 
  group.colors = event_colors, size = 4, angle = 45
) +
  scale_fill_gradient2(low = "#3C8DBC", mid = "white", high = "#E64B35") +
  labs(title = "NK Activation/Inhibition Receptor Imbalance")

print(p_nk_ht)
dev.off()

message("✅ Cytotoxicity scoring and NK receptor analysis completed successfully!")
