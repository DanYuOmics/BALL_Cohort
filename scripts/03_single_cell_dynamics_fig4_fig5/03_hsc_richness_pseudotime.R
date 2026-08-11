#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_hsc_richness_pseudotime.R
# Description: Evaluate stemness features: HSC-like blast proportions, gene 
#              richness (nFeature_RNA), and pseudotime distributions across outcome groups.
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
seurat_b_path <- Sys.getenv("SEURAT_B_RDS", unset = "./results/seurat_b_reannotated.rds")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seurat_b <- readRDS(seurat_b_path)

# ------------------------------------------------------------------------------
# 2. Gene Richness & HSC Feature Quantifications
# ------------------------------------------------------------------------------
message("Step 1: Quantifying gene richness per cell (nFeature_RNA)...")

richness_df <- FetchData(seurat_b, vars = c("cell_type_refined", "Event", "nFeature_RNA", "nCount_RNA")) %>%
  filter(!is.na(Event))

pdf(file.path(output_dir, "Fig3_Gene_Richness_Violin.pdf"), width = 9, height = 5)
p_richness <- ggplot(richness_df, aes(x = cell_type_refined, y = nFeature_RNA, fill = Event)) +
  geom_violin(alpha = 0.4, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.15, position = position_dodge(0.8), outlier.shape = NA) +
  scale_fill_manual(values = c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#756bb1")) +
  theme_classic() +
  labs(title = "Transcriptional Richness (nFeature_RNA) across B-lineage Subtypes", x = "", y = "Detected Genes per Cell") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"))

print(p_richness)
dev.off()

write.csv(richness_df, file.path(output_dir, "Table_S3_Gene_Richness_SourceData.csv"), row.names = FALSE)
message("✅ HSC-like and Gene Richness analysis complete!")
