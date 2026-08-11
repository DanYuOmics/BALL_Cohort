#!/usr/bin/env Rscript
# ==============================================================================
# Script: 01_reannotate_B_subclusters.R
# Description: Load 182k B-lineage single-cell dataset, perform fine-grained 
#              cluster mapping (0-46 into HSC-like, Pro-B, Pre-B, Cycling, Mature B),
#              validate with lineage markers, and generate high-contrast UMAPs.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readxl)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
seurat_b_path <- Sys.getenv("SEURAT_B_RDS", unset = "./results_subsets/seurat_b_all_and_mature_processed.rds")
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Load Datasets & Clinical Metadata
# ------------------------------------------------------------------------------
message("Step 1: Loading 182k B-lineage Seurat object and clinical metadata...")

if (!file.exists(seurat_b_path)) stop("⚠️ Missing Seurat RDS file at: ", seurat_b_path)
seurat_b <- readRDS(seurat_b_path)

if (file.exists(clinical_file)) {
  clinical_data <- read_xlsx(clinical_file)
  message(sprintf("Loaded clinical metadata for %d patients.", nrow(clinical_data)))
}

# ------------------------------------------------------------------------------
# 3. Fine-Grained Cluster Re-annotation (Clusters 0 to 46)
# ------------------------------------------------------------------------------
message("Step 2: Mapping clusters 0-46 into refined B-lineage developmental stages...")

if ("seurat_clusters" %in% colnames(seurat_b@meta.data)) {
  Idents(seurat_b) <- seurat_b$seurat_clusters
}

refined_idents <- c(
  "1"="HSC-like_Blast", "2"="HSC-like_Blast", "8"="HSC-like_Blast", 
  "17"="HSC-like_Blast", "18"="HSC-like_Blast", "20"="HSC-like_Blast",
  "0"="Pro-B_Blast", "4"="Pro-B_Blast", "9"="Pro-B_Blast", "13"="Pro-B_Blast", 
  "19"="Pro-B_Blast", "27"="Pro-B_Blast", "29"="Pro-B_Blast", "32"="Pro-B_Blast", 
  "39"="Pro-B_Blast", "41"="Pro-B_Blast",
  "5"="Pre-B_Blast", "6"="Pre-B_Blast", "7"="Pre-B_Blast", "11"="Pre-B_Blast", 
  "12"="Pre-B_Blast", "16"="Pre-B_Blast", "21"="Pre-B_Blast", "25"="Pre-B_Blast", 
  "30"="Pre-B_Blast", "33"="Pre-B_Blast", "34"="Pre-B_Blast", "38"="Pre-B_Blast", 
  "43"="Pre-B_Blast", "45"="Pre-B_Blast",
  "3"="Cycling_Blast", 
  "10"="Mature_B", "14"="Mature_B", "15"="Mature_B", "22"="Mature_B", 
  "23"="Mature_B", "24"="Mature_B", "26"="Mature_B", "28"="Mature_B", 
  "31"="Mature_B", "35"="Mature_B", "36"="Mature_B", "37"="Mature_B", 
  "40"="Mature_B", "42"="Mature_B", "44"="Mature_B", "46"="Mature_B"
)

seurat_b <- RenameIdents(seurat_b, refined_idents)
seurat_b$cell_type_refined <- Idents(seurat_b)

refined_cols_vivid <- c(
  "HSC-like_Blast" = "#d95f02", "Pro-B_Blast" = "#36648B", 
  "Pre-B_Blast"     = "#a6cee3", "Cycling_Blast" = "#e31a1c", 
  "Mature_B"        = "#33a02c"
)

# ------------------------------------------------------------------------------
# 4. Marker DotPlot Validation & Vivid UMAP Rendering
# ------------------------------------------------------------------------------
message("Step 3: Generating marker validation DotPlots and high-contrast UMAPs...")

markers_to_check <- c("CD19", "MS4A1", "CD37", "CD27", "CD34", "DNTT", "MME", "VPREB1", "MKI67", "TOP2A")
markers_present  <- intersect(markers_to_check, rownames(seurat_b))

p_markers <- DotPlot(seurat_b, features = markers_present, dot.scale = 4) + 
  RotatedAxis() + 
  theme(axis.text.y = element_text(size = 9)) +
  labs(title = "B-lineage Marker Expression Profile (182k cells)")

pdf(file.path(output_dir, "Fig3_B_Cell_Marker_Check.pdf"), width = 10, height = 6)
print(p_markers)
dev.off()

p_umap <- DimPlot(
  seurat_b, group.by = "cell_type_refined", 
  label = TRUE, label.size = 5, repel = TRUE, 
  pt.size = 0.6, label.box = TRUE, cols = refined_cols_vivid
) + 
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  labs(title = "Refined B-lineage Landscape (182k cells)")

pdf(file.path(output_dir, "Fig3A_B_UMAP_Refined_Vivid.pdf"), width = 10, height = 8)
print(p_umap)
dev.off()

saveRDS(seurat_b, file.path(output_dir, "seurat_b_reannotated.rds"))
message("✅ Successfully saved re-annotated B-lineage Seurat object!")
