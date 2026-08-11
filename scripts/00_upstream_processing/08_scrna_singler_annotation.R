#!/usr/bin/env Rscript
# Description: Automated cell-type annotation using SingleR and Human Primary Cell Atlas Data reference.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
})

out_dir <- Sys.getenv("OUT_DIR", unset = "./Core_Results/singleR")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(msg) {
  message(sprintf("[%s] %s", format(Sys.time(), "%F %T"), msg))
}

log_msg("📥 Loading integrated Seurat object...")
seurat_path <- file.path(dirname(out_dir), "integrated_data.rds")
if (!file.exists(seurat_path)) stop("Missing integrated Seurat object at: ", seurat_path)

seurat_object <- readRDS(seurat_path)

log_msg("🧬 Fetching HumanPrimaryCellAtlasData reference...")
ref <- celldex::HumanPrimaryCellAtlasData()

log_msg("🔍 Running SingleR annotation...")
singler_result <- SingleR(test = GetAssayData(seurat_object, slot = "data"),
                          ref = ref,
                          labels = ref$label.main)

seurat_object$SingleR_label <- singler_result$labels

# Save Metadata & Plots
write.csv(seurat_object@meta.data, file.path(out_dir, "annotated_metadata.csv"))

pdf(file.path(out_dir, "annotation_barplot.pdf"), width = 12)
label_counts <- table(seurat_object$SingleR_label)
barplot(sort(label_counts, decreasing = TRUE), las = 2, col = "steelblue",
        main = "Cell Type Distribution by SingleR", ylab = "Cell Count")
dev.off()

pdf(file.path(out_dir, "umap_by_celltype.pdf"), width = 14)
print(DimPlot(seurat_object, group.by = "SingleR_label", label = TRUE, repel = TRUE) +
  ggtitle("UMAP Colored by Cell Type"))
dev.off()

log_msg("✅ SingleR annotation completed!")
