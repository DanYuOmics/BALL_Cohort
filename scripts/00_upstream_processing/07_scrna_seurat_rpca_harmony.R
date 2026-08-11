#!/usr/bin/env Rscript
# Description: scRNA-seq pipeline including 10X loading, QC filtering, RPCA integration, Harmony batch correction, and UMAP clustering.

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(future)
  library(future.apply)
  library(ggplot2)
})

# -----------------------------
# 0. Configuration & Paths
# -----------------------------
data_dir <- Sys.getenv("SCRNA_MATRIX_DIR", unset = "./data/scRNA_matrix")
out_dir  <- Sys.getenv("OUT_DIR", unset = "./Core_Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("🚀 Starting Seurat RPCA + Harmony pipeline...")

# Parallel execution settings
plan("multicore", workers = as.numeric(Sys.getenv("NTHREADS", unset = "16")))
options(future.globals.maxSize = 1024 * 1024^3)

log_msg <- function(msg) {
  message(sprintf("[%s] %s", format(Sys.time(), "%F %T"), msg))
}

# -----------------------------
# 1. Sample Loading and QC
# -----------------------------
load_sample <- function(id) {
  tryCatch({
    sample_path <- file.path(data_dir, id)
    genes_file <- file.path(sample_path, "genes.tsv")
    
    if (file.exists(genes_file)) {
      genes <- read.table(genes_file, sep = "\t", stringsAsFactors = FALSE)
      colnames(genes) <- c("ENSEMBL", "SYMBOL")
      genes$SYMBOL <- ifelse(duplicated(genes$SYMBOL) | duplicated(genes$SYMBOL, fromLast = TRUE),
                             paste0(genes$SYMBOL, "_", genes$ENSEMBL),
                             genes$SYMBOL)
      counts <- Read10X(sample_path, gene.column = 2)
      rownames(counts) <- genes$SYMBOL
    } else {
      counts <- Read10X(sample_path)
    }

    obj <- CreateSeuratObject(counts, project = id, min.cells = 3, min.features = 200)
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
  }, error = function(e) {
    log_msg(paste("⚠️ Failed loading", id, ":", conditionMessage(e)))
    NULL
  })
}

samples <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
ball.list <- future_lapply(samples, load_sample)
ball.list <- ball.list[!sapply(ball.list, is.null)]
log_msg(sprintf("✅ Loaded %d samples successfully", length(ball.list)))

# -----------------------------
# 2. Sequential Preprocessing & RPCA Integration
# -----------------------------
for (i in seq_along(ball.list)) {
  ball.list[[i]] <- NormalizeData(ball.list[[i]])
  ball.list[[i]] <- FindVariableFeatures(ball.list[[i]], selection.method = "vst", nfeatures = 2000)
  ball.list[[i]] <- ScaleData(ball.list[[i]])
  ball.list[[i]] <- RunPCA(ball.list[[i]], npcs = 50)
}

log_msg("🔍 Finding RPCA anchors...")
features <- SelectIntegrationFeatures(ball.list, nfeatures = 2500)
anchors  <- FindIntegrationAnchors(ball.list, anchor.features = features, reduction = "rpca", dims = 1:40)

log_msg("🧬 Integrating datasets...")
integrated <- IntegrateData(anchors, dims = 1:40)
DefaultAssay(integrated) <- "integrated"

# -----------------------------
# 3. Harmony Batch Correction & Clustering
# -----------------------------
log_msg("⚙️ Running Harmony correction and UMAP clustering...")
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 50)
integrated <- RunHarmony(integrated, group.by.vars = "orig.ident", max.iter.harmony = 20)
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindClusters(integrated, resolution = c(0.3, 0.6, 0.9))

# -----------------------------
# 4. Save Outputs & QC Plots
# -----------------------------
saveRDS(integrated, file.path(out_dir, "integrated_data.rds"))
write.csv(integrated@meta.data, file.path(out_dir, "integration_metadata.csv"))

pdf(file.path(out_dir, "core_clustering_plots.pdf"), width = 14)
print(DimPlot(integrated, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clustering Overview"))
print(VlnPlot(integrated, features = c("nFeature_RNA", "percent.mt"), pt.size = 0))
dev.off()

log_msg("✅ RPCA + Harmony Seurat pipeline completed successfully!")
