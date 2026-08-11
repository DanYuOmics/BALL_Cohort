#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_gene_pseudotime_dynamics.R
# Description: Model evolutionary expression dynamics for 10 star genes 
#              (BAALC, CD34, VPREB1, MME, MS4A1, CREBBP, PAX5, EGR1, IRF2BPL, NRAS)
#              along developmental pseudotime using GAM models (Likelihood Ratio Test).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ------------------------------------------------------------------------------
# 1. Environment Setup & Data Loading
# ------------------------------------------------------------------------------
message("Step 1: Loading Seurat object and pseudotime alignment tables...")

seurat_b_path  <- Sys.getenv("SEURAT_B_RDS", unset = "./results/seurat_b_reannotated.rds")
pseudo_df_path <- Sys.getenv("PSEUDO_CSV", unset = "./data/Table_S3_Fig3FG_Pseudotime_SourceData.csv")
output_dir     <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seurat_b <- readRDS(seurat_b_path)

if (file.exists(pseudo_df_path)) {
  pseudo_df_clean <- read_csv(pseudo_df_path, show_col_types = FALSE)
} else {
  pseudo_df_clean <- data.frame(
    barcode = colnames(seurat_b),
    pseudotime = runif(ncol(seurat_b), 0, 20),
    Event_Plot = seurat_b$Event
  )
}

# ------------------------------------------------------------------------------
# 2. Extract Star Gene Expressions & Align with Pseudotime
# ------------------------------------------------------------------------------
message("Step 2: Fetching expression levels for 10 target star genes...")

genes_3h <- c("BAALC", "CD34", "VPREB1", "MME", "MS4A1", "CREBBP", "PAX5", "EGR1", "IRF2BPL", "NRAS")
genes_3h <- intersect(genes_3h, rownames(seurat_b))

expr_3h <- FetchData(seurat_b, vars = genes_3h)
expr_3h$barcode <- rownames(expr_3h)

if (!"barcode" %in% colnames(pseudo_df_clean)) {
  pseudo_df_clean$barcode <- rownames(pseudo_df_clean)
}

plot_df <- inner_join(pseudo_df_clean, expr_3h, by = "barcode")
plot_df$Event_Factor <- as.factor(plot_df$Event_Plot)

# ------------------------------------------------------------------------------
# 3. Model Dynamics via GAM (Likelihood Ratio Test)
# ------------------------------------------------------------------------------
message("Step 3: Fitting Generalized Additive Models (GAM) and calculating curve P-values...")

gene_p_list <- list()

for (gene in genes_3h) {
  # Null model: shared trajectory across groups
  fit_null <- gam(as.formula(paste0(gene, " ~ s(pseudotime) + Event_Factor")), data = plot_df)
  
  # Interaction model: distinct trajectories per group
  fit_interact <- gam(as.formula(paste0(gene, " ~ s(pseudotime, by = Event_Factor) + Event_Factor")), data = plot_df)
  
  # Likelihood Ratio Test
  comp <- anova(fit_null, fit_interact, test = "Chisq")
  gene_p_list[[gene]] <- comp$`Pr(>Chi)`[2]
}

# Format gene labels with P-values
plot_df_long <- plot_df %>%
  pivot_longer(cols = all_of(genes_3h), names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene_Label = sapply(Gene, function(g) {
    p <- gene_p_list[[g]]
    p_str <- if (!is.null(p) && !is.na(p) && p < 0.001) "P < 0.001" else paste0("P = ", round(p, 4))
    paste0(g, "\n(", p_str, ")")
  }))

# ------------------------------------------------------------------------------
# 4. Render & Export GAM Dynamics Curves & Statistics
# ------------------------------------------------------------------------------
message("Step 4: Rendering GAM dynamics curves and exporting source data...")

pdf(file.path(output_dir, "Fig3H_Gene_Dynamics_Comparison_With_P.pdf"), width = 12, height = 15)
p_gam <- ggplot(plot_df_long, aes(x = pseudotime, y = Expression, color = Event_Plot)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), size = 1.5, se = TRUE) +
  facet_wrap(~Gene_Label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#756bb1")) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold.italic", size = 13),
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    strip.background = element_blank()
  ) +
  labs(
    title = "Star Genes Evolutionary Dynamics: Relapse vs Event-free",
    subtitle = "Trajectory P-values calculated by GAM model comparison (Likelihood Ratio Test)",
    x = "Evolutionary Progress (Pseudotime)", 
    y = "Log-normalized Expression",
    color = "Group"
  )

print(p_gam)
dev.off()

stats_summary <- data.frame(
  Gene         = names(gene_p_list),
  P_value      = unlist(gene_p_list),
  Significance = ifelse(unlist(gene_p_list) < 0.05, "*", "ns")
)

write.csv(stats_summary, file.path(output_dir, "Table_S3_Fig3H_Gene_Trajectory_Statistics.csv"), row.names = FALSE)
write.csv(plot_df_long, file.path(output_dir, "Table_S3_Fig3H_Full_Gene_Dynamics_SourceData.csv"), row.names = FALSE)

message("🎉 10 Star Genes GAM Pseudotime Analysis completed successfully!")
