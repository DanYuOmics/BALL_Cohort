#!/usr/bin/env Rscript
# ==============================================================================
# Script: 02_cell_proportion_comparison.R
# Description: Compare B-lineage cell proportions between clinical outcome groups 
#              (sample-level with zero-filling), compute single-cell chemoresistance 
#              module scores via AddModuleScore, and export source data tables.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
  library(readxl)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
seurat_b_path <- Sys.getenv("SEURAT_B_RDS", unset = "./results/seurat_b_reannotated.rds")
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

seurat_b <- readRDS(seurat_b_path)
clinical_data <- read_xlsx(clinical_file)

# ------------------------------------------------------------------------------
# 2. Patient-level Cell Proportion Comparison (with Zero-filling)
# ------------------------------------------------------------------------------
message("Step 1: Calculating patient-level cell proportions with complete zero-filling...")

prop_master <- seurat_b@meta.data %>%
  filter(!is.na(Event) & Event != "Unknown") %>%
  mutate(Event_Clean = case_when(
    Event %in% c("Relapse/Deceased", "1", "Relapse/Died") ~ "Relapse/Deceased",
    Event %in% c("Event-free", "0") ~ "Event-free",
    TRUE ~ as.character(Event)
  )) %>%
  # Count cells per patient per cell type
  count(orig.ident, Event_Clean, cell_type_refined, name = "cell_count") %>%
  # Zero-filling: ensure patients with 0 cells in a subtype get 0% proportion
  complete(nesting(orig.ident, Event_Clean), cell_type_refined, fill = list(cell_count = 0)) %>%
  group_by(orig.ident) %>%
  mutate(
    total_cells = sum(cell_count),
    Proportion  = cell_count / total_cells
  ) %>%
  ungroup()

# Wilcoxon rank-sum test across patient groups
stat_final <- prop_master %>%
  group_by(cell_type_refined) %>%
  filter(n_distinct(Event_Clean) == 2) %>%
  wilcox_test(Proportion ~ Event_Clean) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  add_xy_position(x = "cell_type_refined", dodge = 0.8)

n_ef  <- n_distinct(prop_master$orig.ident[prop_master$Event_Clean == "Event-free"])
n_rel <- n_distinct(prop_master$orig.ident[prop_master$Event_Clean == "Relapse/Deceased"])

# Render Proportion Comparison Plot
pdf(file.path(output_dir, "Fig3C_B_Proportion_Comparison.pdf"), width = 10, height = 4.5)
p_prop <- ggplot(prop_master, aes(x = cell_type_refined, y = Proportion, fill = Event_Clean)) +
  geom_violin(alpha = 0.2, position = position_dodge(0.8), color = NA, trim = TRUE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.25, position = position_dodge(0.8), color = "grey30", linewidth = 0.4) +
  geom_jitter(aes(color = Event_Clean), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1.2, alpha = 0.7) +
  stat_pvalue_manual(stat_final, label = "p.format", hide.ns = FALSE, tip.length = 0) +
  scale_fill_manual(
    values = c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#756bb1"),
    labels = c(paste0("Event-free (n=", n_ef, ")"), paste0("Relapse/Deceased (n=", n_rel, ")"))
  ) +
  scale_color_manual(values = c("Event-free" = "#147a5c", "Relapse/Deceased" = "#5a528a")) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Cell Type Proportion Comparison across Outcome Groups", x = "", y = "Proportion per Patient (%)", fill = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 10, color = "black", face = "bold"), legend.position = "top") +
  guides(color = "none")

print(p_prop)
dev.off()

# ------------------------------------------------------------------------------
# 3. Drug Resistance Scoring via AddModuleScore
# ------------------------------------------------------------------------------
message("Step 2: Computing single-cell drug resistance module scores (AddModuleScore)...")

res_genes <- list(DrugResistance = c("ABCB1", "ABCC1", "BCL2", "BIRC5", "MYC", "TYMS", "DHFR"))
res_genes_present <- list(DrugResistance = intersect(res_genes$DrugResistance, rownames(seurat_b)))

seurat_b <- AddModuleScore(seurat_b, features = res_genes_present, name = "DrugRes_Score")
res_data_cell_level <- FetchData(seurat_b, vars = c("cell_type_refined", "Event", "DrugRes_Score1"))

pdf(file.path(output_dir, "Fig3E_CellLevel_Resistance_Comparison.pdf"), width = 10, height = 6)
p_res <- ggviolin(
  res_data_cell_level, x = "cell_type_refined", y = "DrugRes_Score1", 
  fill = "Event", palette = c("Event-free" = "#1b9e77", "Relapse/Deceased" = "#756bb1"),
  add = "boxplot", add.params = list(fill = "white", width = 0.1)
) +
  stat_compare_means(aes(group = Event), label = "p.signif", method = "wilcox.test") +
  labs(title = "Single-cell Chemoresistance Potential across B-lineage Subtypes", x = "", y = "Drug Resistance Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))

print(p_res)
dev.off()

# ------------------------------------------------------------------------------
# 4. Export Source Data Tables
# ------------------------------------------------------------------------------
message("Step 3: Exporting source data tables...")
write.csv(prop_master, file.path(output_dir, "Table_S3_Fig3C_B_Proportion_SourceData.csv"), row.names = FALSE)
write.csv(res_data_cell_level, file.path(output_dir, "Table_S3_Fig3E_Resistance_SourceData.csv"), row.names = FALSE)

message("✅ Cell proportion and drug resistance analysis completed successfully!")
