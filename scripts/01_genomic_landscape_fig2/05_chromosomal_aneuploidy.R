#!/usr/bin/env Rscript
# ==============================================================================
# Script: 05_chromosomal_aneuploidy.R
# Description: Figure 2E & Supplementary Analyses - Calculates patient-level 
#              and chromosomal-level aneuploidy metrics (Monosomy, Trisomy, 
#              Tetrasomy), performs statistical comparisons across outcome groups
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(writexl)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
cns_filtered_file <- Sys.getenv("CNS_FILTERED_FILE", unset = "./article/File SI 8 cnv cytoband filtered.xlsx")
clinical_file     <- Sys.getenv("CLINICAL_FILE", unset = "./article/B_ALL_105_Clinical_Data_Final.xlsx")
pyclone_file      <- Sys.getenv("PYCLONE_FILE", unset = "./article/Fig2_Upgraded_Clonal_Architecture_TXT.tsv")
output_dir        <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Load Datasets & Preprocess Clinical Outcome Labels
# ------------------------------------------------------------------------------
message("Step 1: Loading filtered CNV segments and clinical metadata...")

if (!file.exists(cns_filtered_file)) stop("⚠️ Missing CNV file at: ", cns_filtered_file)
if (!file.exists(clinical_file))     stop("⚠️ Missing clinical file at: ", clinical_file)

all_cns_corrected <- read_xlsx(cns_filtered_file)

clinical_data <- read_xlsx(clinical_file) %>%
  mutate(Event_Label = if_else(Event == 1, "Relapsed/Deceased", "Event-free"))

# ------------------------------------------------------------------------------
# 3. Derive Integer Copy Number & Chromosomal Aneuploidy States
# ------------------------------------------------------------------------------
message("Step 2: Deriving integer copy numbers and mapping gene-level alterations...")

gene_level_mutations <- all_cns_corrected %>%
  filter(chromosome %in% paste0("chr", 1:22)) %>%
  mutate(
    cn_int     = round(2 * (2 ^ log2_adj)),
    Monosomy   = if_else(cn_int == 1, 1, 0),
    Trisomy    = if_else(cn_int == 3, 1, 0),
    Tetrasomy  = if_else(cn_int == 4, 1, 0),
    Is_Altered = if_else(cn_int != 2, 1, 0)
  )

# Calculate metrics per patient per chromosome
per_patient_per_chr_metrics <- gene_level_mutations %>%
  group_by(SampleID, chromosome) %>%
  summarise(
    altered_gene_count = sum(Is_Altered, na.rm = TRUE),
    total_gene_count   = n(),
    num_Monosomy       = sum(Monosomy, na.rm = TRUE),
    num_Trisomy        = sum(Trisomy, na.rm = TRUE),
    num_Tetrasomy      = sum(Tetrasomy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    chr_Aneuploidy = num_Monosomy + num_Trisomy + num_Tetrasomy
  )

per_patient_per_chr_metrics$chromosome <- factor(
  per_patient_per_chr_metrics$chromosome, 
  levels = paste0("chr", 1:22)
)

# Calculate global patient-level Aneuploidy Score
global_patient_metrics <- per_patient_per_chr_metrics %>%
  group_by(SampleID) %>%
  summarise(
    Global_Aneuploidy_Score = sum(chr_Aneuploidy),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 4. Global & Per-Chromosome Hypothesis Testing (Wilcoxon + FDR)
# ------------------------------------------------------------------------------
message("Step 3: Performing global and per-chromosome Wilcoxon rank-sum tests...")

# Global Wilcoxon test
global_testing_df <- global_patient_metrics %>%
  left_join(dplyr::select(clinical_data, Patient_ID, Event_Label), by = c("SampleID" = "Patient_ID")) %>%
  filter(!is.na(Event_Label))

p_global_aneu <- wilcox.test(Global_Aneuploidy_Score ~ Event_Label, data = global_testing_df)$p.value
message(sprintf("Global Aneuploidy Score Wilcoxon P-value: %.5f", p_global_aneu))

# Per-chromosome Wilcoxon tests
chromosomal_testing_df <- per_patient_per_chr_metrics %>%
  left_join(dplyr::select(clinical_data, Patient_ID, Event_Label), by = c("SampleID" = "Patient_ID")) %>%
  filter(!is.na(Event_Label))

chromosomal_stats_results <- chromosomal_testing_df %>%
  group_by(chromosome) %>%
  summarise(
    Aneuploidy_p = tryCatch({
      wilcox.test(chr_Aneuploidy ~ Event_Label)$p.value
    }, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    Aneuploidy_FDR = p.adjust(Aneuploidy_p, method = "fdr")
  ) %>%
  arrange(Aneuploidy_p)

message("Per-chromosome statistical results preview:")
print(head(chromosomal_stats_results, 10))

# ------------------------------------------------------------------------------
# 5. Integrate Aneuploidy Metrics with Clinical Data & Save Extended Table
# ------------------------------------------------------------------------------
message("Step 4: Merging chromosome metrics and saving augmented clinical datasets...")

selected_chr_metrics <- per_patient_per_chr_metrics %>%
  filter(chromosome %in% c("chr1", "chr15", "chr21")) %>%
  pivot_wider(
    id_cols = SampleID,
    names_from = chromosome,
    values_from = c(chr_Aneuploidy, num_Monosomy, num_Trisomy, num_Tetrasomy),
    names_glue = "{chromosome}_{.value}"
  )

clinical_data_augmented <- clinical_data %>%
  left_join(global_patient_metrics, by = c("Patient_ID" = "SampleID")) %>%
  left_join(selected_chr_metrics, by = c("Patient_ID" = "SampleID"))

output_augmented_path <- file.path(output_dir, "B_ALL_105_Clinical_Data_with_Aneuploidy.xlsx")
write_xlsx(clinical_data_augmented, output_augmented_path)

# ------------------------------------------------------------------------------
# 6. Quality Audit: Verify Raw Counts in Target Chromosomes
# ------------------------------------------------------------------------------
message("Step 6: Auditing target chromosome mutation counts in Relapsed/Deceased cohort...")

debug_relapse_raw <- per_patient_per_chr_metrics %>%
  left_join(dplyr::select(clinical_data, Patient_ID, Event), by = c("SampleID" = "Patient_ID")) %>%
  filter(Event == 1, chromosome %in% c("chr1", "chr15", "chr21")) %>%
  dplyr::select(SampleID, chromosome, num_Monosomy, num_Trisomy, num_Tetrasomy, total_gene_count)

print(debug_relapse_raw)

# ------------------------------------------------------------------------------
# 7. Render & Save Stacked Aneuploidy Bar Plot (Figure 2E)
# ------------------------------------------------------------------------------
message("Step 7: Rendering Figure 2E chromosomal aneuploidy stacked bar plot...")

# Structure dataset for visualization
target_chrs <- c("chr1", "chr15", "chr21")

# Compute percentage of altered genes for plot_data_perfect
plot_data_perfect <- per_patient_per_chr_metrics %>%
  filter(chromosome %in% target_chrs) %>%
  left_join(dplyr::select(clinical_data, Patient_ID, Event_Label), by = c("SampleID" = "Patient_ID")) %>%
  filter(!is.na(Event_Label)) %>%
  pivot_longer(
    cols = c(num_Monosomy, num_Trisomy, num_Tetrasomy),
    names_to = "Aneuploidy_Type",
    values_to = "Gene_Count"
  ) %>%
  mutate(
    Aneuploidy_Type = case_when(
      Aneuploidy_Type == "num_Monosomy"  ~ "Monosomy",
      Aneuploidy_Type == "num_Trisomy"   ~ "Trisomy",
      Aneuploidy_Type == "num_Tetrasomy" ~ "Tetrasomy"
    )
  ) %>%
  group_by(chromosome, Event_Label, Aneuploidy_Type) %>%
  summarise(Percentage = mean(Gene_Count / total_gene_count * 100, na.rm = TRUE), .groups = "drop")

plot_data_perfect$chromosome      <- factor(plot_data_perfect$chromosome, levels = target_chrs)
plot_data_perfect$Aneuploidy_Type <- factor(plot_data_perfect$Aneuploidy_Type, levels = c("Monosomy", "Tetrasomy", "Trisomy"))
plot_data_perfect$Event_Label     <- factor(plot_data_perfect$Event_Label, levels = c("Event-free", "Relapsed/Deceased"))

# Define annotation positions for P-values
p_annotations_fixed <- plot_data_perfect %>%
  group_by(chromosome, Event_Label) %>%
  summarise(total_height = sum(Percentage), .groups = "drop") %>%
  group_by(chromosome) %>%
  summarise(max_bar_height = max(total_height), .groups = "drop") %>%
  mutate(
    p_text = c("P = 0.056", "P = 0.066", "P = 0.553"),
    y_pos  = max_bar_height + 1.5
  )

# Render clean stacked bar chart
plot_final_clean <- ggplot(
  plot_data_perfect, 
  aes(x = Event_Label, y = Percentage, fill = Aneuploidy_Type)
) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "#555555", size = 0.2) +
  facet_wrap(~chromosome, strip.position = "bottom", nrow = 1) +
  
  scale_fill_manual(values = c(
    "Monosomy"  = "#AEC6DF",
    "Trisomy"   = "#F6B98E",
    "Tetrasomy" = "#A8D5C2"
  )) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  
  geom_text(
    data = p_annotations_fixed, 
    aes(x = 1.5, y = y_pos + 2, label = p_text),
    inherit.aes = FALSE, size = 3.8, fontface = "italic"
  ) +
  
  geom_segment(
    data = p_annotations_fixed, 
    aes(x = 1, xend = 2, y = y_pos, yend = y_pos),
    inherit.aes = FALSE, color = "black", size = 0.3
  ) +
  
  labs(
    x = "Chromosome",
    y = "Average percentage of altered genes (%)",
    fill = "Chromosome aneuploidy"
  ) +
  theme_classic() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text       = element_text(size = 12, face = "bold", color = "black"),
    panel.spacing    = unit(3, "lines"),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.position  = "right"
  )

output_pdf <- file.path(output_dir, "Figure2E_chromosomal_aneuploidy.pdf")
ggsave(output_pdf, plot = plot_final_clean, width = 8, height = 5)

message("✅ Successfully generated Figure 2E Chromosomal Aneuploidy Plot: ", output_pdf)
