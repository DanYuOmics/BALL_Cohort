#!/usr/bin/env Rscript
# ==============================================================================
# Script: 03_subclone_evolution_figures.R
# Description: Module 3 - Comprehensive PyClone downstream analysis (excluding
#              diversity indices):
#              1. Figure 3A: Population-level clonal architecture stratification (Violin)
#              2. Figure 3B: Driver gene clonal niche preferences (Trunk vs Branch)
#              3. Figure 3C: Multi-omics threshold collapse analysis (Hits Score 0-3 vs >=4)
#              4. Figure 3D: Extreme index case study 
#              5. Exports multi-sheet statistical tables for cluster prevalence and niche stats.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# ------------------------------------------------------------------------------
# 1. Environment Setup & Clinical Metadata Ingestion
# ------------------------------------------------------------------------------
message("Step 1: Setting up environment and loading clinical metadata...")

clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_data_augmented_fixed_20260620_perfect.xlsx")
pyclone_txt   <- Sys.getenv("PYCLONE_TXT_OUT", unset = "./article/Fig2_Upgraded_Clonal_Architecture_TXT.tsv")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical metadata file at: ", clinical_file)
}

clinical_data_raw_updated <- read_excel(clinical_file) %>%
  mutate(
    Patient_ID_Clean = Patient_ID,
    event_str        = as.character(Event),
    event_binary     = case_when(
      event_str %in% c("1", "Relapse/Died", "RELAPSED") ~ 1,
      TRUE ~ 0
    )
  )

message(sprintf("Loaded clinical metadata for %d patients.", nrow(clinical_data_raw_updated)))

# ------------------------------------------------------------------------------
# 2. Ingest & De-multiplex PyClone Subclonal Outputs
# ------------------------------------------------------------------------------
message("Step 2: Parsing normalized PyClone output matrix and extracting gene symbols...")

if (!file.exists(pyclone_txt)) {
  stop("⚠️ Error: Missing PyClone output file at: ", pyclone_txt)
}

pyclone_results <- read.table(pyclone_txt, header = TRUE, sep = "\t")

fig2_ready_df <- pyclone_results %>%
  mutate(
    real_sample_id   = sapply(strsplit(as.character(mutation_id), "_"), `[`, 1),
    True_Gene_Symbol = sapply(strsplit(as.character(mutation_id), "_"), `[`, 2)
  ) %>%
  dplyr::select(mutation_id, sample_id = real_sample_id, hugo_symbol = True_Gene_Symbol, cluster_id, cellular_fraction)

pyclone_joint_df <- fig2_ready_df %>%
  inner_join(clinical_data_raw_updated, by = c("sample_id" = "Patient_ID_Clean"))

# ------------------------------------------------------------------------------
# 3. Figure 3A: Clonal Architecture Stratification Across Population (Violin Plot)
# ------------------------------------------------------------------------------
message("Step 3: Rendering Figure 3A - Population-level clonal architecture stratification...")

cluster_summary <- fig2_ready_df %>%
  group_by(cluster_id) %>%
  summarise(
    mutation_count = n(),
    mean_ccf       = mean(cellular_fraction),
    median_ccf     = median(cellular_fraction),
    .groups        = "drop"
  ) %>%
  arrange(desc(median_ccf))

plot_df_3a <- fig2_ready_df %>%
  mutate(cluster_id = factor(cluster_id, levels = c(7, 6, 5, 2, 1, 4, 3, 0)))

fig3a <- ggplot(plot_df_3a, aes(x = cluster_id, y = cellular_fraction, fill = cluster_id)) +
  geom_violin(trim = FALSE, alpha = 0.4, color = "black", scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8, color = "black", lwd = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8, aes(color = cluster_id)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "Clonal Architecture Stratification Across the B-ALL Population",
    subtitle = "Identification of a dominant truncal clone (Cluster 7) and 7 subclonal waves",
    x = "Inferred Population Evolutionary Cluster",
    y = "Cellular Clonal Fraction (CCF)"
  ) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5, face = "italic"),
    axis.title      = element_text(face = "bold", size = 10),
    axis.text       = element_text(color = "black", size = 9)
  )

# ------------------------------------------------------------------------------
# 4. Figure 3B: Clonal Niche Preference of Top Recurrent Drivers (Trunk vs Branch)
# ------------------------------------------------------------------------------
message("Step 4: Analyzing Trunk vs Branch niche preferences for top driver genes...")

top_prevalence_drivers <- fig2_ready_df %>%
  group_by(hugo_symbol) %>%
  summarise(Patient_Freq = n_distinct(sample_id), .groups = "drop") %>%
  slice_max(order_by = Patient_Freq, n = 20, with_ties = FALSE) %>%
  pull(hugo_symbol)

ordered_genes <- fig2_ready_df %>%
  filter(hugo_symbol %in% top_prevalence_drivers) %>%
  group_by(hugo_symbol) %>%
  summarise(patient_count = n_distinct(sample_id), .groups = "drop") %>%
  arrange(desc(patient_count)) %>%
  pull(hugo_symbol)

gene_niche_stats <- fig2_ready_df %>%
  filter(hugo_symbol %in% top_prevalence_drivers) %>%
  mutate(Hierarchy = ifelse(cluster_id == 7, "Truncal (Trunk)", "Subclonal (Branch)")) %>%
  group_by(hugo_symbol, Hierarchy) %>%
  summarise(patient_n = n_distinct(sample_id), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Hierarchy, values_from = patient_n, values_fill = 0) %>%
  mutate(
    total_patients = `Truncal (Trunk)` + `Subclonal (Branch)`,
    truncal_pct    = `Truncal (Trunk)` / total_patients * 100,
    branch_pct     = `Subclonal (Branch)` / total_patients * 100
  ) %>%
  rowwise() %>%
  mutate(binom_p = binom.test(`Truncal (Trunk)`, total_patients, p = 0.5)$p.value) %>%
  ungroup() %>%
  mutate(sig_label = if_else(binom_p < 0.05, "*", ""))

gene_niche_stats$hugo_symbol <- factor(gene_niche_stats$hugo_symbol, levels = rev(ordered_genes))

panel_b_niche <- fig2_ready_df %>%
  filter(hugo_symbol %in% top_prevalence_drivers) %>%
  mutate(Hierarchy = ifelse(cluster_id == 7, "Truncal (Trunk)", "Subclonal (Branch)")) %>%
  group_by(hugo_symbol, Hierarchy) %>%
  summarise(n_mut = n(), .groups = "drop") %>%
  group_by(hugo_symbol) %>%
  mutate(Percentage = n_mut / sum(n_mut) * 100) %>%
  ungroup()

panel_b_niche$hugo_symbol <- factor(panel_b_niche$hugo_symbol, levels = rev(ordered_genes))

fig3b <- ggplot(panel_b_niche, aes(x = hugo_symbol, y = Percentage, fill = Hierarchy)) +
  geom_bar(stat = "identity", width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(data = gene_niche_stats, aes(x = hugo_symbol, y = 105, label = sig_label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Truncal (Trunk)" = "#6BAED6", "Subclonal (Branch)" = "#FDAE6B")) +
  theme_classic() +
  coord_flip() +
  labs(
    title = "Clonal Niche Preference of Top Recurrent Drivers",
    x = "Top 20 Recurrent B-ALL Driver Genes",
    y = "Cluster Mapping Proportion (%)"
  ) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    plot.title  = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position = "top",
    legend.title    = element_blank()
  )

# ------------------------------------------------------------------------------
# 5. Figure 3C: Multi-Omics Synergy & Threshold Collapse (Hits Score 0-3 vs >=4)
# ------------------------------------------------------------------------------
message("Step 5: Testing multi-omics synergy and threshold collapse effect...")

patient_star_synergy <- pyclone_joint_df %>%
  group_by(Patient_ID) %>%
  summarise(
    Hits_Score = sum(c(
      any(CREBBP == "Yes" | CREBBP == 1, na.rm = TRUE),
      any(IRF2BPL == "Yes" | IRF2BPL == 1, na.rm = TRUE),
      any(EGR1 == "Yes" | EGR1 == 1, na.rm = TRUE),
      any(PAX5 == "Yes" | PAX5 == 1, na.rm = TRUE),
      any(NRAS == "Yes" | NRAS == 1, na.rm = TRUE),
      any(`14p13.3_Del` == "Yes" | `14p13.3_Del` == 1, na.rm = TRUE),
      any(`11q13.3_Amp` == "Yes" | `11q13.3_Amp` == 1, na.rm = TRUE),
      any(`17p11.2_Amp` == "Yes" | `17p11.2_Amp` == 1, na.rm = TRUE),
      any(`12q22.3_Amp` == "Yes" | `12q22.3_Amp` == 1, na.rm = TRUE)
    )),
    Active_Subclonal_Waves = n_distinct(cluster_id[cluster_id != 7]),
    .groups = "drop"
  )

panel_c_binary_df <- patient_star_synergy %>%
  mutate(
    Hits_Group = ifelse(Hits_Score >= 4, "High Co-occurrence\n(Score >= 4)", "Low/Med Co-occurrence\n(Score 0 - 3)"),
    Hits_Group = factor(Hits_Group, levels = c("Low/Med Co-occurrence\n(Score 0 - 3)", "High Co-occurrence\n(Score >= 4)"))
  )

binary_wilcox_p <- wilcox.test(Active_Subclonal_Waves ~ Hits_Group, data = panel_c_binary_df)$p.value
c_p_annotation  <- sprintf("Wilcoxon P = %.4f *", binary_wilcox_p)

fig3c <- ggplot(panel_c_binary_df, aes(x = Hits_Group, y = Active_Subclonal_Waves, fill = Hits_Group)) +
  geom_violin(alpha = 0.4, color = "black", width = 0.5, trim = TRUE) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.15, color = "black", linewidth = 0.6) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 2.2, color = "gray20") +
  scale_fill_manual(values = c("#A2D2FF", "#FFD23F")) +
  theme_classic() +
  labs(
    title = "Multi-Omics Threshold Collapse Drives Subclonal Wave Explosion",
    x = "Cumulative High-Risk Hits",
    y = "Number of Active Subclonal Waves Per Patient"
  ) +
  annotate("text", x = 1.5, y = max(panel_c_binary_df$Active_Subclonal_Waves) * 0.95,
           label = c_p_annotation, fontface = "italic", size = 4, color = "black") +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.title      = element_text(face = "bold", size = 10),
    axis.text       = element_text(color = "black", size = 9)
  )

# ------------------------------------------------------------------------------
# 6. Figure 3D: Extreme Case Study Analysis 
# ------------------------------------------------------------------------------
message("Step 6: Extracting extreme index patient case study ")

extreme_case_id <- "Here the case id was de-identified"

panel_d_extreme_case <- pyclone_joint_df %>%
  filter(Patient_ID == extreme_case_id) %>%
  dplyr::select(mutation_id, hugo_symbol, cluster_id, cellular_fraction, Event, EFS_Month) %>%
  arrange(desc(cellular_fraction))

message(sprintf("Found %d mutations for index case %s:", nrow(panel_d_extreme_case), extreme_case_id))
print(panel_d_extreme_case)

# ------------------------------------------------------------------------------
# 7. Save Figures and Export Multi-Sheet Excel Reports
# ------------------------------------------------------------------------------
message("Step 7: Exporting Figure 3 PDF panels and Excel statistical sheets...")

output_pdf <- file.path(output_dir, "Figure3_PyClone_Evolution_Main.pdf")
pdf(output_pdf, width = 12, height = 10)
print(ggpubr::ggarrange(fig3a, fig3b, fig3c, ncol = 2, nrow = 2, labels = c("A", "B", "C")))
dev.off()

# Compute cluster prevalence across outcome groups for Excel report
cluster_cols <- paste0("Cluster_", 0:7)
cluster_matrix <- fig2_ready_df %>%
  mutate(cluster_flag = 1) %>%
  distinct(sample_id, cluster_id, cluster_flag) %>%
  pivot_wider(names_from = cluster_id, values_from = cluster_flag, values_fill = 0, names_prefix = "Cluster_")

clinical_with_clusters <- clinical_data_raw_updated %>%
  left_join(cluster_matrix, by = c("Patient_ID" = "sample_id"))

cluster_prevalence <- map_df(cluster_cols, function(cl) {
  if (!cl %in% colnames(clinical_with_clusters)) return(NULL)
  
  df <- clinical_with_clusters %>%
    select(Event, !!sym(cl)) %>%
    filter(!is.na(Event)) %>%
    mutate(Event_Label = ifelse(Event == 1, "Relapsed/Deceased", "Event-free"))
  
  tab <- table(df$Event_Label, df[[cl]])
  if (ncol(tab) < 2) {
    tab <- cbind(tab, rep(0, nrow(tab)))
    colnames(tab) <- c("0", "1")
  }
  
  event_total     <- sum(tab["Relapsed/Deceased", ])
  eventfree_total <- sum(tab["Event-free", ])
  
  event_pct     <- ifelse(event_total > 0, tab["Relapsed/Deceased", "1"] / event_total * 100, 0)
  eventfree_pct <- ifelse(eventfree_total > 0, tab["Event-free", "1"] / eventfree_total * 100, 0)
  pval          <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA)
  
  tibble(Cluster_ID = cl, Event_pct = event_pct, EventFree_pct = eventfree_pct, p_value = pval)
}) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

# Save multi-sheet Excel report
excel_out <- file.path(output_dir, "fig3_Cluster_Prevalence_Results.xlsx")
write.xlsx(
  list(
    Cluster_Prevalence = cluster_prevalence,
    Gene_Niche_Stats   = gene_niche_stats,
    Synergy_Threshold  = patient_star_synergy,
    Extreme_Case_Study = panel_d_extreme_case
  ),
  file     = excel_out,
  rowNames = FALSE
)

message("🎉 PyClone downstream analysis completed successfully!")
message("  - PDF Panels : ", output_pdf)
message("  - Excel Data : ", excel_out)
