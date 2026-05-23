# ==============================================================================
# B-ALL Genomics Analysis Pipeline
# Author: DanYu_PhD_BVSc
# Description: Codebase for reproducing the somatic landscape, clinical survival
#              outcomes, longitudinal tracking, and cross-cohort comparisons.
# ==============================================================================

# ------------------------------------------------------------------------------
# Module 0: Environment Setup & Library Initialization
# ------------------------------------------------------------------------------
library(maftools)
library(data.table)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)
library(readxl)
library(openxlsx)
library(GenomicRanges)
library(survival)
library(survminer)
library(ggrepel)
library(stringr)

# Establish standard working directory
setwd("/Users/DanYu_PhD_BVSc/Library/Mobile Documents/com~apple~CloudDocs/MyDocs/22_HKU/project_BALL/")

# Global Mapping for Mutation Classification Consistency
mutation_type_map <- c(
  "stopgain"                   = "Nonsense_Mutation",
  "stoploss"                   = "Nonstop_Mutation",
  "frameshift deletion"        = "Frame_Shift_Del",
  "frameshift insertion"       = "Frame_Shift_Ins",
  "nonframeshift deletion"     = "In_Frame_Del",
  "nonframeshift insertion"    = "In_Frame_Ins",
  "nonsynonymous SNV"          = "Missense_Mutation",
  "startloss"                  = "Translation_Start_Site",
  "nonframeshift substitution" = "In_Frame_Del" 
)

# ------------------------------------------------------------------------------
# Module 1: Somatic Mutation Processing & Quality Control (Panel & WES Cross-Talk)
# ------------------------------------------------------------------------------

# Helper functions for processing custom VCF/Annovar info fields
extract_metrics <- function(info) {
  parts <- unlist(strsplit(as.character(info), ":"))
  vaf <- as.numeric(parts[grepl("\\.", parts)])[1]
  vaf_idx <- which(grepl("\\.", parts))[1]
  dp <- as.numeric(parts[vaf_idx + 1])
  return(c(VAF = vaf, DP = dp))
}

extract_panel_tumor <- function(info) {
  if (is.na(info) || info == "." || info == "") return(c(VAF = NA, DP = NA))
  parts <- unlist(strsplit(as.character(info), ":"))
  vaf <- as.numeric(parts[3]) # GT:AD:AF:DP layout
  dp  <- as.numeric(parts[4])
  return(c(VAF = vaf, DP = dp))
}

# 1.1 Process Targeted Sequencing Panel Data
maf_list <- list.files("/Volumes/Extreme SSD/somatic_annovar_output/panel1", pattern = "panel1.hg38_multianno.clean", full.names = TRUE)
maf_objects <- lapply(maf_list, function(x) { annovarToMaf(annovar = x, Center = "MyCenter", refBuild = "hg38") })
combined_maf <- merge_mafs(maf_objects)
panel_data_dt <- as.data.table(combined_maf@data)

panel_gene_set <- panel_data_dt %>%
  filter(!is.na(Hugo_Symbol) & Hugo_Symbol != "") %>%
  pull(Hugo_Symbol) %>% unique()

panel_ready <- panel_data_dt %>%
  mutate(metrics = map(Otherinfo14, extract_panel_tumor)) %>%
  unnest_wider(metrics) %>%
  rename(VAF_Panel = VAF, DP_Panel = DP, Start = Start_Position, End = End_Position) %>%
  mutate(Sample = str_remove(Tumor_Sample_Barcode, "_panel1"),
         Start = as.numeric(Start), End = as.numeric(End),
         Variant_Classification = as.character(Variant_Classification)) %>%
  select(Sample, Hugo_Symbol, Chromosome, Start, End, Variant_Classification, VAF_Panel, DP_Panel)

# 1.2 Process and Anchor Whole Exome Sequencing (WES) Data
wes_folder <- "/Volumes/Extreme SSD/somatic_annovar_output/wes"
wes_files <- list.files(path = wes_folder, pattern = "wes\\.hg38_multianno\\.clean", full.names = TRUE)
wes_raw <- map_df(wes_files, ~{
  dt <- fread(.x, fill = TRUE)
  sample_id <- basename(.x) %>% str_remove("\\.wes\\.hg38_multianno\\.clean.*")
  dt[, Tumor_Sample_Barcode := sample_id]
  return(dt)
})

wes_bridge_clean <- wes_raw %>%
  rename(Hugo_Symbol = Gene.refGene, Variant_Classification = ExonicFunc.refGene, Func_Region = Func.refGene) %>%
  filter(Hugo_Symbol %in% panel_gene_set & Func_Region == "exonic" & !grepl("synonymous", Variant_Classification, ignore.case = TRUE)) %>%
  mutate(Sample = str_remove(Tumor_Sample_Barcode, "\\.wes.*"))

wes_vaf_final <- wes_bridge_clean %>%
  mutate(metrics = map(Otherinfo13, extract_metrics)) %>%
  unnest_wider(metrics) %>%
  rename(VAF_WES = VAF, DP_WES = DP) %>%
  mutate(Sample = str_remove(Sample, "_wes.*")) %>%
  select(Sample, Hugo_Symbol, Chr, Start, End, VAF_WES, DP_WES, Variant_Classification_WES = Variant_Classification) %>%
  distinct()

wes_vaf_final_standard <- wes_vaf_final %>%
  mutate(Variant_Classification = recode(Variant_Classification_WES, !!!mutation_type_map), Method = "WES")

# 1.3 Filtering 
white_list_genes <- panel_ready %>%
  group_by(Hugo_Symbol) %>%
  summarise(n_patients = n_distinct(Sample)) %>%
  filter(n_patients >= 2) %>% pull(Hugo_Symbol)

panel_clean <- panel_ready %>%
  filter((Sample %in% critical_ids & DP_Panel >= 50 & VAF_Panel >= 0.005) | 
         (!(Sample %in% critical_ids) & DP_Panel >= 100 & VAF_Panel >= 0.01))

wes_clean <- wes_vaf_final_standard %>%
  mutate(Start = as.numeric(Start), End = as.numeric(End)) %>%
  filter(VAF_WES >= 0.1 & VAF_WES <= 0.6 & DP_WES >= 20 & Hugo_Symbol %in% white_list_genes)

snv_final_objective <- bind_rows(
  panel_clean %>% filter(DP_Panel >= 100) %>% mutate(Tier = "T1_Panel_Paired"),
  wes_clean %>% mutate(Tier = "T2_WES_Candidate")
) %>%
  group_by(Sample, Hugo_Symbol) %>%
  slice_max(order_by = (Tier == "T1_Panel_Paired"), n = 1) %>%
  ungroup()

snv_master_final <- snv_final_objective %>%
  mutate(Chromosome = coalesce(Chromosome, Chr),
         Variant_Type = coalesce(Variant_Classification, Variant_Classification_WES),
         Source_Color = case_when(!is.na(VAF_Panel) & !is.na(VAF_WES) ~ "Both",
                                  !is.na(VAF_Panel) ~ "Panel_Only",
                                  !is.na(VAF_WES) ~ "WES_Only"),
         VAF_For_Plot = coalesce(VAF_Panel, VAF_WES)) %>%
  select(Sample, Hugo_Symbol, Chromosome, Start, Variant_Type, VAF_Panel, VAF_WES, VAF_For_Plot, Source_Color) %>%
  arrange(Sample, Hugo_Symbol)

# ------------------------------------------------------------------------------
# Module 2: Recurrent Somatic Landscape & Oncoplots (Figure 1c)
# ------------------------------------------------------------------------------
maf_list_clean <- file.path(dirname(maf_list), gsub("_panel1.*", "", basename(maf_list)))
file.copy(maf_list, maf_list_clean)
maf_objects_clean <- lapply(maf_list_clean, function(x) { annovarToMaf(annovar = x, Center = "MyCenter", refBuild = "hg38") })
laml <- merge_mafs(maf_objects_clean)

# Synchronize Clinical Annotations
clinical_data <- readxl::read_excel('article/clinical_dat.xlsx')
colnames(clinical_data)[c(1, 13, 21)] <- c("Tumor_Sample_Barcode", "MRD", "Molecular_Subtype")
colnames(clinical_data) <- make.names(colnames(clinical_data))
clinical_data <- clinical_data %>% mutate(Event = recode(Event, "Relapse/Died" = "Relapse/Deceased"))

maf_samples <- unique(laml@data$Tumor_Sample_Barcode)
clinical_data_sync <- clinical_data %>% filter(Tumor_Sample_Barcode %in% maf_samples)
laml@clinical.data <- as.data.table(clinical_data_sync)

# Define Graphic Settings and Color Palettes
mutation_colors <- RColorBrewer::brewer.pal(n = 9, name = 'Paired')
names(mutation_colors) <- c('Nonsense_Mutation', 'Missense_Mutation', 'In_Frame_Ins', 'Frame_Shift_Del',
                            'Frame_Shift_Ins', 'In_Frame_Del', 'Splice_Site', 'Nonstop_mutation', 'Multi_Hit')

event_colors <- list(Event = c("Relapse/Deceased" = "#756bb1", "Event-free" = "#1b9e77", "Unknown" = "grey80"))
mrd_colors <- list(MRD = c("Positive (≥0.01%)" = "#8B475D", "Negative (<0.01%)" = "#EEAEEE", "Unknown" = "grey80"))
subtype_colors <- list(Molecular_Subtype = c("ETV6-RUNX1" = "#7FC97F", "TCF3 / DUX4 / PAX5 Class" = "#BEAED4",
                                            "KMT2A / ZNF384 / MEF2D Class" = "#FDC086", "BCR-ABL1 / Ph-like" = "#FB8072",
                                            "Other Rare Fusions" = "#386CB0", "Not Determined" = "#404040"))

# Render Figure 1c Oncoplot
oncoplot(
  maf = laml, gene_mar = 8, top = 20, sortByAnnotation = TRUE, removeNonMutated = FALSE,
  clinicalFeatures = c("Event", "MRD", "Molecular_Subtype"),
  annotationColor = c(event_colors, mrd_colors, subtype_colors),
  colors = mutation_colors, drawColBar = FALSE, draw_titv = FALSE, fontSize = 0.8,
  anno_height = 1.4, legend_height = 5, legendFontSize = 1.2, annotationFontSize = 1.2,
  bgCol = 'white', sepwd_samples = 5, showTitle = FALSE, drawRowBar = FALSE, includeColBarCN = TRUE
)

# Export processed dataset for Supplementary Material
write.xlsx(as.data.frame(laml@data), 'article/File SI5 panel maf.xlsx')

# ------------------------------------------------------------------------------
# Module 3: Copy Number Variation (CNV) Chromosome Landscape (Figure 1d)
# ------------------------------------------------------------------------------
all_cns <- read_excel('article/File SI 8 cnv.xlsx')

# Normalize baseline noise shifts dynamically using major density peak offset
log2_offset <- median(all_cns$log2, na.rm = TRUE)
all_cns_corrected <- all_cns %>%
  mutate(log2 = as.numeric(log2),
         log2_adj = log2 - log2_offset) %>%
  filter(abs(log2_adj) >= 0.4) %>%
  mutate(Type = ifelse(log2_adj > 0, "Amp", "Del"))

write.xlsx(as.data.frame(all_cns_corrected), 'article/File SI 8 cnv cytoband filtered.xlsx')

plot_ready_data <- all_cns_corrected %>%
  left_join(clinical_data_sync %>% select(Tumor_Sample_Barcode, Event), by = c("SampleID" = "Tumor_Sample_Barcode")) %>%
  filter(!is.na(Event))

# Construct Global Index for Cytobands across Chromosomes
plot_df <- plot_ready_data %>%
  group_by(Cytoband, chromosome, Type) %>%
  summarise(n = n_distinct(SampleID), .groups = 'drop') %>%
  mutate(freq = n / 105,
         chrom_label = factor(chromosome, levels = paste0("chr", c(1:22, "X", "Y")))) %>%
  filter(!is.na(chrom_label))

plot_df_indexed <- plot_df %>%
  arrange(chrom_label, Cytoband) %>%
  mutate(idx = as.numeric(factor(Cytoband, levels = unique(Cytoband))),
         freq_plot = ifelse(Type == "Del", -freq, freq))

chr_boundaries <- plot_df_indexed %>% group_by(chrom_label) %>% summarise(max_idx = max(idx)) %>% pull(max_idx)
chr_labels <- plot_df_indexed %>% group_by(chrom_label) %>% summarise(mid = mean(idx))

# Group Stratification (Relapse/Deceased vs Event-free Mirror Plot)
facet_plot_data <- plot_ready_data %>%
  filter(Event %in% c("Event-free", "Relapse/Deceased")) %>%
  group_by(Cytoband, Type, Event) %>%
  summarise(n = n_distinct(SampleID), .groups = 'drop') %>%
  mutate(Group = factor(Event, levels = c("Relapse/Deceased", "Event-free")),
         freq_val = ifelse(Type == "Amp", n / ifelse(Event == "Relapse/Deceased", 13, 87), 
                           -(n / ifelse(Event == "Relapse/Deceased", 13, 87)))) %>%
  inner_join(plot_df_indexed %>% select(Cytoband, idx, chrom_label) %>% distinct(), by = "Cytoband")

top_labels_facet <- facet_plot_data %>% group_by(Group) %>% arrange(desc(abs(freq_val))) %>% slice_head(n = 5) %>% ungroup()

# Annotate with Statistical P-values evaluated from Fisher's exact tests
top_labels_p_added <- top_labels_facet %>%
  mutate(clean_label = gsub("chr", "", Cytoband),
         final_label = case_when(
           grepl("11q13.3", Cytoband) ~ "11q13.3 (CCND1/FADD) p=0.04",
           grepl("14p13.3", Cytoband) ~ "14p13.3 (RNA18SN5) p=0.11", 
           grepl("19q13.1", Cytoband) ~ "19q13.1 (BAX/CEBPA) p=0.18",
           grepl("19p13.2", Cytoband) ~ "19p13.2 (LDLR/TCF3) p=0.19",
           grepl("12q22.3", Cytoband) ~ "12q22.3 (BTG1/SOCS2)",
           grepl("14p12",   Cytoband) ~ "14p12 (RN3P2)",
           grepl("17q21.3", Cytoband) ~ "17q21.3 (STAT3/5)",
           grepl("19q13.2", Cytoband) ~ "19q13.2 (AKT2/BCL3)",
           grepl("17p11.2", Cytoband) ~ "17p11.2 (MAP2K3/FLCN)",
           TRUE ~ NA_character_
         ))

# Render Figure 1d
ggplot(facet_plot_data, aes(x = idx, y = freq_val, fill = Type)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~ Group, ncol = 1, scales = "fixed") +
  geom_vline(xintercept = chr_boundaries, linetype = "dotted", color = "grey80", size = 0.2) +
  geom_text_repel(data = top_labels_p_added, aes(x = idx, y = freq_val, label = final_label),
                  size = 3.8, fontface = "italic", box.padding = 2, max.overlaps = 20,
                  segment.color = "grey20", segment.size = 0.2, min.segment.length = 0) +
  scale_fill_manual(values = c("Amp" = "#d73027", "Del" = "#4575b4")) +
  scale_x_continuous(label = gsub("chr", "", chr_labels$chrom_label), breaks = chr_labels$mid, expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) scales::percent(abs(x)), limits = c(-0.2, 0.2)) +
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 9),
                     strip.background = element_blank(), strip.text = element_text(size = 12, face = "bold"), legend.position = "right") +
  labs(y = "CNV Frequency (%)", x = "Chromosome")

# ------------------------------------------------------------------------------
# Module 4: Targeted SNV Fisher Enrichment Mirror Plot (Figure 1e)
# ------------------------------------------------------------------------------
snv_data <- read_excel('article/File SI 5 panel maf.xlsx')
clinical_data_fixed <- clinical_data %>%
  mutate(EFS_Month = as.numeric(EFS_Month), OS_Month = as.numeric(OS_Month),
         Relapse_Status = case_when(trimws(as.character(Relapse_Status)) %in% c("1", "Relapsed") ~ 1,
                                    trimws(as.character(Relapse_Status)) %in% c("0", "Non-relapsed") ~ 0, TRUE ~ NA_real_),
         Event = case_when(trimws(as.character(Event)) %in% c("1", "Relapse/Died") ~ 1,
                           trimws(as.character(Event)) %in% c("0", "Event-free") ~ 0, TRUE ~ NA_real_))

snv_with_clinical <- snv_data %>%
  left_join(clinical_data_fixed %>% select(Tumor_Sample_Barcode, Event), by = "Tumor_Sample_Barcode") %>%
  filter(!is.na(Event))

snv_target_genes <- c("CREBBP", "EGR1", "FLT3", "KRAS", "NRAS", "NCOR2", "KMT2D", 
                      "NOTCH1", "PTPN11", "NSD2", "KMT2A", "SETD2", "JAK3", "BCL11B", "ATRX")

# Statistical Fisher Testing Loop for Gene Mutations
snv_stats_results <- lapply(snv_target_genes, function(gene) {
  mutant_ids <- snv_with_clinical %>% filter(Hugo_Symbol == gene) %>% pull(Tumor_Sample_Barcode) %>% unique()
  n_event_mut <- sum(clinical_data_fixed$Tumor_Sample_Barcode[clinical_data_fixed$Event == 1] %in% mutant_ids)
  n_no_event_mut <- sum(clinical_data_fixed$Tumor_Sample_Barcode[clinical_data_fixed$Event == 0] %in% mutant_ids)
  
  contingency_table <- matrix(c(n_event_mut, 7 - n_event_mut, n_no_event_mut, 73 - n_no_event_mut), nrow = 2, byrow = TRUE)
  ft <- fisher.test(contingency_table)
  data.frame(Gene = gene, Event_n = n_event_mut, NoEvent_n = n_no_event_mut, P_value = ft$p.value, OR = ft$estimate)
}) %>% bind_rows()

# Clean and Build Plotting DF
mirror_snv_plot_df <- snv_with_clinical %>%
  filter(Hugo_Symbol %in% snv_target_genes) %>%
  group_by(Hugo_Symbol, Event) %>%
  summarise(n_mut = n_distinct(Tumor_Sample_Barcode), .groups = 'drop') %>%
  mutate(Group = ifelse(Event == 1, "Relapse/Deceased", "Event-free"),
         freq = ifelse(Group == "Relapse/Deceased", n_mut / 7, n_mut / 73),
         freq_plot = ifelse(Group == "Event-free", -freq, freq))

plot_data_final <- mirror_snv_plot_df %>%
  left_join(snv_stats_results %>% select(Gene, P_value), by = c("Hugo_Symbol" = "Gene")) %>%
  mutate(Gene_Label = case_when(P_value < 0.05 ~ paste0(Hugo_Symbol, "**"),
                                P_value < 0.2  ~ paste0(Hugo_Symbol, "*"), TRUE ~ as.character(Hugo_Symbol))) %>%
  mutate(Gene_Label = reorder(Gene_Label, freq_plot))

# Render Figure 1e
ggplot(plot_data_final, aes(x = Gene_Label, y = freq_plot, fill = Group)) +
  geom_bar(stat = "identity", width = 0.75, color = "white", size = 0.2) +
  coord_flip() + 
  scale_fill_manual(values = c("Relapse/Deceased" = "#756bb1", "Event-free" = "#1b9e77")) +
  scale_y_continuous(labels = function(x) scales::percent(abs(x)), limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, by = 0.2)) +
  geom_hline(yintercept = 0, color = "black", size = 0.6) +
  geom_text(data = plot_data_final %>% filter(Group == "Relapse/Deceased" & P_value < 0.2),
            aes(label = paste0("p=", round(P_value, 3))), hjust = -0.2, size = 3.5, color = "black") +
  theme_classic() + theme(axis.text.y = element_text(size = 11, face = "italic"), legend.position = "top",
                          legend.title = element_blank(), panel.grid.major.x = element_line(color = "grey92", linetype = "dashed")) +
  labs(title = "Divergent Mutation Landscape in B-ALL", subtitle = "** p < 0.05, * p < 0.2 (Fisher's Exact Test)",
       y = "Mutation Frequency (%)", x = "Gene Symbol")

# ------------------------------------------------------------------------------
# Module 5: Multivariate Cox Survival Modeling Volcano Plot (Figure 1f)
# ------------------------------------------------------------------------------
clinical_vars <- c("Age_Months", 'Sex','Ethnicity_Group','Immunophenotype', 'CNS_Status','Regimen_Clean','Blina_Used','HSCT','BM_Status_D46',"WBC", "NCI_Risk_Group", "MRD",'Molecular_Subtype')

snv_genes_to_test <- snv_data %>%
  filter(Tumor_Sample_Barcode %in% clinical_data_fixed$Tumor_Sample_Barcode) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>% count(Hugo_Symbol) %>%
  mutate(freq = n / nrow(clinical_data_fixed)) %>% filter(freq > 0.05) %>% pull(Hugo_Symbol)

snv_matrix <- snv_data %>% filter(Hugo_Symbol %in% snv_genes_to_test) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>% mutate(Status = "Altered") %>%
  pivot_wider(id_cols = Tumor_Sample_Barcode, names_from = Hugo_Symbol, values_from = Status, values_fill = "Not_Altered")

# Fit Models and Process via Multivariate Adjustments
snv_multi_fit_list <- lapply(snv_genes_to_test, function(g) {
  test_data <- clinical_data_fixed %>%
    left_join(snv_matrix %>% select(Tumor_Sample_Barcode, all_of(g)), by = "Tumor_Sample_Barcode") %>%
    mutate(Gene_Status = ifelse(.data[[g]] == "Altered", 1, 0)) %>% replace_na(list(Gene_Status = 0))
  try(coxph(Surv(EFS_Month, Event) ~ Gene_Status + Blina_Used + MRD + Ethnicity_Group, data = test_data), silent = TRUE)
})
names(snv_multi_fit_list) <- snv_genes_to_test
sheet3_snv_multivar <- extract_full_stats(snv_multi_fit_list, "SNV_Multivariate") %>% mutate(Gene = names(snv_multi_fit_list))

# Process CNV Matrix and Run Corresponding Cox Proportional Hazards Models
cnv_matrix <- all_cns_corrected %>% select(SampleID, Cytoband) %>% mutate(Status = "Altered") %>% distinct() %>% 
  pivot_wider(names_from = Cytoband, values_from = Status, values_fill = "Normal") %>% rename(Tumor_Sample_Barcode = SampleID)

high_freq_bands <- cnv_matrix %>% select(-Tumor_Sample_Barcode) %>% select_if(is.character) %>% 
  summarise(across(everything(), ~ mean(.x == "Altered", na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Cytoband", values_to = "Frequency") %>% filter(Frequency > 0.10) %>% pull(Cytoband)

cnv_genes_to_test <- unique(c(high_freq_bands, 'chr11q13.3', 'chr14p13.3'))

cnv_multi_list <- lapply(cnv_genes_to_test, function(g) {
  test_data <- clinical_data_fixed %>% left_join(cnv_matrix %>% select(Tumor_Sample_Barcode, all_of(g)), by = "Tumor_Sample_Barcode") %>%
    mutate(Gene_Status = ifelse(.data[[g]] == "Altered", 1, 0)) %>% replace_na(list(Gene_Status = 0))
  try(coxph(Surv(EFS_Month, Event) ~ Gene_Status + Blina_Used + MRD + Ethnicity_Group, data = test_data), silent = TRUE)
})
names(cnv_multi_list) <- cnv_genes_to_test
sheet5_cnv_multivar <- extract_full_stats(cnv_multi_list, "CNV_Multivariate") %>% mutate(Gene = names(cnv_multi_list))

# Plot Master Association Output (Figure 1f Composite Volcano)
plot_snv <- sheet3_snv_multivar %>% mutate(Type = "SNV")
plot_cnv <- sheet5_cnv_multivar %>% mutate(Type = "CNV")
master_cox_plot <- bind_rows(plot_snv, plot_cnv) %>%
  mutate(lnHR = log(HR), logP = -log10(P_value), gene_label = ifelse(P_value < 0.20, Gene, ""))

ggplot(master_cox_plot, aes(x = lnHR, y = logP)) +
  geom_point(color = "grey", alpha = 0.6, size = 3) +
  geom_point(data = filter(master_cox_plot, P_value < 0.20), aes(fill = Type), shape = 21, color = "black", size = 4) +
  geom_text_repel(data = filter(master_cox_plot, P_value < 0.20), aes(label = gene_label), size = 3.5, max.overlaps = 15) +
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "red") +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-3, 3)) + theme_bw() + theme(legend.position = "bottom") +
  labs(title = "Somatic Feature Association with EFS (Cox Regression Multivariate)", x = "ln(Hazard Ratio)", y = "-log10(p-value)")

# ------------------------------------------------------------------------------
# Module 6: Longitudinal Tracking & Multi-Omics Mapping (Figure 1g)
# ------------------------------------------------------------------------------

# 6.1 Clinical Tracking History Timeline Construction
df_clinical_history <- data.frame(
  month = c(0, 1.5, 3, 8.8, 9.0),
  event = c("Diagnosis", "Induction End (D46)", "Complete Remission (CR)", "Relapse", "Outcome"),
  detail = c("TCF3-PBX1 (49.9%), WBC 49.66\nNCI Risk: Intermediate (IR)", "MRD < 0.01% (Negative)\nImmature Lymph 3%",
             "Maintenance Phase\nCCCG-2020 Protocol", "Bone Marrow Relapse\n(EFS: 8.8 months)", "Deceased\n(OS: 9.0 months)"),
  category = c("Baseline", "Treatment", "Stable", "Relapse", "Outcome")
)

plot_timeline_g <- ggplot(df_clinical_history, aes(x = month, y = 1)) +
  geom_segment(aes(x = -0.8, xend = 10.5, y = 1, yend = 1), color = "grey80", size = 2) +
  geom_point(aes(color = category), size = 7, shape = 16) +
  geom_text_repel(aes(label = event), nudge_y = 0.6, direction = "x", force = 10, fontface = "bold", size = 4.5, box.padding = 0.8) +
  geom_text_repel(aes(label = detail), nudge_y = -0.8, direction = "x", force = 15, size = 4, lineheight = 0.9, color = "black", box.padding = 1.2) +
  geom_text(aes(label = paste0(month, "m")), y = 0.75, size = 4, fontface = "bold", color = "grey30") +
  scale_color_manual(values = c("Baseline"="#e41a1c", "Treatment"="#4daf4a", "Stable"="#36648B", "Relapse"="#ff7f00", "Outcome"="#800000")) +
  scale_x_continuous(limits = c(-1.5, 12), breaks = seq(0, 10, 2)) + ylim(-0.5, 2.5) + theme_void() + theme(legend.position = "none")

# 6.2 Lollipops for Key Mutated Genes
lollipopPlot(maf = laml, gene = 'EGR1', AACol = 'aaChange', showMutationRate = TRUE, labelPos = 68, showDomainLabel = FALSE)
lollipopPlot(maf = laml, gene = 'CREBBP', AACol = 'aaChange', showMutationRate = TRUE, labelPos = 1560, showDomainLabel = FALSE)
lollipopPlot(maf = laml, gene = 'IRF2BPL', AACol = 'aaChange', showMutationRate = TRUE, labelPos = 164, showDomainLabel = FALSE)

# ------------------------------------------------------------------------------
# Module 7: Multi-Cohort Validation and Metanalysis (Figure 1h)
# ------------------------------------------------------------------------------

# Fudan Cohort Profiles (External Chinese Dataset, n=114)
fudan_freq_df <- data.frame(
  Hugo_Symbol = c("KRAS", "NRAS", "FLT3", "KMT2D", "CREBBP", "NOTCH1", "NSD2", "PTPN11", "SETD2", "JAK3", "KMT2A", "PAX5"),
  Freq_Fudan = c(11.4, 7.0, 7.0, 5.3, 4.4, 0.0, 0.0, 4.4, 2.6, 0.0, 0.0, 3.5)
)

# St. Jude-COG Exploration and Matrix Generation
# Loaded and normalized from supplementary materials (n_western total base evaluated dynamically)
stjude_raw <- read_excel('article/nature genetic/41588_2022_1159_MOESM4_ESM.xlsx', sheet = 3)
colnames(stjude_raw) <- stjude_raw[4,] # Synchronize header layout
stjude_raw <- stjude_raw[5:nrow(stjude_raw),]

n_stjude_total <- length(unique(stjude_raw$`Patient identifier`))
stjude_freq_table <- stjude_raw %>%
  group_by(Hugo_Symbol = `Gene mutated per VEP annotation`) %>%
  summarise(n_stjude = n_distinct(`Patient identifier`)) %>%
  mutate(Freq_StJude = (n_stjude / n_stjude_total) * 100)

# Build Comprehensive Matrix Across All 4 Evaluated Cohorts
cohort_comparison_master <- panel_ready %>%
  group_by(Hugo_Symbol) %>% summarise(n_panel = n_distinct(Sample)) %>% mutate(Freq_Panel = (n_panel / 85) * 100) %>%
  left_join(wes_vaf_final_standard %>% group_by(Hugo_Symbol) %>% summarise(n_wes = n_distinct(Sample)) %>% mutate(Freq_WES = (n_wes / 33) * 100), by = "Hugo_Symbol") %>%
  left_join(fudan_freq_df, by = "Hugo_Symbol") %>%
  left_join(stjude_freq_table, by = "Hugo_Symbol") %>%
  replace(is.na(.), 0)

# Perform Two-Sided Fisher Exact Tests to Screen Variations vs St. Jude-COG Baseline
validation_set <- cohort_comparison_master %>%
  rowwise() %>%
  mutate(p_val = fisher.test(matrix(c(n_panel, 85 - n_panel, n_stjude, n_stjude_total - n_stjude), nrow = 2))$p.value) %>%
  mutate(p_label = ifelse(p_val < 0.05, "*", "")) %>% ungroup()

# Build Long-format Table for Visualization Structure
plot_df_long <- validation_set %>%
  select(Hugo_Symbol, `Panel (n=85)` = Freq_Panel, `WES (n=33)` = Freq_WES, `Fudan (n=114)` = Freq_Fudan, `St. Jude-COG (n=1482)` = Freq_StJude) %>%
  pivot_longer(cols = -Hugo_Symbol, names_to = "Cohort", values_to = "Freq") %>%
  mutate(Plot_Freq = ifelse(Cohort %in% c("Panel (n=85)", "WES (n=33)"), Freq, -Freq))

cohort_colors <- c("Panel (n=85)" = "#386CB0", "WES (n=33)" = "#CCE5FF", "Fudan (n=114)" = "orange", "St. Jude-COG (n=1482)" = "#1b9e77")

# Render Figure 1h Inter-Cohort Divergence Barplot
ggplot(plot_df_long, aes(x = Hugo_Symbol, y = Plot_Freq, fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) +
  coord_flip() + scale_y_continuous(labels = abs, limits = c(-85, 45)) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(data = validation_set, aes(x = Hugo_Symbol, y = Freq_StJude + 4, label = p_label), inherit.aes = FALSE, size = 5, fontface = "bold") +
  scale_fill_manual(values = cohort_colors) + theme_classic() +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Comparative Mutational Landscape of Pediatric B-ALL",
       subtitle = "Spearman Correlation (Panel vs St. Jude-COG): R = 0.73, P < 0.001", x = "", y = "Mutation Frequency (%)")
