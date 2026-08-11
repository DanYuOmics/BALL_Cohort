#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_cox_high_risk_alterations.R
# Description: Figure 2D & Supplementary Tables - Univariate and multivariate 
#              Cox proportional hazards regression for clinical factors, SNVs, 
#              and CNVs associated with Event-Free Survival (EFS). Exports multi-sheet 
#              Excel reports and generates combined association volcano plots.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(survival)
  library(survminer)
  library(openxlsx)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./article/clinical_dat.xlsx")
maf_file      <- Sys.getenv("MAF_FILE", unset = "./article/File SI 5 panel maf.xlsx")
cnv_file      <- Sys.getenv("CNV_FILE", unset = "./article/File SI 8 cnv cytoband filtered.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Helper Function: Extract Comprehensive Cox Regression Statistics
# ------------------------------------------------------------------------------
extract_full_stats <- function(fit_list, type_label = "SNV") {
  purrr::map_dfr(names(fit_list), function(gene_name) {
    x <- fit_list[[gene_name]]
    if (is.null(x) || inherits(x, "try-error")) return(NULL)
    s <- summary(x)
    
    # Extract the first variable coefficient (primary feature)
    coef_mat <- s$coefficients[1, , drop = FALSE]
    conf_mat <- s$conf.int[1, , drop = FALSE]

    data.frame(
      Gene          = gene_name,
      Variable      = rownames(coef_mat),
      HR            = conf_mat[1, "exp(coef)"],
      Lower_95      = conf_mat[1, "lower .95"],
      Upper_95      = conf_mat[1, "upper .95"],
      Std_Error     = coef_mat[1, "se(coef)"],
      z_value       = coef_mat[1, "z"],
      P_value       = coef_mat[1, "Pr(>|z|)"],
      N_total       = s$n,
      Events        = s$nevent,
      AIC           = AIC(x),
      Concordance   = s$concordance[1],
      Analysis_Type = type_label,
      stringsAsFactors = FALSE
    )
  })
}

# ------------------------------------------------------------------------------
# 3. Load and Preprocess Clinical Metadata
# ------------------------------------------------------------------------------
message("Step 1: Loading and harmonizing clinical metadata...")
if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical dataset at: ", clinical_file)
}

clinical_data <- read_excel(clinical_file)

clinical_data_fixed <- clinical_data %>%
  mutate(
    EFS_Month = as.numeric(EFS_Month),
    OS_Month  = as.numeric(OS_Month),
    Relapse_Status = case_when(
      trimws(as.character(Relapse_Status)) %in% c("1", "Relapsed", "RELAPSED") ~ 1,
      trimws(as.character(Relapse_Status)) %in% c("0", "Non-relapsed", "NON-RELAPSED") ~ 0,
      TRUE ~ NA_real_
    ),
    Event = case_when(
      trimws(as.character(Event)) %in% c("1", "Relapse/Died") ~ 1,
      trimws(as.character(Event)) %in% c("0", "Event-free") ~ 0,
      TRUE ~ NA_real_
    )
  )

# Standardize key column names
colnames(clinical_data_fixed)[1]  <- "Tumor_Sample_Barcode"
colnames(clinical_data_fixed)[21] <- "Molecular_Subtype"

# ------------------------------------------------------------------------------
# 4. Univariate Cox Analysis for Clinical Factors
# ------------------------------------------------------------------------------
message("Step 2: Running Univariate Cox regression for clinical factors...")
clinical_vars <- c("Age_Months", "Sex", "Ethnicity_Group", "Immunophenotype", 
                   "CNS_Status", "Regimen_Clean", "Blina_Used", "HSCT", 
                   "BM_Status_D46", "WBC", "NCI_Risk_Group", "MRD", "Molecular_Subtype")

clinical_fit_list <- lapply(clinical_vars, function(x) {
  formula_str <- paste0("Surv(EFS_Month, Event) ~ ", x)
  try(coxph(as.formula(formula_str), data = clinical_data_fixed), silent = TRUE)
})
names(clinical_fit_list) <- clinical_vars

sheet1_clinical_univ <- extract_full_stats(clinical_fit_list, "Clinical_Univariate")

# ------------------------------------------------------------------------------
# 5. SNV Cox Regression Analysis (Univariate & Multivariate)
# ------------------------------------------------------------------------------
message("Step 3: Constructing SNV matrix and running Cox models...")
if (!file.exists(maf_file)) {
  stop("⚠️ Error: Missing MAF dataset at: ", maf_file)
}

maf_data <- read_excel(maf_file)

# Select SNVs with alteration frequency > 5%
snv_genes_to_test <- maf_data %>%
  filter(Tumor_Sample_Barcode %in% clinical_data_fixed$Tumor_Sample_Barcode) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  count(Hugo_Symbol) %>%
  mutate(freq = n / nrow(clinical_data_fixed)) %>%
  filter(freq > 0.05) %>%
  pull(Hugo_Symbol)

snv_matrix <- maf_data %>%
  filter(Hugo_Symbol %in% snv_genes_to_test) %>%
  distinct(Hugo_Symbol, Tumor_Sample_Barcode) %>%
  mutate(Status = "Altered") %>%
  pivot_wider(
    id_cols     = Tumor_Sample_Barcode,
    names_from  = Hugo_Symbol,
    values_from = Status,
    values_fill = "Not_Altered"
  )

# SNV Univariate Cox
snv_uni_fit_list <- lapply(snv_genes_to_test, function(g) {
  test_data <- clinical_data_fixed %>%
    left_join(snv_matrix %>% dplyr::select(Tumor_Sample_Barcode, all_of(g)), by = "Tumor_Sample_Barcode") %>%
    mutate(Gene_Status = ifelse(.data[[g]] == "Altered", 1, 0)) %>%
    replace_na(list(Gene_Status = 0))
  
  try(coxph(Surv(EFS_Month, Event) ~ Gene_Status, data = test_data), silent = TRUE)
})
names(snv_uni_fit_list) <- snv_genes_to_test
sheet2_snv_univ <- extract_full_stats(snv_uni_fit_list, "SNV_Univariate")

# SNV Multivariate Cox (Adjusting for Blina_Used, MRD, and Ethnicity_Group)
snv_multi_fit_list <- lapply(snv_genes_to_test, function(g) {
  test_data <- clinical_data_fixed %>%
    left_join(snv_matrix %>% dplyr::select(Tumor_Sample_Barcode, all_of(g)), by = "Tumor_Sample_Barcode") %>%
    mutate(Gene_Status = ifelse(.data[[g]] == "Altered", 1, 0)) %>%
    replace_na(list(Gene_Status = 0))
  
  try(coxph(Surv(EFS_Month, Event) ~ Gene_Status + Blina_Used + MRD + Ethnicity_Group, data = test_data), silent = TRUE)
})
names(snv_multi_fit_list) <- snv_genes_to_test
sheet3_snv_multivar <- extract_full_stats(snv_multi_fit_list, "SNV_Multivariate")

# ------------------------------------------------------------------------------
# 6. CNV Cox Regression Analysis (Univariate & Multivariate)
# ------------------------------------------------------------------------------
message("Step 4: Constructing CNV matrix and running Cox models...")
if (!file.exists(cnv_file)) {
  stop("⚠️ Error: Missing CNV dataset at: ", cnv_file)
}

cnv_band <- read_excel(cnv_file)

cnv_matrix <- cnv_band %>%
  dplyr::select(SampleID, Cytoband, Type) %>%
  dplyr::mutate(Status = "Altered") %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from  = Cytoband, 
    values_from = Status, 
    values_fill = "Normal"
  ) %>%
  dplyr::rename(Tumor_Sample_Barcode = SampleID)

# Select high-frequency Cytobands (> 10% frequency) + target bands
band_frequencies <- cnv_matrix %>%
  dplyr::select(-Tumor_Sample_Barcode) %>%
  dplyr::select_if(is.character) %>%
  summarise(across(everything(), ~ mean(.x == "Altered", na.rm = TRUE))) %>%
  tidyr::pivot_longer(everything(), names_to = "Cytoband", values_to = "Frequency")

high_freq_bands <- band_frequencies %>%
  dplyr::filter(Frequency > 0.10) %>%
  dplyr::pull(Cytoband)

cnv_genes_to_test <- unique(c(high_freq_bands, "chr11q13.3", "chr14p13.3"))

# CNV Univariate Cox
cnv_uni_list <- lapply(cnv_genes_to_test, function(g) {
  test_data <- clinical_data_fixed %>%
    dplyr::left_join(cnv_matrix %>% dplyr::select(Tumor_Sample_Barcode, all_of(g)), by = "Tumor_Sample_Barcode") %>%
    mutate(Gene_Status = ifelse(.data[[g]] == "Altered", 1, 0)) %>%
    replace_na(list(Gene_Status = 0))
  
  try(coxph(Surv(EFS_Month, Event) ~ Gene_Status, data = test_data), silent = TRUE)
})
names(cnv_uni_list) <- cnv_genes_to_test
sheet4_cnv_univ <- extract_full_stats(cnv_uni_list, "CNV_Univariate")

# CNV Multivariate Cox
cnv_multi_list <- lapply(cnv_genes_to_test, function(g) {
  test_data <- clinical_data_fixed %>%
    dplyr::left_join(cnv_matrix %>% dplyr::select(Tumor_Sample_Barcode, all_of(g)), by = "Tumor_Sample_Barcode") %>%
    mutate(Gene_Status = ifelse(.data[[g]] == "Altered", 1, 0)) %>%
    replace_na(list(Gene_Status = 0))
  
  try(coxph(Surv(EFS_Month, Event) ~ Gene_Status + Blina_Used + MRD + Ethnicity_Group, data = test_data), silent = TRUE)
})
names(cnv_multi_list) <- cnv_genes_to_test
sheet5_cnv_multivar <- extract_full_stats(cnv_multi_list, "CNV_Multivariate")

# ------------------------------------------------------------------------------
# 7. Calculate Alteration Counts & Export Consolidated Excel Report
# ------------------------------------------------------------------------------
message("Step 5: Exporting consolidated Cox regression results to Excel with FDR adjustments...")

snv_counts_df <- snv_matrix %>%
  tidyr::pivot_longer(cols = -Tumor_Sample_Barcode, names_to = "Gene", values_to = "Status") %>%
  group_by(Gene) %>%
  summarise(n_altered = sum(Status == "Altered", na.rm = TRUE), .groups = "drop")

cnv_counts_df <- cnv_matrix %>%
  tidyr::pivot_longer(cols = -Tumor_Sample_Barcode, names_to = "Gene", values_to = "Status") %>%
  group_by(Gene) %>%
  summarise(n_altered = sum(Status == "Altered", na.rm = TRUE), .groups = "drop")

sheet2_final <- sheet2_snv_univ %>% left_join(snv_counts_df, by = "Gene") %>% mutate(Type = "SNV", FDR_adj_P = p.adjust(P_value, method = "fdr"))
sheet3_final <- sheet3_snv_multivar %>% left_join(snv_counts_df, by = "Gene") %>% mutate(Type = "SNV", FDR_adj_P = p.adjust(P_value, method = "fdr"))
sheet4_final <- sheet4_cnv_univ %>% left_join(cnv_counts_df, by = "Gene") %>% mutate(Type = "CNV", FDR_adj_P = p.adjust(P_value, method = "fdr"))
sheet5_final <- sheet5_cnv_multivar %>% left_join(cnv_counts_df, by = "Gene") %>% mutate(Type = "CNV", FDR_adj_P = p.adjust(P_value, method = "fdr"))

excel_output <- file.path(output_dir, paste0("File_SI_8_COX_", Sys.Date(), ".xlsx"))
wb_final <- createWorkbook()

sheet_list <- list(
  "Clinical_Univ" = sheet1_clinical_univ,
  "SNV_Univ"      = sheet2_final,
  "SNV_Multivar"  = sheet3_final,
  "CNV_Univ"      = sheet4_final,
  "CNV_Multivar"  = sheet5_final
)

for (s_name in names(sheet_list)) {
  addWorksheet(wb_final, s_name)
  writeData(wb_final, s_name, sheet_list[[s_name]])
}

saveWorkbook(wb_final, excel_output, overwrite = TRUE)
message("✅ Successfully saved Cox summary Excel report to: ", excel_output)

# ------------------------------------------------------------------------------
# 8. Render Joint Volcano Plot (SNV + CNV High-Risk Alterations)
# ------------------------------------------------------------------------------
message("Step 6: Rendering joint volcano scatter plot (Figure 2D)...")

combined_plot_data <- bind_rows(sheet3_final, sheet5_final) %>%
  arrange(P_value) %>%
  mutate(
    lnHR = log(HR),
    logP = -log10(P_value),
    gene_label = ifelse(P_value < 0.20, Gene, "")
  )

p2d <- ggplot(combined_plot_data, aes(x = lnHR, y = logP)) +
  geom_point(aes(size = n_altered), color = "grey", alpha = 0.6) +
  geom_point(
    data = filter(combined_plot_data, P_value < 0.20),
    aes(size = n_altered), shape = 21, color = "black", fill = "#377eb8"
  ) +
  geom_text_repel(
    data = filter(combined_plot_data, P_value < 0.20),
    aes(label = gene_label), size = 3.5, max.overlaps = 15
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "red") +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  coord_cartesian(xlim = c(-3, 3)) +
  scale_size(name = "Altered Patients (n)") +
  labs(
    title = "Combined SNV & CNV Associations with EFS (Multi-variable Cox)",
    subtitle = "Solid red line = p < 0.05, Dashed line = p < 0.10",
    x = "ln(Hazard Ratio)",
    y = "-log10(p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

volcano_pdf <- file.path(output_dir, "Figure2D_combined_cox_volcano.pdf")
ggsave(volcano_pdf, plot = p2d, width = 8, height = 6)

message("✅ Successfully generated Figure 2D Joint Cox Volcano Plot: ", volcano_pdf)
