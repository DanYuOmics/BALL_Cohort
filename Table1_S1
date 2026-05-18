library(dplyr)
library(gtsummary)
library(flextable)
library(readxl)
library(labelled)
library(openxlsx)

# Set working directory
setwd("/Users/DanYu_PhD_BVSc/Library/Mobile Documents/com~apple~CloudDocs/MyDocs/22_HKU/project_BALL/")

# ==============================================================================
# 1. Data Import and Preprocessing
# ==============================================================================
clinical_data <- readxl::read_excel('article/SI 1 Clinical_Data_105_Patients.xlsx')

reset_gtsummary_theme()
my_stats_format <- list(
  all_continuous() ~ 1,
  all_categorical() ~ c(0, 1)
)

# Prepare dataset for Table 1 and stratified analyses
table1_prepared <- clinical_data %>%
  mutate(
    EFS_Month = as.numeric(EFS_Month),
    OS_Month = as.numeric(OS_Month),
    WBC = as.numeric(WBC),
    Event = case_when(
      Event %in% c("1", 1, "Relapsed", "RELAPSED") ~ "Relapse/Died",
      Event %in% c("0", 0, "Non-relapsed", "NON-RELAPSED") ~ "Event-free",
      TRUE ~ "Unknown"
    ),
    # Molecular subtype grouping
    Table_1_Group = case_when(
      Subtype %in% c("Not_Determined", "Negative", "Blank") | is.na(Subtype) ~ "Not Determined",
      Subtype == "ETV6-RUNX1" ~ "ETV6-RUNX1",
      Subtype %in% c("TCF3-PBX1", "DUX4-rearranged", "PAX5-rearranged") ~ "TCF3 / DUX4 / PAX5 Class",
      Subtype %in% c("KMT2A-rearranged", "ZNF384-rearranged", "MEF2D-rearranged") ~ "KMT2A / ZNF384 / MEF2D Class",
      Subtype %in% c("BCR-ABL1", "ABL1-rearranged", "CRLF2-rearranged") ~ "BCR-ABL1 / Ph-like",
      TRUE ~ "Other Rare Fusions"
    ),
    Table_1_Group = factor(Table_1_Group, levels = c(
      "ETV6-RUNX1", "TCF3 / DUX4 / PAX5 Class", "KMT2A / ZNF384 / MEF2D Class",
      "BCR-ABL1 / Ph-like", "Other Rare Fusions", "Not Determined"
    ))
  ) %>%
  set_variable_labels(
    Age_Months = "Age at Diagnosis (Months)",
    Sex = "Sex",
    Ethnicity_Group = "Ethnicity",
    WBC = "WBC at Diagnosis (x10^9/L)",
    Immunophenotype = "Immunophenotype",
    CNS_Status = "CNS Status",
    NCI_Risk_Group = "NCI Risk Group",
    Regimen_Clean = "Treatment",
    Blina_Used = "Blinatumomab Treatment",
    HSCT = "HSCT Receipt",
    BM_Status_D46 = "BM Status at D46",
    MRD_Flow = "Induction MRD (Flow)²",
    Relapse_Status = "Relapse Status",
    Vital_Status = "Vital Status",
    EFS_Month = "EFS (Months)",
    OS_Month = "OS (Months)",
    Table_1_Group = "Molecular Subtype"
  )

# ==============================================================================
# 2. Table 1: Overall Cohort Summary
# ==============================================================================
table1_obj <- table1_prepared %>%
  select(Age_Months, Sex, Ethnicity_Group, WBC, Immunophenotype, CNS_Status,
         NCI_Risk_Group, Table_1_Group, Regimen_Clean, Blina_Used, HSCT,
         BM_Status_D46, MRD_Flow, Relapse_Status, Vital_Status, EFS_Month, OS_Month) %>%
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{median} ({p25} - {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = my_stats_format,
    missing = "no"
  ) %>%
  modify_header(all_stat_cols() ~ "**Total Cohort** (N={n})") %>%
  modify_footnote(all_stat_cols() ~ "Statistics: median (IQR) or n (%). \n2 MRD was assessed via flow cytometry.") %>%
  bold_labels()

# ==============================================================================
# 3. Table S1: Stratified by Molecular Subtype
# ==============================================================================
table_s1_obj <- table1_prepared %>%
  filter(Table_1_Group != "Not Determined") %>%
  mutate(Table_1_Group = droplevels(Table_1_Group)) %>%
  tbl_summary(by = Table_1_Group, missing = "no", digits = my_stats_format) %>%
  add_p(test = list(all_categorical() ~ "fisher.test", all_continuous() ~ "kruskal.test"),
        test.args = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 2000))) %>%
  bold_labels()

# ==============================================================================
# 4. Tables S2–S4: Stratified by Relapse, Vital Status, and Event
# ==============================================================================
table_relapse_obj <- table1_prepared %>%
  tbl_summary(by = Relapse_Status, missing = "no", digits = my_stats_format) %>%
  add_p(test = list(all_categorical() ~ "fisher.test", all_continuous() ~ "kruskal.test"),
        test.args = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 2000))) %>%
  bold_labels()

table_died_obj <- table1_prepared %>%
  tbl_summary(by = Vital_Status, missing = "no", digits = my_stats_format) %>%
  add_p(test = list(all_categorical() ~ "fisher.test", all_continuous() ~ "kruskal.test"),
        test.args = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 2000))) %>%
  bold_labels()

table_event_obj <- table1_prepared %>%
  tbl_summary(by = Event, missing = "no", digits = my_stats_format) %>%
  add_p(test = list(all_categorical() ~ "fisher.test", all_continuous() ~ "kruskal.test"),
        test.args = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 2000))) %>%
  bold_labels()

# ==============================================================================
# 5. Tables S5–S7: Stratified by Gene Mutations (CREBBP, EGR1, BRAF)
# ==============================================================================
table_EGR1_Mut <- stratified_data %>%
  tbl_summary(by = EGR1_Mut, missing = "no", digits = my_stats_format) %>%
  add_p(test = list(all_categorical() ~ "fisher.test", all_continuous() ~ "kruskal.test"),
        test.args = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 2000))) %>%
  bold_labels()

table_CREBBP_Mut <- stratified_data %>%
  tbl_summary(by = CREBBP_Mut, missing = "no", digits = my_stats_format) %>%
  add_p(test = list(all_categorical() ~ "fisher.test", all_continuous() ~ "kruskal.test"),
        test.args = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 2000))) %>%
  bold_labels()

table_BRAF_Mut <- stratified_data %>%
  tbl_summary(by = BRAF_Mut, missing = "no", digits = my_stats_format) %>%
  add_p(test = list(all_categorical() ~ "fisher.test", all_continuous() ~ "kruskal.test"),
        test.args = list(all_categorical() ~ list(simulate.p.value = TRUE, B = 2000))) %>%
  bold_labels()

# ==============================================================================
# 6. Export All Tables
# ==============================================================================
apply_style <- function(gt_obj) {
  gt_obj %>% as_flex_table() %>%
    font(fontname = "Arial", part = "all") %>%
    fontsize(size = 8, part = "all") %>%
    autofit()
}

save_as_docx(
  "Table 1 Summary" = apply_style(table1_obj),
  "Table S1 Comparison" = apply_style(table_s1_obj),
  "Table S2 Relapse" = apply_style(table_relapse_obj),
  "Table S3 Died" = apply_style(table_died_obj),
  "Table S4 Event" = apply_style(table_event_obj),
  "Table S5 EGR1_Mut" = apply_style(table_EGR1_Mut),
  "Table S6 CREBBP_Mut" = apply_style(table_CREBBP_Mut),
  "Table S7 BRAF_Mut" = apply_style(table_B
