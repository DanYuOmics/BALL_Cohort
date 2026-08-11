#!/usr/bin/env Rscript
# ==============================================================================
# Script: 01_oncoplot_landscape.R
# Description: Figure 2A - Load ANNOVAR outputs into a MAF object and plot 
#              the mutational landscape (Oncoplot) with clinical annotations.
# ==============================================================================

suppressPackageStartupMessages({
  library(maftools)
  library(data.table)
  library(readxl)
  library(dplyr)
  library(RColorBrewer)
})

# 1. Directory and Path Configurations
annovar_dir   <- Sys.getenv("ANNOVAR_PANEL_DIR", unset = "/Volumes/Extreme SSD/somatic_annovar_output/panel1")
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./data/clinical_dat.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Batch Loading ANNOVAR Clean Output Files into MAF Objects
maf_list <- list.files(annovar_dir, pattern = "panel1.hg38_multianno.clean", full.names = TRUE)

if (length(maf_list) == 0) {
  stop("⚠️ Error: No ANNOVAR clean output files found in directory: ", annovar_dir)
}

message(sprintf("Loading and converting %d ANNOVAR files to MAF objects...", length(maf_list)))

maf_objects <- lapply(maf_list, function(x) {
  annovarToMaf(annovar = x, Center = "MyCenter", refBuild = "hg38")
})
combined_maf <- merge_mafs(maf_objects)
laml <- combined_maf

# 3. Load and Synchronize Clinical Metadata
if (file.exists(clinical_file)) {
  message("Loading clinical metadata from: ", clinical_file)
  clinical_data <- read_excel(clinical_file)
  
  # Standardize key column names
  colnames(clinical_data)[1]  <- "Tumor_Sample_Barcode"
  colnames(clinical_data)[13] <- "MRD"
  colnames(clinical_data)[21] <- "Molecular_Subtype"

  clinical_data_dt <- as.data.table(clinical_data)
  laml@clinical.data <- clinical_data_dt
} else {
  warning("⚠️ Clinical file not found. Proceeding without clinical annotations.")
}

# 4. Define Color Schemes for Mutation Types and Clinical Features
mutation_colors <- brewer.pal(n = 9, name = 'Paired')
names(mutation_colors) <- c(
  'Nonsense_Mutation', 'Missense_Mutation', 'In_Frame_Ins', 'Frame_Shift_Del',
  'Frame_Shift_Ins', 'In_Frame_Del', 'Splice_Site', 'Nonstop_mutation', 'Multi_Hit'
)

event.col <- list(Event = c(
  "Relapse/Deceased" = "#756bb1", 
  "Event-free"       = "#1b9e77", 
  "Unknown"          = "grey80"
))

MRD.col <- list(MRD = c(
  "Positive (≥0.01%)" = "#8B475D", 
  "Negative (<0.01%)" = "#EEAEEE", 
  "Unknown"          = "grey80"
))

molecular.col <- list(Molecular_Subtype = c(
  "ETV6-RUNX1"                   = "#7FC97F",
  "TCF3 / DUX4 / PAX5 Class"     = "#BEAED4",
  "KMT2A / ZNF384 / MEF2D Class" = "#FDC086",
  "BCR-ABL1 / Ph-like"           = "#FB8072",
  "Other Rare Fusions"           = "#386CB0",
  "Not Determined"               = "#404040"
))

annotation_colors <- c(event.col, MRD.col, molecular.col)

# 5. Render and Save Oncoplot
output_pdf <- file.path(output_dir, "Figure2A_oncoplot.pdf")
pdf(output_pdf, width = 12, height = 8)

oncoplot(
  maf = laml,
  top = 20,
  clinicalFeatures = c("Event", "MRD", "Molecular_Subtype"),
  sortByAnnotation = TRUE,
  annotationColor = annotation_colors,
  colors = mutation_colors,
  drawColBar = FALSE,
  draw_titv = FALSE,
  fontSize = 0.8,
  anno_height = 1.4,
  bgCol = 'white',
  sepwd_samples = 5,
  showTitle = FALSE
)

dev.off()
message("✅ Successfully generated Figure 2A Oncoplot: ", output_pdf)
