#!/usr/bin/env Rscript
# ==============================================================================
# Script: 02_cnv_cytoband_landscape.R
# Description: Figure 2C - Process raw CNVkit .cns files, map to cytobands,
#              perform baseline log2 ratio adjustment, load clinical metadata,
#              and plot the mirrored CNV landscape across clinical outcomes.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(GenomicRanges)
  library(scales)
  library(openxlsx)
})

# ------------------------------------------------------------------------------
# 1. Environment & Path Configurations
# ------------------------------------------------------------------------------
cns_dir       <- Sys.getenv("CNVKIT_CNS_DIR", unset = "/Volumes/Extreme SSD/cnvkit_panel")
clinical_file <- Sys.getenv("CLINICAL_FILE", unset = "./data/clinical_data_augmented_fixed_20260620_perfect.xlsx")
output_dir    <- Sys.getenv("OUTPUT_DIR", unset = "./results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. Raw CNS File Processing & Cytoband Mapping
# ------------------------------------------------------------------------------
message("Step 1: Reading raw CNVkit .cns files...")
cns_files <- list.files(cns_dir, pattern = "TUMOUR\\.sort\\.dedup\\.cns$", full.names = TRUE)

if (length(cns_files) == 0) {
  stop("⚠️ Error: No CNVkit .cns files found in directory: ", cns_dir)
}

# Merge all CNS segment files
all_cns <- map_df(cns_files, function(f) {
  read_tsv(f, col_types = cols(.default = "c", start = "i", end = "i", log2 = "d"), show_col_types = FALSE) %>%
    mutate(SampleID = gsub("_TUMOUR.*", "", basename(f)))
})

message(sprintf("Successfully loaded %d CNS segments across %d samples.", 
                nrow(all_cns), n_distinct(all_cns$SampleID)))

# Map CNS segments to Cytoband using GenomicRanges
message("Step 2: Mapping genomic coordinates to cytobands...")
cns_gr <- makeGRangesFromDataFrame(all_cns, seqnames.field = "chromosome")

# Note: Ensure 'band_gr' object is loaded or fetched from UCSC hg38 cytoband data
if (!exists("band_gr")) {
  message("Fetching hg38 cytoband data from UCSC...")
  session <- rtracklayer::browserSession("UCSC")
  rtracklayer::genome(session) <- "hg38"
  query <- rtracklayer::ucscTableQuery(session, table = "cytoBand")
  cyto_table <- rtracklayer::getTable(query)
  band_gr <- GRanges(
    seqnames = cyto_table$chrom,
    ranges = IRanges(start = cyto_table$chromStart + 1, end = cyto_table$chromEnd),
    band = cyto_table$name
  )
}

band_hits <- findOverlaps(cns_gr, band_gr)
all_cns$Cytoband <- "Unknown"
all_cns$Cytoband[queryHits(band_hits)] <- paste0(
  seqnames(band_gr)[subjectHits(band_hits)], 
  band_gr$band[subjectHits(band_hits)]
)

# Export raw mapped CNV table
raw_cnv_export <- file.path(output_dir, "File_SI_8_cnv.xlsx")
write.xlsx(as.data.frame(all_cns), raw_cnv_export)
message("Saved raw mapped CNV table to: ", raw_cnv_export)

# ------------------------------------------------------------------------------
# 3. Log2 Baseline Adjustment & Denoising Filter
# ------------------------------------------------------------------------------
message("Step 3: Calculating median log2 ratio offset for baseline normalization...")
all_cns$log2 <- as.numeric(all_cns$log2)
log2_offset <- median(all_cns$log2, na.rm = TRUE)
message(sprintf("Detected baseline log2 offset: %.3f", log2_offset))

all_cns_corrected <- all_cns %>%
  mutate(log2_adj = log2 - log2_offset) %>%
  filter(abs(log2_adj) >= 0.4) %>%  # Apply stringent threshold for focal/recurrent CNVs
  mutate(Type = ifelse(log2_adj > 0, "Amp", "Del"))

filtered_cnv_export <- file.path(output_dir, "File_SI_8_cnv_cytoband_filtered.xlsx")
write.xlsx(as.data.frame(all_cns_corrected), filtered_cnv_export)
message("Saved filtered CNV table to: ", filtered_cnv_export)

# ------------------------------------------------------------------------------
# 4. Load Clinical Metadata & Integrate
# ------------------------------------------------------------------------------
message("Step 4: Loading clinical metadata...")
if (!file.exists(clinical_file)) {
  stop("⚠️ Error: Missing clinical metadata file at: ", clinical_file)
}

clinical_data_sync <- read_xlsx(clinical_file)

# Merge CNV records with clinical outcome labels
plot_ready_data <- all_cns_corrected %>%
  left_join(
    clinical_data_sync %>% select(Patient_ID, Event),
    by = c("SampleID" = "Patient_ID")
  ) %>%
  filter(!is.na(Event)) %>%
  mutate(Event = ifelse(Event == 1, "Relapse/Deceased", "Event-free"))

# ------------------------------------------------------------------------------
# 5. Structure Dataset for Genome-Wide Mirrored Landscape
# ------------------------------------------------------------------------------
message("Step 5: Structuring genome-wide cytoband coordinate indices...")
cohort_sample_size <- n_distinct(plot_ready_data$SampleID)

plot_df_indexed <- plot_ready_data %>%
  group_by(Cytoband, chromosome, Type) %>%
  summarise(n = n_distinct(SampleID), .groups = 'drop') %>%
  mutate(freq = n / cohort_sample_size) %>%
  mutate(chrom_label = factor(chromosome, levels = paste0("chr", c(1:22, "X", "Y")))) %>%
  filter(!is.na(chrom_label)) %>%
  arrange(chrom_label, Cytoband) %>%
  mutate(idx = as.numeric(factor(Cytoband, levels = unique(Cytoband)))) %>%
  mutate(freq_plot = ifelse(Type == "Del", -freq, freq))

# Calculate chromosome boundaries and midpoints for x-axis formatting
chr_boundaries <- plot_df_indexed %>% 
  group_by(chrom_label) %>% 
  summarise(max_idx = max(idx), .groups = "drop") %>% 
  pull(max_idx)

chr_labels <- plot_df_indexed %>% 
  group_by(chrom_label) %>% 
  summarise(mid = mean(idx), .groups = "drop")

# Structure subgroup frequency dataset (Relapse/Deceased vs. Event-free)
facet_plot_data <- plot_ready_data %>%
  group_by(Cytoband, Type, Event) %>%
  summarise(n = n_distinct(SampleID), .groups = 'drop') %>%
  left_join(
    plot_ready_data %>% group_by(Event) %>% summarise(group_size = n_distinct(SampleID), .groups = "drop"),
    by = "Event"
  ) %>%
  mutate(
    Group = factor(Event, levels = c("Relapse/Deceased", "Event-free")),
    freq_val = ifelse(Type == "Amp", n / group_size, -(n / group_size))
  ) %>%
  inner_join(plot_df_indexed %>% select(Cytoband, idx, chrom_label) %>% distinct(), by = "Cytoband")

# Select top candidate cytobands for gene annotation
target_labels <- facet_plot_data %>%
  filter(grepl("11q13.3|14p13.3|17p11.2|12q22.3|19p13.2|19q13.1", Cytoband)) %>%
  group_by(Group, Cytoband) %>%
  slice_max(order_by = abs(freq_val), n = 1) %>%
  ungroup() %>%
  mutate(Gene_Label = case_when(
    grepl("11q13.3", Cytoband) ~ "11q13.3 (CCND1/FADD)",
    grepl("14p13.3", Cytoband) ~ "14p13.3 (RNA18SN5)",
    grepl("17p11.2", Cytoband) ~ "17p11.2 (MAP2K3)",
    grepl("12q22.3", Cytoband) ~ "12q22.3 (BTG1)",
    TRUE ~ Cytoband
  ))

# ------------------------------------------------------------------------------
# 6. Plot & Export Mirrored CNV Landscape
# ------------------------------------------------------------------------------
message("Step 6: Rendering mirrored CNV landscape plot...")

p2c <- ggplot(facet_plot_data, aes(x = idx, y = freq_val, fill = Type)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~ Group, ncol = 1, scales = "fixed") +
  geom_vline(xintercept = chr_boundaries, linetype = "dotted", color = "grey80", size = 0.2) +
  
  geom_text_repel(
    data = target_labels,
    aes(x = idx, y = freq_val, label = Gene_Label),
    size = 3.5, fontface = "italic", box.padding = 1.5,
    segment.color = "grey20", segment.size = 0.2, min.segment.length = 0
  ) +
  
  scale_fill_manual(values = c("Amp" = "#d73027", "Del" = "#4575b4")) +
  scale_x_continuous(label = gsub("chr", "", chr_labels$chrom_label), breaks = chr_labels$mid, expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) percent(abs(x)), limits = c(-0.3, 0.3)) +
  geom_hline(yintercept = 0, color = "black", size = 0.4) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    y = "CNV Frequency (%)", 
    x = "Chromosome", 
    title = "B-ALL CNV Landscape: Relapse/Deceased vs. Event-Free"
  )

output_pdf <- file.path(output_dir, "Figure2C_cnv_landscape_mirrored.pdf")
ggsave(output_pdf, plot = p2c, width = 10, height = 6)

message("✅ Successfully executed CNV workflow and generated plot: ", output_pdf)
