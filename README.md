# B-ALL project

This repository contains the core computational pipelines, analytical scripts, and code required to reproduce the key findings and main figures for the manuscript:

> **Genomic and Subclonal Complexity with B Cell Maturation Associated with Poor Outcome in Paediatric B ALL**


[![R-version](https://img.shields.io/badge/R-%3E%3D_4.2.0-blue.svg)](https://www.r-project.org/)
[![Python-version](https://img.shields.io/badge/Python-%3E%3D_3.9-green.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


---

## 📂 Repository Structure

The code is organized sequentially into modular directories reflecting the computational and analytical workflows:

```text
.
├── scripts/
│   ├── 00_upstream_processing/
│   │   ├── 01_bwa_alignment.sh
│   │   ├── 02_bam_qc.sh
│   │   ├── 03_gatk_markduplicates.sh
│   │   ├── 04_mutect2_calling.sh
│   │   ├── 05_filter_mutect_calls.sh
│   │   ├── 06_annovar_annotation.sh
│   │   ├── 07_scrna_seurat_rpca_harmony.R
│   │   ├── 08_scrna_singler_annotation.R
│   │   └── 09_vdj_mixcr_pipeline.sh
│   │
│   ├── 01_genomic_landscape_fig2/
│   │   ├── 01_oncoplot_landscape.R
│   │   ├── 02_cnv_cytoband_landscape.R
│   │   ├── 03_comparative_mutational_map.R
│   │   ├── 04_cox_high_risk_alterations.R
│   │   └── 05_chromosomal_aneuploidy.R
│   │
│   ├── 02_subclonal_evolution_fig3/
│   │   ├── 01_prepare_subclone_input.R
│   │   ├── 02_pyclone_deconvolution.py
│   │   ├── 03_subclone_evolution_figures.R
│   │   └── 04_subclone_diversity_survival.R
│   │
│   ├── 03_single_cell_dynamics_fig4_fig5/
│   │   ├── 01_reannotate_B_subclusters.R
│   │   ├── 02_cell_proportion_comparison.R
│   │   ├── 03_hsc_richness_pseudotime.R
│   │   └── 04_gene_pseudotime_dynamics.R
│   │
│   ├── 04_microenvironment/
│   │   ├── 01_cytotoxicity_and_nk_receptors.R
│   │   ├── 02_immune_efficacy_and_lr_pairs.R
│   │   ├── 03_pseudotime_ridges_and_density.R
│   │   └── 04_virtual_drug_screening_lgals9.R
│   │
│   └── 05_prognostic_model/
│       ├── 01_univariate_cox_forest.R
│       ├── 02_genomic_score_subclone_correlation.R
│       ├── 03_discovery_cohort_survival.R
│       ├── 04_prepare_target_validation_data.R
│       └── 05_validation_cohort_and_time_roc.R
│
├── README.md
└── LICENSE
