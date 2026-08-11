# B-ALL project

This repository contains the core computational pipelines, analytical scripts, and code required to reproduce the key findings and main figures for the manuscript:

> **Genomic and Subclonal Complexity with B Cell Maturation Associated with Poor Outcome in Paediatric B ALL**

---

## 📌 Overview

This study integrates ultra-deep targeted panel sequencing (~1,500×), whole-exome sequencing (WES), single-cell RNA sequencing (scRNA-seq), and single-cell V(D)J profiling across 105 paediatric B-ALL patients to dissect genomic instability, subclonal evolutionary topologies, and chemoresistant cell reservoirs.

This repository provides:
1. **Upstream Processing Frameworks**: Scripts for somatic variant calling, copy number variation (CNV) profiling, and single-cell RNA/VDJ data integration.
2. **Main Figure Analytical Code**: Scripts to reproduce main figures (Figures 2–7), including subclonal deconvolution, pseudotime trajectory mapping, module scoring, and prognostic model evaluation.

---

## 🛠 Repository Structure

### `scripts/00_upstream_processing/`
- `01_panel_WES_variant_calling.sh`: Somatic variant calling using BWA-MEM, Picard, and GATK Mutect2 in paired and tumour-only modes.
- `02_cnvkit_pipeline.sh`: Copy number segmentation and reference normalization using CNVkit.
- `03_scRNA_seurat_harmony.R`: Quality control, RPCA integration, Harmony batch correction, and UMAP dimensionality reduction using Seurat.

### `scripts/01_genomic_landscape_fig2/`
- `fig2a_oncoplot_landscape.R`: Visualization of recurrent mutational landscapes using `maftools`.
- `fig2c_outcome_associated_variants.R`: Multi-variable Cox regression and alteration enrichment comparisons between outcome groups.
- `fig2e_chromosomal_aneuploidy.R`: Quantifying chromosomal copy number alterations and aneuploidy levels.

### `scripts/02_subclonal_evolution_fig3/`
- `fig3a_pyclone_vi_clustering.py`: Bayesian variational inference for subclonal deconvolution using PyClone-VI.
- `fig3bc_evolutionary_topology.R`: Analysis of truncal vs. branching mutation preferences and testing the threshold effect ($\ge$4 high-risk alterations).
- `fig3_shannon_diversity.R`: Intra-tumour subclonal diversity calculation via CCF-weighted Shannon index ($H'$).

### `scripts/03_single_cell_dynamics_fig4_fig5/`
- `fig4c_hsc_like_proportions.R`: Quantifying HSC-like blast proportions and transcriptional diversity (gene richness per cell).
- `fig4e_monocle3_pseudotime.R`: Single-cell developmental pseudotime trajectory inference using `Monocle3`.
- `fig4g_drug_resistance_score.R`: Stage-specific drug resistance module scoring using `AddModuleScore`.
- `fig5a_mature_B_subclustering.R`: Fine-grained subclustering of mature B-cell reservoirs (*BANK1*, *BCL2*, *CD44*).

### `scripts/04_microenvironment_fig6/`
- `fig6cd_cytotoxicity_scoring.R`: Cytotoxic effector scoring for host T and NK lymphocytes.
- `fig6e_cellchat_interactions.R`: Cell-cell communication and immune checkpoint interaction profiling using `CellChat`.
- `fig6gh_vdj_repertoire.R`: Immune repertoire diversity and top clonotype frequency analysis using `MiXCR` outputs.

### `scripts/05_prognostic_classifier_fig7/`
- `fig7c_multivariable_cox_model.R`: Construction of the evolution-enhanced prognostic classifier combining NCI risk, MRD, molecular subtypes, and subclonal mutation count.
- `fig7e_time_dependent_roc.R`: Time-dependent ROC (24-month AUC) analysis and independent validation on the TARGET B-ALL Phase 2 cohort.

---

## 💻 System Requirements & Environment

The scripts were developed and tested on Linux (Ubuntu 22.04 / CentOS 7) and R (v4.3.3) / Python (v3.12.2).

### Key Dependencies:
- **R packages**: `Seurat` (v5.0), `Monocle3` (v1.3.1), `CellChat`, `maftools`, `timeROC`, `survival`, `survminer`, `ComplexHeatmap`, `tidyverse`.
- **Python libraries**: `PyClone-VI`, `numpy`, `pandas`, `h5py`, `scikit-learn`.
- **Command-line tools**: `GATK` (v4.x), `BWA`, `CNVkit`, `MiXCR` (v3.0.13).
