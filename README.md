```markdown
# Genomic and Subclonal Complexity with B Cell Maturation Associated with Poor Outcome in Paediatric B-ALL

This repository contains upstream computational pipelines, statistical models, figure-generation scripts, and demonstration workflows for reproducing the findings of the manuscript:

> **Genomic and Subclonal Complexity with B Cell Maturation Associated with Poor Outcome in Paediatric B-ALL**

---

## 1. System Requirements

### Hardware Requirements
* **Standard Downstream Analysis & Reproduction**:
  * Operating System: macOS (tested on macOS 14 Sonoma), Linux (tested on Ubuntu 20.04/22.04 LTS), or Windows 10/11 (via WSL2).
  * Processor: Multi-core x86-64 or Apple Silicon processor (4 cores minimum; 8 cores recommended).
  * Memory (RAM): Minimum 16 GB for standard statistical modeling and demo execution; 32 GB recommended for integrated single-cell Seurat object manipulation.
  * Storage: ~2 GB free disk space for cloned repository, demo datasets, and environment dependencies.
* **High-Throughput Upstream Processing (FASTQ to BAM/VCF)**:
  * Recommended: High-Performance Computing (HPC) Linux cluster with Slurm/PBS schedulers.
  * Allocation: ≥16 CPU cores and ≥64 GB RAM per node.

### Software Dependencies & Versions
* **Core Languages**:
  * R (tested on versions >= 4.3.0)
  * Python (tested on versions >= 3.9.0)
  * Bash / GNU Shell (>= 4.4)
* **Key R Packages**:
  * Single-cell analysis: `Seurat` (v4.4.0 / v5.0.0), `SingleR` (v2.4.0), `celldex` (v1.12.0), `monocle3` (v1.3.1)
  * Genomics & Variant analysis: `maftools` (v2.18.0)
  * Survival & Risk modeling: `survival` (v3.5-8), `survminer` (v0.4.9), `timeROC` (v0.4)
  * Data manipulation & Plotting: `tidyverse` (v2.0.0), `pheatmap` (v1.0.12), `openxlsx` (v4.2.5.2), `cowplot` (v1.1.3)
* **Python Environment & Tools**:
  * `PyClone-VI` (v0.1.1) or `PyClone` (v0.13.1)
  * `pandas` (v1.5.3 / v2.0.0), `numpy` (v1.24.3), `scipy` (v1.10.1)
* **Upstream Bioinformatics Tools**:
  * `BWA` (v0.7.17), `SAMtools` (v1.17), `GATK` (v4.4.0.0), `ANNOVAR`, `CNVkit` (v0.9.10), `MiXCR` (v3.0.13 / v4.0.0)

---

## 2. Installation Guide

Setting up the analytical environment via Conda typically takes **10–15 minutes** on a standard desktop computer.

```bash
# 1. Clone repository
git clone [https://github.com/DanYuOmics/BALL_Cohort.git](https://github.com/DanYuOmics/BALL_Cohort.git)
cd BALL_Cohort

# 2. Create and activate a conda environment for R and dependencies
conda create -n ball_omics r-base=4.3.3 python=3.10 -y
conda activate ball_omics

# 3. Install required R packages
R -e 'install.packages(c("tidyverse", "survival", "survminer", "timeROC", "pheatmap", "openxlsx", "cowplot", "BiocManager"), repos="[https://cloud.r-project.org](https://cloud.r-project.org)")'
R -e 'BiocManager::install(c("Seurat", "SingleR", "celldex", "maftools", "fgsea"))'

# 4. Install Python dependencies
pip install pyclone-vi pandas numpy scipy

```

---

## 3. Demo Dataset & Expected Output

To facilitate independent review and rapid code reproduction without downloading multi-terabyte raw sequencing files, a de-identified, lightweight demo dataset is provided in `data/demo/`.

### Demo Files Provided:

* `data/demo/demo_clinical_survival.tsv`: Downsampled clinical features and risk scores (n = 30).
* `data/demo/demo_maf_subset.maf`: Minimal somatic mutation table for testing oncoplot rendering.
* `data/demo/demo_sc_subset.rds`: Downsampled Seurat object containing pre-annotated B-cell subpopulations.

### Running the Demo:

Execute the prognostic model validation and Kaplan-Meier curve generation on the demo cohort:

```bash
# Execute demo survival pipeline
Rscript scripts/05_prognostic_model/03_discovery_cohort_survival.R \
  --input data/demo/demo_clinical_survival.tsv \
  --output_dir results/demo_output

```

### Expected Output & Verification:

* **Output files**:
* `results/demo_output/KM_survival_curve.pdf`
* `results/demo_output/TimeROC_24month.pdf`


* **Expected Result**: Time-dependent ROC curve showing dynamic AUC metrics and Kaplan-Meier stratified survival curves separating risk categories with log-rank statistical outputs.
* **Expected Runtime**: **< 2 minutes** on an ordinary desktop computer (e.g., Apple M-series or Intel Core i7, 16 GB RAM).

---

## 4. Repository Structure

```
.
├── data/
│   └── demo/                               # Lightweight demo datasets for evaluation
├── scripts/
│   ├── 00_upstream_processing/             # Pipeline scripts for raw sequencing data
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
│   ├── 01_genomic_landscape_fig2/          # WES / targeted panel mutation analyses
│   │   ├── 01_oncoplot_landscape.R
│   │   ├── 02_cnv_cytoband_landscape.R
│   │   ├── 03_comparative_mutational_map.R
│   │   ├── 04_cox_high_risk_alterations.R
│   │   └── 05_chromosomal_aneuploidy.R
│   │
│   ├── 02_subclonal_evolution_fig3/        # Subclonal architecture and PyClone workflows
│   │   ├── 01_prepare_subclone_input.R
│   │   ├── 02_pyclone_deconvolution.py
│   │   ├── 03_subclone_evolution_figures.R
│   │   └── 04_subclone_diversity_survival.R
│   │
│   ├── 03_single_cell_dynamics_fig4_fig5/  # scRNA-seq clustering and differentiation trajectories
│   │   ├── 01_reannotate_B_subclusters.R
│   │   ├── 02_cell_proportion_comparison.R
│   │   ├── 03_hsc_richness_pseudotime.R
│   │   └── 04_gene_pseudotime_dynamics.R
│   │
│   ├── 04_microenvironment/                # Immune microenvironment, NK/T cell receptor analysis
│   │   ├── 01_cytotoxicity_and_nk_receptors.R
│   │   ├── 02_immune_efficacy_and_lr_pairs.R
│   │   ├── 03_pseudotime_ridges_and_density.R
│   │   └── 04_virtual_drug_screening_lgals9.R
│   │
│   └── 05_prognostic_model/                # Multi-omics prognostic score and TARGET validation
│       ├── 01_univariate_cox_forest.R
│       ├── 02_genomic_score_subclone_correlation.R
│       ├── 03_discovery_cohort_survival.R
│       ├── 04_prepare_target_validation_data.R
│       └── 05_validation_cohort_and_time_roc.R
│
├── results/                                # Output directory for generated plots and tables
├── LICENSE
└── README.md

```

---

## 5. Instructions for Reproduction on Full Datasets

1. **Controlled Raw Data Access**:
* Raw sequencing reads (targeted capture, WES, scRNA-seq, and V(D)J amplicon profiling) are archived under controlled access in the Genome Sequence Archive for Human (GSA-Human) at the China National Center for Bioinformation (CNCB-NGDC) under BioProject **PRJCA074593** (Submission ID: **subHRA033221**).
* Authorized researchers can apply for data access under ethical and legal compliance via the GSA-Human Data Access Committee (DAC).


2. **Upstream Processing**:
* Upstream Bash scripts in `scripts/00_upstream_processing/` illustrate parameters, alignment against the human reference genome (`GRCh38 / hg38`), MuTect2 variant calling, and SingleR references.


3. **Reproducing Paper Figures**:
* Once processed matrices are placed into the respective working directories, run scripts in folders `01_` through `05_` sequentially to generate the analytical plots, Cox regression summaries, and subclonal evolution trees.



---

## 6. Computational Costs and Parallelization

* **Upstream High-Throughput Processing**:
* Whole-exome and targeted capture alignment, deduplication, and somatic variant calling required approximately 10–14 CPU hours per sample with multi-threading enabled (`-t 8`).
* Total computational allocation for the entire cohort of 105 patients was ~2,500 core-hours on an academic Slurm cluster.


* **Downstream Statistical Modeling & Figure Reproduction**:
* Downstream R and Python analyses (directories `01_` to `05_`) require minimal computational overhead and can be fully executed on a personal desktop computer within 1–2 hours.



---

## License

This project is licensed under the MIT License - see the [LICENSE](https://www.google.com/search?q=LICENSE&utm_source=gemini) file for details.

```

```
