# B-ALL Multi-Omics Project

## Overview
This README documents data availability, analysis strategies, and codes used for figure generation in the study of B-ALL (acute lymphoblastic leukemia).  
The project integrates **genomic data**, **single-cell RNA sequencing (scRNA)**, and **VDJ repertoire analysis** to provide a comprehensive view of disease biology.

---

## Part 1: Genomic Data

### Tier Overview
| Tier | Definition     | Count | Relapse Cases | Strategy for Paper |
|------|----------------|-------|---------------|--------------------|
| 1    | Panel + WES    | 29    | 4             | Validation Set. Use Panel for SNVs, WES for CNVs. |
| 2    | WES Only       | 4     | 4             | Rescue Set. |
| 3    | Panel Only     | 63    | 2             | Standard Set. |
| 4    | No DNA         | 9     | 1             | Exclude from Figs 1 & 7. |
| **Total** |         | **105** | **11**        |                    |

### Subtyping Hierarchy Strategy
- **Tier 1 (Panel + WES): Gold Standard**
  - Panel: fusions (ETV6-RUNX1, TCF3-PBX1, BCR-ABL1)
  - WES: ploidy changes (High Hyperdiploidy, Hypodiploidy), IKZF1 deletions
  - Confidence: High
- **Tier 2 (WES Only): Rescue Set**
  - WES for all drivers; focus on relapse cases
  - Confidence: Medium-High
- **Tier 3 (Panel Only): Clinical Standard**
  - Panel for fusions and SNVs (PAX5)
  - Limitations: may miss High Hyperdiploidy and Ph-like drivers
  - Confidence: Medium
- **Tier 4 (No DNA)**
  - Do not subtype; label as "Unknown"

### Analysis Assignment Matrix (Per Figure)
- **Figure 1:** Genomic Landscape (Tiers 1–3)
- **Figure 2:** Clonal Evolution (Tiers 1–2 only)
- **Figure 3–4:** scRNA Atlas & Immune TME (all scRNA samples including Tier 4)
- **Figure 5:** Pharmacogenomics (Tiers 1–3)
- **Figure 6:** Multi-Omics Integration (intersection of DNA + scRNA/VDJ)
- **Figure 7:** Predictive Model (Tiers 1–3, validation with scRNA subset)

---

## Part 2: Single-Cell RNA (scRNA) Data

### Data Availability
- 52 scRNA samples across Tiers 1, 2, 3, and 4.

### Contributions
- Immune microenvironment characterization (Figures 3 & 4).
- Validation subset for predictive modeling (Figure 7).
- Integration with genomic drivers (Figure 6).

### Analysis Highlights
- Cell type annotation (T-cells, B-blasts, myeloid subsets).
- Differential expression and pathway enrichment.
- Immune TME profiling.

---

## Part 3: VDJ Repertoire Data

### Data Availability
- Subset of patients with VDJ sequencing.

### Contributions
- Clonal correlation with genomic drivers (Figure 6C).
- Immune repertoire diversity analysis.

### Analysis Highlights
- Clonal expansion patterns.
- Correlation with relapse and treatment response.

---

## Codes and Pipelines
All scripts for preprocessing, analysis, and figure generation are introduced in this README.  
- **Genomic analysis:** Panel/WES harmonization, CNV calling, driver classification.  
- **scRNA analysis:** Cell clustering, annotation, immune TME profiling.  
- **VDJ analysis:** Clonal assignment, diversity metrics, integration with scRNA.  

---

## Figures Overview
- **Figure 1–2:** Genomic landscape & clonal evolution.  
- **Figure 3–4:** Single-cell atlas & immune TME.  
- **Figure 5:** Pharmacogenomics.  
- **Figure 6:** Multi-omics integration (Genomic + scRNA + VDJ).  
- **Figure 7:** Predictive modeling.  

