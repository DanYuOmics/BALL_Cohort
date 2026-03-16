好的，我帮你整理一个最终版的 **README**，只保留 **Genomic Data 表格**，并更新了 WES、Panel、scRNA、VDJ、Relapse、Died 的情况。

---

# B-ALL Multi-Omics Project

## Overview
This repository documents data availability, analysis strategies, and codes used for figure generation in the study of pediatric B-ALL (acute lymphoblastic leukemia).  
The project integrates **genomic data (Panel/WES)**, **single-cell RNA sequencing (scRNA)**, and **VDJ repertoire analysis** to provide a comprehensive view of disease biology, relapse mechanisms, and clinical translation.

---

## Part 1: Genomic Data

### Tier Overview
| Tier | Definition        | WES | Panel | Count | Relapse | Died | Notes |
|------|------------------|-----|-------|-------|---------|------|-------|
| 1    | Panel + WES      | 29  | 29    | 29    | 4       | 0    | Validation Set. Use Panel for SNVs, WES for CNVs. |
| 2    | WES Only         | 4   | 0     | 4     | 4       | 0    | Rescue Set. Focus on relapse cases. |
| 3    | Panel Only       | 0   | 71    | 71    | 2       | 2    | Clinical Standard. Use Panel for fusions and SNVs. |
| 4    | No DNA           | 0   | 0     | 9     | 1       | 0    | Excluded from genomic figures. |
| **Total** |              | 33  | 100   | 113   | 11      | 2    | Matches dataset |

- **Relapse patients**: 11 total (Tier 1 = 4, Tier 2 = 4, Tier 3 = 2, Tier 4 = 1).  
- **Deaths**: 2 patients (both in Tier 3).  

### Subtyping Hierarchy Strategy
- **Tier 1 (Panel + WES)**:  
  - Panel → fusions (ETV6-RUNX1, TCF3-PBX1, BCR-ABL1)  
  - WES → ploidy changes (Hyperdiploidy, Hypodiploidy), IKZF1 deletions  
  - Confidence: High  

- **Tier 2 (WES Only)**:  
  - WES → SNVs, CNVs, SVs for all drivers  
  - Focus on relapse cases  
  - Confidence: Medium-High  

- **Tier 3 (Panel Only)**:  
  - Panel → fusions and SNVs (e.g., PAX5)  
  - Limitation: may miss Hyperdiploidy and Ph-like drivers  
  - Confidence: Medium  

- **Tier 4 (No DNA)**:  
  - Do not subtype; label as "Unknown"  
  - Confidence: NA  

---

## Part 2: Single-Cell RNA (scRNA) Data

- **52 scRNA samples** across Tiers 1–4.  
- Includes relapse and non-relapse cases.  
- Contributions: immune microenvironment characterization, validation subset for predictive modeling, integration with genomic drivers.  

---

## Part 3: VDJ Repertoire Data

- **~20 patients** with VDJ sequencing.  
- Overlaps with scRNA samples (e.g., DDN24018220, DDN24020340, DDN24020515, DDN24021572, DDN24022580, DDN24023693, DDN24026188, DDN24026291, DDN24026340, DDN24028415, DDN25001852, DDN25003835, DDN25004233).  
- Contributions: clonal correlation with genomic drivers, immune repertoire diversity analysis.  

---

## Figures Overview
- **Figure 1–2**: Genomic landscape & clonal evolution.  
- **Figure 3–4**: Single-cell atlas & immune TME.  
- **Figure 5**: Pharmacogenomics.  
- **Figure 6**: Multi-omics integration (Genomic + scRNA + VDJ).  
- **Figure 7**: Predictive modeling.  

