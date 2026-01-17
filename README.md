## BALL_leukemia
## This README describes data availability and its contribution to figures.

## Genomic Data Summary
## Tier Overview
| Tier | Definition     | Count | Relapse Cases | Strategy for Paper |
|------|----------------|-------|---------------|--------------------|
| 1    | Panel + WES    | 29    | 4             | Validation Set. Use Panel for SNVs, WES for CNVs. |
| 2    | WES Only       | 4     | 4             | Rescue Set. |
| 3    | Panel Only     | 63    | 2             | Standard Set. |
| 4    | No DNA         | 9     | 1             | Exclude from Figs 1 & 7. |
| **Total** |         | **105** | **11**        |                    |

---

## Subtyping Hierarchy Strategy
- **Tier 1 (Panel + WES): Gold Standard**
  - Panel: fusions (ETV6-RUNX1, TCF3-PBX1, BCR-ABL1)
  - WES: ploidy changes (High Hyperdiploidy, Hypodiploidy), IKZF1 deletions
  - Confidence: High

- **Tier 2 (WES Only): Rescue Set (Critical for Relapse)**
  - WES for all drivers
  - Challenge: fusions harder to detect; check SVs in ETV6, KMT2A, TCF3
  - If fusion missed, check CNVs (High Hyperdiploidy)
  - Confidence: Medium-High

- **Tier 3 (Panel Only): Clinical Standard**
  - Panel for fusions and SNVs (PAX5)
  - Limitation: may miss High Hyperdiploidy and Ph-like drivers
  - Action: classify by available drivers; if none, label "B-other"
  - Confidence: Medium

- **Tier 4 (No DNA)**
  - Do not subtype
  - Label as "Unknown" in Table 1

---

## Analysis Assignment Matrix (Per Figure)

### Figure 1: Genomic Landscape
- **Oncoplot (Panel B):** Tiers 1 + 2 + 3 (add Data Source track)
- **CNV Landscape (Panel C):** Tiers 1 + 2 (WES), Tier 3 (targeted CNVs only)
- **Volcano Plots (Panel D/E):** Tiers 1 + 2 + 3 (maximize statistical power)

### Figure 2: Clonal Evolution (Deep Genomics)
- **Clonal Architecture (Panel A):** Tiers 1 + 2 (WES required)
- **Evolutionary Trees (Panel B):** Select 2 patients from Tier 1 or 2
- **Mutational Signatures (Panel C):** Tiers 1 + 2 (COSMIC signatures need WES)

### Figures 3 & 4: Single-Cell Atlas & Immune TME
- Use all scRNA samples (N=52) across Tiers 1, 3, and 4
- Tier 4 contributes via scRNA (T-cells, B-blasts)

### Figure 5: Pharmacogenomics
- **Germline Heatmap (Panel B):** Tiers 1 + 2 + 3
- If Tier 3 lacks NUDT15/TPMT, use Tiers 1 + 2 only

### Figure 6: Multi-Omics Integration
- **Driver vs. Phenotype (Panel A/B):** Intersection of (Tiers 1+2+3) AND (scRNA samples)
- **Genomic-Clonal Correlation (Panel C):** Intersection of (Tiers 1+2+3) AND (VDJ samples)

### Figure 7: Predictive Model
- **Nomogram/LASSO:** Tiers 1 + 2 + 3 (N=96)
- Exclude Tier 4 (missing genomic subtype)
- **Validation (Panel E):** scRNA subset (N=52)

