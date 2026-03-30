# R Pipeline for PTM Proteomics of Age- and Variant-Associated Modification Dynamics in Drosophila

Analytical pipeline for an unpublished PTM proteomics study examining how post-translational modifications change across two experimental dimensions in *Drosophila melanogaster*: aging (Group 1 vs. Group 3) and genotype (Variant vs. Control). This extends the analytical framework used in a published study from the same lab to a new PTM dataset with a more complex two-factor design.

> **Note:** This repository contains scripts only. Raw mass spectrometry data and processed quantification files are not included, as the underlying manuscript is in preparation. Variant identity, specific modification types, tissue identity, and experimental conditions have been anonymized pending publication.

---

## Related Publication

Curley M., Rai M., Chuang C-L., et al. "Transgenic sensors reveal compartment-specific effects of aggregation-prone proteins on subcellular proteostasis during aging." *Cell Reports Methods*, October 2024.

The current dataset represents additional PTM proteomics analyses associated with this research program. The analytical approach — modification site identification, fraction-level quantification, limma-based statistical testing, and GO enrichment — is adapted from prior work and extended here to accommodate the two-factor experimental design.

---

## Study Overview

This project quantifies post-translational modification changes across two experimental dimensions in *Drosophila melanogaster*:

**Dimension 1 — Age:** Group 1 (baseline) vs. Group 3 (aged), analyzing how modifications accumulate or are lost with aging in a variant genetic background.

**Dimension 2 — Genotype:** Variant vs. Control at the same age group, identifying modifications specifically associated with variant expression.

**Six comparisons analyzed:**

| Comparison label | Biological meaning |
|---|---|
| `Variant_Group3_vs_Group1_Ratio` | Age effect on fraction ratio (Variant background) |
| `Variant_vs_Control_Group3_Ratio` | Genotype effect on fraction ratio |
| `Variant_Group3_vs_Group1_SOL` | Age effect on Fraction A (soluble) |
| `Variant_vs_Control_Group3_SOL` | Genotype effect on Fraction A |
| `Variant_Group3_vs_Group1_INSOL` | Age effect on Fraction B (insoluble) |
| `Variant_vs_Control_Group3_INSOL` | Genotype effect on Fraction B |

**Statistical approach:** limma with Benjamini–Hochberg FDR correction. Significance threshold: adjusted p < 0.05. No fold-change threshold applied in enrichment analyses.

---

## Repository Structure

```
R-Pipeline-Drosophila-PTM-Variant-Aging/
├── 01_data-processing/
│   └── 01_data_processing_visualization.R
└── 02_GO-enrichment/
    └── 02_GO_enrichment.R
```

---

## Scripts

### `01_data-processing/`

**`01_data_processing_visualization.R`** — End-to-end pipeline covering six steps:

**Step 1 — Data loading and contaminant removal.** Loads raw limma results from two Excel workbooks (fraction ratio comparisons and fraction-level comparisons). Identifies and removes contaminant proteins (keratins, albumin, trypsin, immunoglobulins, and common bovine proteins) while retaining legitimate *Drosophila* proteins that may match generic contaminant patterns. Clean datasets are cached to `clean_data/` to avoid reprocessing on subsequent runs.

**Step 2 — Long-format organization.** Merges all six comparison sheets into a single long-format data frame (`organized_PTM_data.csv`), with columns for gene name, protein ID, modification type, modification site, fold change, p-value, and comparison identity. This file serves as the primary input for Script 02.

**Step 3 — Volcano plots.** Generates labeled and unlabeled volcano plots for each of the six comparisons. Points are colored by modification type; significance defined by adj.p < 0.05 and |log2FC| > 1. Top 20 sites labeled by combined significance and fold-change rank.

**Step 4 — Global modification pattern plots.** Summarizes the percentage of sites changed per modification type per comparison, displayed as a faceted barplot (colored by mean log2FC) and a tile heatmap. These global views show which modification types are most responsive to age vs. genotype differences.

**Step 5 — Overview UpSet plot.** UpSet plot of all significant sites across all six comparisons, showing how sites are shared or unique across comparisons. Exported with a summary CSV identifying proteins in each intersection.

**Step 6 — Per-modification-type UpSet plots.** Individual UpSet plots for each modification type, using `ComplexHeatmap::UpSet` for publication-quality output. Each plot shows site-level sharing across comparisons, with modification-specific colors and gene name + modification site labels. Summary CSVs list the specific proteins in each intersection pattern.

**Step 7 — Publication-quality scatter plots.** Three scatter plots with PTM-type color coding and R² annotations: (A) Fraction A FC vs. Fraction B FC for the age comparison, (B) Fraction A FC vs. Fraction B FC for the genotype comparison, (C) Fraction A FC vs. fraction ratio FC for the age comparison. Non-significant sites shown as grey background; significant sites colored by modification type and sized by -log10(adj.p) of the ratio comparison. Top 15 sites labeled with gene name and modification site.

### `02_GO-enrichment/`

**`02_GO_enrichment.R`** — GO overrepresentation analysis via `clusterProfiler::enrichGO`. Runs BP, MF, and CC ontologies for each of the six comparisons. Includes a gene name cleaning step that handles common *Drosophila* annotation quirks (e.g., "Dmel\" prefixes, isoform suffixes like -RA/-RB) before ID mapping. Exports unmapped significant proteins as a QC file, enrichment result tables (CSV), and dotplots (PDF) for all significant comparisons.

---

## Dependencies

```r
# CRAN
install.packages(c(
  "tidyverse", "readxl", "ggrepel", "UpSetR", "RColorBrewer"
))

# Bioconductor
BiocManager::install(c(
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Dm.eg.db",
  "enrichplot"
))
```

**R version:** Developed under R 4.3+.

---

## Data Requirements

Scripts expect input files in the working directory. The anonymized file names below reflect what is used in the scripts.

| File | Used by | Description |
|---|---|---|
| `dataA_ratio_comparisons.xlsx` | Script 01 | Fraction ratio limma results (2 sheets) |
| `dataA_levels_comparisons.xlsx` | Script 01 | Fraction-level limma results (4 sheets) |
| `organized_PTM_data.csv` | Script 02 | Output of Script 01 Step 2 |

**Expected sheet names** in `dataA_ratio_comparisons.xlsx`: `Variant_Group3_vs_Group1`, `Variant_vs_Control_Group3`

**Expected sheet names** in `dataA_levels_comparisons.xlsx`: `INSOL_Variant_Group3_vs_Group1`, `INSOL_Variant_vs_Control_Group3`, `SOL_Variant_Group3_vs_Group1`, `SOL_Variant_vs_Control_Group3`

Update these sheet names in Script 01 Step 1 if your files use different conventions.

---

## Key Analytical Decisions

- **Two-factor design:** Age and genotype comparisons are run in parallel rather than as an interaction model, allowing independent assessment of each effect while using the same analytical framework.
- **Contaminant retention logic:** Common *Drosophila* housekeeping proteins (Act5C, Gapdh1, Gapdh2, etc.) are explicitly retained even if their names overlap with contaminant patterns, preventing false removal of biologically relevant proteins.
- **Clean data caching:** Script 01 checks for existing cleaned CSVs before reprocessing raw Excel files, saving substantial runtime on iterative analysis.
- **Gene name cleaning:** Script 02 handles three common *Drosophila* annotation artifacts before Entrez ID mapping: "Dmel\" prefixes, isoform suffixes (-RA, -RB), and "Isoform" text in gene names.
- **No fold-change threshold** in GO enrichment: significance is determined by adjusted p-value alone (adj.p < 0.05), consistent with the other repositories in this portfolio.
- **Background for ORA:** All proteins detected in the dataset, not the full proteome.

---

## Citation

Manuscript in preparation. Repository will be updated with full citation upon publication.

---

## Author

**Chia-Lung Chuang**
Postdoctoral Research Associate
Department of Developmental Neurobiology
St. Jude Children's Research Hospital
Memphis, TN

PI: Dr. Fabio Demontis
