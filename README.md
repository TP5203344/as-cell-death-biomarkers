# Cell Death-Related Diagnostic Biomarkers in Human Atherosclerosis

**Identification of Cell Death-related Diagnostic Biomarkers and Immune Infiltration Features in Human Atherosclerosis: A Multi-stage Bioinformatics Study**

This repository contains the complete R/Python analysis code for the manuscript. All analyses are fully reproducible using the public datasets described below.

## Datasets

All datasets are publicly available from the Gene Expression Omnibus (GEO):

| Dataset | Platform | Samples | Role |
|---------|----------|---------|------|
| GSE100927 | Agilent GPL17077 | 69 AS vs 35 control | Discovery (DEG) |
| GSE43292 | Affymetrix GPL6244 | 32 early vs 32 advanced plaques | WGCNA + training |
| GSE28829 | Affymetrix GPL5175 | 13 early vs 16 advanced plaques | External validation |
| GSE115469 | scRNA-seq | 8,439 cells (carotid plaque) | Single-cell localization |

GWAS summary statistics for Mendelian randomization:
- Exposure: eQTLGen consortium cis-eQTLs
- Outcome: UK Biobank CAD GWAS (ebi-a-GCST005194) via OpenGWAS/IEU

## Pipeline Overview

```
Figure 1: DE analysis (limma) → volcano/heatmap/sankey-bubble plots
    ↓
40 cell-death-related candidate genes (DEG ∩ cell death gene sets)
    ↓
Figure 3: LASSO + RF-RFE → 17 core genes
    ↓
Figure 3: ROC / Nomogram / Calibration / DCA
    ↓
Figure 4: CIBERSORT immune infiltration → consensus clustering (K=3) → 17-gene × 22-immune-cell correlation → DCA
    ↓
Figure 4/6: Two-sample MR (eQTLGen → CAD) — AIM2 protective, MPO risk
    ↓
Figure 5: scRNA-seq localization (Seurat) — pyroptosis/NETosis genes in macrophages/monocytes
```

## Directory Structure

```
Figure1/  DE analysis and visualization (volcano, heatmap, sankey-bubble)
Figure2/  WGCNA module identification
Figure3/  Machine learning (LASSO + RF-RFE), validation, nomogram
Figure4/  Immune infiltration, immune subtype, MR analysis, DCA
Figure5/  Single-cell RNA-seq analysis (Seurat)
```

## Requirements

- R ≥ 4.4.0 with packages: limma, WGCNA, glmnet, caret, pROC, rms, Seurat, clusterProfiler, org.Hs.eg.db, CIBERSORT (LM22 signature)
- Python ≥ 3.11 with: pandas, numpy, matplotlib, scipy
- OpenGWAS JWT token for MR API access (register at https://api.opengwas.io/)

## Citation

Tao Y, Rong N, Shi S, Ding Y. Identification of Cell Death-related Diagnostic Biomarkers and Immune Infiltration Features in Human Atherosclerosis: A Multi-stage Bioinformatics Study. *Gene* (under review).

## License

MIT License
