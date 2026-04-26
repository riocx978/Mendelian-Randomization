# Mendelian Randomization Analysis of Periodontal Disease

> Causal inference study identifying genetic risk factors for chronic and acute periodontitis — published as a Master's thesis at the University of South Florida (2024).

📄 [Read the full thesis](https://digitalcommons.usf.edu/etd/10605)

---

## Overview

This study uses **two-sample Mendelian Randomization (MR)** to evaluate whether five modifiable exposures causally influence periodontal disease risk. By using SNPs as genetic instrumental variables, MR sidesteps the confounding and reverse causation problems that plague observational studies.

**Exposures tested:**
- Smoking status
- Maternal smoking after birth
- Plasminogen levels
- Creatinine levels
- Absence of psychosocial stress

**Outcomes:**
- Chronic periodontitis
- Acute periodontitis

---

## Key Findings

| Exposure | Outcome | Direction | Significance |
|----------|---------|-----------|--------------|
| Creatinine levels | Acute periodontitis | Positive ↑ | Significant |
| Maternal smoking | Periodontitis risk | Positive ↑ | Significant |
| Plasminogen levels | Periodontitis risk | Negative ↓ | Significant |
| Smoking status | Chronic periodontitis | Positive ↑ | Borderline |

**Highlights:**
- Creatinine → acute periodontitis link supports the hypothesis that renal dysfunction amplifies systemic inflammation
- Maternal smoking effect aligns with prenatal nicotine exposure literature
- Plasminogen's negative association diverges from prior pro-inflammatory models — flagged as a priority for follow-up research

---

## Methods

All analyses were conducted in **R** using a two-sample MR framework.

**MR methods:**
- Inverse Variance Weighted (IVW) — primary estimate
- MR-Egger — tests and corrects for directional pleiotropy
- Weighted Median — robust when up to 50% of IVs are invalid

**Instrumental variable selection:**
- Genome-wide significance threshold (P < 5×10⁻⁸)
- LD pruning (r² threshold applied)
- Allele harmonization across exposure and outcome GWAS

**Sensitivity analyses:**
- Leave-one-out analysis
- Cochran's Q heterogeneity test
- Funnel plots and scatter plots for visual diagnostics

---

## Repository Structure

```
Mendelian-Randomization/
└── TwoSampleMR.R       # Full MR pipeline: IV selection, harmonization, analysis, sensitivity checks
```

---

## Requirements

```r
install.packages(c("TwoSampleMR", "MRInstruments", "tidyverse"))
```

---

## Citation

Charles, R. (2024). *Unravelling the Impact of Blood Metabolites and Lifestyle Factors on Periodontal Disease Using Mendelian Randomization.* Master's Thesis, University of South Florida. https://digitalcommons.usf.edu/etd/10605

---

## Contact

**Rhea Charles** · [riocx1997@gmail.com](mailto:riocx1997@gmail.com) · [LinkedIn](https://www.linkedin.com/in/rhea-charles/) · [GitHub](https://github.com/riocx978)
