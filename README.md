# Mendelian Randomization Analysis of Periodontal Disease

## Overview

This repository contains the code and analysis for a Mendelian Randomization (MR) study investigating the causal effects of multiple exposures on periodontal disease.

The study evaluates the following exposures:

- Smoking status  
- Maternal smoking after birth  
- Plasminogen levels  
- Creatinine levels  
- Absence of psychosocial stress  

Genetic variants (SNPs) were used as instrumental variables (IVs) to assess causal relationships with:

- Chronic periodontitis  
- Acute periodontitis  

---

## Data Sources

### Exposure Data
- Smoking status  
- Maternal smoking after birth  
- Plasminogen levels  
- Creatinine levels  
- Absence of psychosocial stress  

### Outcome Data
- Chronic periodontal disease  
- Acute periodontal disease  

### Instrumental Variables
- SNPs selected based on:
  - Genome-wide significance threshold (P-value)  
  - Linkage disequilibrium threshold (r²)  

---

## Methods

This analysis was conducted in **R** using two-sample Mendelian Randomization approaches.

MR methods implemented:

- Inverse Variance Weighted (IVW)  
- MR-Egger  
- Weighted Median  

Sensitivity analyses included:

- Leave-one-out analysis  
- Heterogeneity testing  
- Graphical diagnostics  

---

## Workflow

1. Selection of instrumental variables based on genome-wide significant SNPs  
2. Harmonization of exposure and outcome datasets to align effect alleles  
3. MR analysis using multiple complementary methods  
4. Sensitivity analyses to assess robustness  

---

## Code Structure

- `TwoSampleMR.R` — Main MR analysis pipeline  

---

## Requirements

Install the required R packages:

```r
install.packages(c("TwoSampleMR", "MRInstruments", "tidyverse"))
```

## Key Findings
This study provides insights into the causal relationships between selected risk factors and periodontitis:
### Creatinine Levels and Acute Periodontitis  
A significant positive association was observed, supporting evidence that renal dysfunction may exacerbate inflammatory conditions.
### Maternal Smoking and Periodontitis Risk  
A significant positive association was identified, aligning with prior research on prenatal nicotine exposure and adverse health outcomes.
### Plasminogen Levels and Periodontitis Risk  
A negative association was detected, differing from earlier studies suggesting a pro-inflammatory role of plasminogen. This highlights the need for further investigation.
### Smoking and Chronic Periodontitis  
A borderline significant positive association was observed, reinforcing smoking as a risk factor.
Sensitivity analyses and graphical evaluations supported the robustness of these findings.

## Contact  
For any questions, please contact  
Rhea Charles at riocx1997@gmail.com or  
open an issue on this repository.
