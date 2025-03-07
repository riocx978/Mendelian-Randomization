# Mendelian-Randomization
#Overview
This repository contains the code and analysis for a Mendelian Randomization (MR) study investigating the association of multiple exposures with periodontal disease. The goal of this study is to assess [causal relationships between exposures: Smoking Status, Maternal Smoking After Birth effects, Plasminogen, Creatinine Levels and Absence of Psychosocial Stress with Periodontal Disease using genetic variants as instrumental variables (IVs).

Data Sources
Exposure Data: Smoking Status, Maternal Smoking After Birth effects, Plasminogen, Creatinine Levels and Absence of Psychosocial Stress
Outcome Data: Chronic and Acute Periodontal Disease
Instrumental Variables: SNPs selected based on P value threshold and r2 

Methods
This analysis was performed using Two-Sample MR, MR-Egger, IVW, Weighted Median] implemented in R. 

The workflow includes:

Selection of instrumental variables based on genome-wide significant SNPs.
Harmonization of exposure and outcome datasets to align effect alleles.

MR Analysis using methods such as:
Inverse Variance Weighted (IVW)
MR-Egger
Weighted Median
Sensitivity Analyses (leave-one-out, heterogeneity tests)

Code Structure
TwoSampleMR.R - Main MR pipeline

To reproduce the analysis, install the following packages:
install.packages(c("TwoSampleMR", "MRInstruments", "tidyverse"))

Key findings include:
This study provides nuanced insights into the causal relationships between selected risk factors and periodontitis, revealing both expected and novel associations:
Creatinine levels and acute periodontitis: A significant positive association was observed, supporting prior findings that renal dysfunction may exacerbate inflammatory conditions, including periodontal disease.
Maternal smoking and periodontitis risk: A significant positive association was found, aligning with previous research linking prenatal nicotine exposure to adverse health outcomes in offspring.
Plasminogen levels and periodontitis risk: Contrary to expectations, a negative association was detected, diverging from earlier studies suggesting a positive role of plasminogen in inflammatory processes. This discrepancy underscores the need for further investigation, considering the complex interplay of genetic and environmental factors.
Smoking levels and chronic periodontitis: A borderline significant positive effect was observed, reinforcing smoking as a risk factor but with variations in the strength of association across studies.
Sensitivity analyses and graphical evaluations confirmed the robustness of these findings, refining the understanding of each exposure’s role in periodontitis risk.

#Contact
For any questions, please contact Rhea Charles at riocx1997@gmail.com or open an issue on this repository.
