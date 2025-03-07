# Mendelian-Randomization
#Overview
This repository contains the code and analysis for a Mendelian Randomization (MR) study investigating the association of multiple exposures with periodontal disease. The goal of this study is to assess [causal relationships between exposures: Smoking Status, Maternal Smoking After Birth effects, Plasminogen, Creatinine Levels and Absence of Psychosocial Stress with Periodontal Disease using genetic variants as instrumental variables (IVs).
\\
Data Sources
Exposure Data: Smoking Status, Maternal Smoking After Birth effects, Plasminogen, Creatinine Levels and Absence of Psychosocial Stress
Outcome Data: Chronic and Acute Periodontal Disease
Instrumental Variables: SNPs selected based on P value threshold and r2 

Methods
This analysis was performed using Two-Sample MR, MR-Egger, IVW, Weighted Median] implemented in R. 
\\
The workflow includes:

Selection of instrumental variables based on genome-wide significant SNPs.
Harmonization of exposure and outcome datasets to align effect alleles.
\\
MR Analysis using methods such as:
Inverse Variance Weighted (IVW)
MR-Egger
Weighted Median
Sensitivity Analyses (leave-one-out, heterogeneity tests)
\\
Code Structure
TwoSampleMR.R - Main MR pipeline
\\

To reproduce the analysis, install the following packages:
``r
install.packages(c("TwoSampleMR", "MRInstruments", "tidyverse"))
For Python (if applicable):
\\
Key findings include:

[Summarize main findings, e.g., significant causal relationship between X and Y]
[Interpretation and implications of the results]
\\
To run the MR analysis, execute the following script:
``r
Rscript scripts/mr_analysis.R
\\
#Contact
For any questions, please contact Rhea Charles at riocx1997@gmail.com or open an issue on this repository.
