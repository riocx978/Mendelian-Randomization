# =============================================================================
# Mendelian Randomization Analysis: Periodontal Disease
# Author: Rhea Charles
# Institution: University of South Florida
# Published thesis: https://digitalcommons.usf.edu/etd/10605
#
# Description:
#   Two-sample MR analysis evaluating causal effects of five exposures on
#   chronic and acute periodontitis. Includes IV selection, harmonization,
#   MR estimation (IVW, MR-Egger, Weighted Median), sensitivity analyses,
#   and visualization (forest plots, heatmap, funnel plots).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

options(scipen = 999)  # Suppress scientific notation
options(ieugwasr_api = "gwas-api.mrcieu.ac.uk/")

library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(gt)
library(meta)
library(ComplexHeatmap)
library(reshape2)
library(grid)


# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# Update these paths to point to your local data directory

data_path <- "data/"   # <-- set your data directory here

# Outcome GWAS summary statistics
perio_acute_ukbb <- vroom(paste0(data_path, "AcutePerioUKBB.txt"))
perio_acute      <- vroom(paste0(data_path, "AcutePerio.txt"))
perio_chronic    <- vroom(paste0(data_path, "ChronicPerio.txt"))

# Exposure GWAS summary statistics
exposure_creatinine   <- vroom(paste0(data_path, "exposure3.txt"))  # Creatinine levels
exposure_mat_smoking  <- vroom(paste0(data_path, "exposure4.txt"))  # Maternal smoking after birth
exposure_plasminogen  <- vroom(paste0(data_path, "exposure5.txt"))  # Plasminogen levels
exposure_smoking      <- vroom(paste0(data_path, "exposure6.txt"))  # Smoking status
exposure_stress       <- vroom(paste0(data_path, "exposure7.txt"))  # Absence of psychosocial stress


# -----------------------------------------------------------------------------
# 2. Format outcome data
# -----------------------------------------------------------------------------

format_outcome <- function(dat, outcome_label) {
  dat %>%
    na.omit() %>%
    format_data(
      type               = "outcome",
      snp_col            = "SNP",
      beta_col           = "beta_EUR",
      se_col             = "se_EUR",
      effect_allele_col  = "ref",
      other_allele_col   = "alt",
      eaf_col            = "af_cases_EUR",
      pval_col           = "neglog10_pval_EUR"
    ) %>%
    mutate(outcome = outcome_label)
}

perio_acute_fmt <- format_outcome(perio_acute_ukbb, "Acute Perio")
perio_chronic_fmt <- format_outcome(perio_chronic, "Chronic Perio")

vroom_write(perio_acute_fmt, "output/AcutePerioUKBB_Adj.txt")


# -----------------------------------------------------------------------------
# 3. Format and clump exposure data
# -----------------------------------------------------------------------------
# Clumping selects independent, genome-wide significant SNPs as IVs
# Thresholds: p < 5e-8, r2 < 0.001, window = 10,000 kb

format_and_clump <- function(dat, exposure_label) {
  dat %>%
    na.omit() %>%
    format_data(
      type                    = "exposure",
      snp_col                 = "variant_id",
      beta_col                = "beta",
      se_col                  = "standard_error",
      effect_allele_col       = "effect_allele",
      other_allele_col        = "other_allele",
      eaf_col                 = "effect_allele_frequency",
      pval_col                = "p_value"
    ) %>%
    mutate(exposure = exposure_label) %>%
    filter(pval.exposure <= 5e-8) %>%
    clump_data(
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p1 = 5e-8,
      clump_p2 = 1,
      pop      = "EUR"
    )
}

tophits_creatinine  <- format_and_clump(exposure_creatinine,  "Creatinine")
tophits_mat_smoking <- format_and_clump(exposure_mat_smoking, "Maternal Smoking")
tophits_plasminogen <- format_and_clump(exposure_plasminogen, "Plasminogen")
tophits_smoking     <- format_and_clump(exposure_smoking,     "Smoking")
tophits_stress      <- format_and_clump(exposure_stress,      "Absence of Psychosocial Stress")


# -----------------------------------------------------------------------------
# 4. Harmonize exposure and outcome data
# -----------------------------------------------------------------------------
# Aligns effect alleles between exposure and outcome GWAS

harmonized_creatinine  <- harmonise_data(tophits_creatinine,  perio_acute_fmt)
harmonized_mat_smoking <- harmonise_data(tophits_mat_smoking, perio_acute_fmt)
harmonized_plasminogen <- harmonise_data(tophits_plasminogen, perio_acute_fmt)
harmonized_smoking     <- harmonise_data(tophits_smoking,     perio_chronic_fmt)
harmonized_stress      <- harmonise_data(tophits_stress,      perio_acute_fmt)


# -----------------------------------------------------------------------------
# 5. MR analysis
# -----------------------------------------------------------------------------
# Three complementary methods:
#   IVW           — primary estimate; assumes all SNPs are valid IVs
#   MR-Egger      — detects and corrects for directional pleiotropy
#   Weighted Median — robust when up to 50% of IVs are invalid

mr_methods <- c(
  "mr_egger_regression",
  "mr_ivw",
  "mr_ivw_radial",
  "mr_ivw_mre",
  "mr_weighted_median"
)

run_mr <- function(harmonized_dat) {
  mr(harmonized_dat, method_list = mr_methods) %>%
    generate_odds_ratios()
}

mr_creatinine  <- run_mr(harmonized_creatinine)
mr_mat_smoking <- run_mr(harmonized_mat_smoking)
mr_plasminogen <- run_mr(harmonized_plasminogen)
mr_smoking     <- run_mr(harmonized_smoking)
mr_stress      <- run_mr(harmonized_stress)

# Combine all results
mr_combined <- bind_rows(
  mr_creatinine,
  mr_mat_smoking,
  mr_plasminogen,
  mr_smoking,
  mr_stress
)

# Label grouping: show exposure name only on first row of each group
mr_combined <- mr_combined %>%
  group_by(exposure) %>%
  mutate(group = ifelse(row_number() == 1, exposure, "")) %>%
  ungroup()


# -----------------------------------------------------------------------------
# 6. Sensitivity analyses
# -----------------------------------------------------------------------------

run_sensitivity <- function(harmonized_dat, label) {
  cat("\n===", label, "===\n")

  cat("\nHeterogeneity test:\n")
  print(mr_heterogeneity(harmonized_dat))

  cat("\nPleiotropy test (MR-Egger intercept):\n")
  print(mr_pleiotropy_test(harmonized_dat))

  loo <- mr_leaveoneout(harmonized_dat)
  mr_leaveoneout_plot(loo)

  return(loo)
}

loo_creatinine  <- run_sensitivity(harmonized_creatinine,  "Creatinine")
loo_mat_smoking <- run_sensitivity(harmonized_mat_smoking, "Maternal Smoking")
loo_plasminogen <- run_sensitivity(harmonized_plasminogen, "Plasminogen")
loo_smoking     <- run_sensitivity(harmonized_smoking,     "Smoking")
loo_stress      <- run_sensitivity(harmonized_stress,      "Absence of Psychosocial Stress")


# -----------------------------------------------------------------------------
# 7. Outlier detection and removal
# -----------------------------------------------------------------------------
# Flags SNPs whose removal shifts the MR estimate by > 2 SD

remove_outliers <- function(harmonized_dat, exposure_dat) {
  loo     <- mr_leaveoneout(harmonized_dat)
  mean_b  <- mean(loo$b)
  sd_b    <- sd(loo$b)

  outlier_snps <- loo %>%
    filter(abs(b - mean_b) > 2 * sd_b) %>%
    pull(SNP)

  if (length(outlier_snps) > 0) {
    cat("Outlier SNPs removed:", paste(outlier_snps, collapse = ", "), "\n")
  } else {
    cat("No outliers detected.\n")
  }

  exposure_dat %>% filter(!SNP %in% outlier_snps)
}

tophits_creatinine_clean <- remove_outliers(harmonized_creatinine, tophits_creatinine)


# -----------------------------------------------------------------------------
# 8. Visualization
# -----------------------------------------------------------------------------

# --- 8a. Forest plot (meta-analysis across all exposures and methods) ---

meta_all <- metagen(
  TE     = mr_combined$b,
  seTE   = mr_combined$se,
  studlab = mr_combined$method,
  data   = mr_combined,
  sm     = "OR"
)

png("output/forest_plot.png", width = 15, height = 7.5, units = "in", res = 300)
forest(
  meta_all,
  studlab              = TRUE,
  comb.fixed           = FALSE,
  comb.random          = TRUE,
  print.tau2           = FALSE,
  xlab                 = "Log Odds Ratio for Periodontitis",
  leftcols             = c("group", "method", "nsnp", "b", "se", "pval"),
  leftlabs             = c("Exposure", "Method", "SNPs", "Effect Size", "SE", "P-value"),
  rightcols            = c("or", "or_lci95", "or_uci95"),
  rightlabs            = c("OR", "Lower CI", "Upper CI"),
  fontsize             = 12,
  colgap.forest.left   = "0.8cm",
  colgap.forest.right  = "0.8cm",
  digits               = 2,
  backtransf           = TRUE,
  col.square           = "red",
  col.diamond.random   = "blue",
  col.diamond.lines    = "blue",
  col.predict          = "purple"
)
dev.off()


# --- 8b. Heatmap of odds ratios across exposures and methods ---

heatmap_data <- mr_combined %>%
  mutate(
    significance = case_when(
      pval < 0.05             ~ "*",
      pval >= 0.05 & pval < 0.06 ~ "+",
      TRUE                    ~ ""
    )
  )

melted <- melt(
  heatmap_data,
  id.vars      = c("exposure", "method", "outcome", "pval", "significance"),
  measure.vars = "or"
)

median_or <- median(melted$value)

heatmap_plot <- ggplot(melted, aes(x = exposure, y = method, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = significance), color = "black", size = 3) +
  scale_fill_gradient2(
    low      = "blue",
    high     = "red",
    mid      = "lightblue",
    midpoint = median_or,
    name     = "Odds Ratio"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.caption = element_text(hjust = 0.5, size = 10)
  ) +
  labs(
    title   = "Odds Ratios: Exposures vs. Periodontitis",
    x       = "Exposure",
    y       = "MR Method",
    caption = "Significance: (*) p < 0.05, (+) p 0.05–0.06"
  )

ggsave("output/heatmap.png", plot = heatmap_plot, width = 12, height = 8, bg = "white")


# --- 8c. Per-exposure diagnostic plots ---

plot_diagnostics <- function(harmonized_dat) {
  sp   <- mr(harmonized_dat)
  mr_scatter_plot(sp, harmonized_dat)
  mr_forest_plot(mr_singlesnp(harmonized_dat))
  mr_funnel_plot(mr_singlesnp(harmonized_dat))
}

plot_diagnostics(harmonized_creatinine)
plot_diagnostics(harmonized_mat_smoking)
plot_diagnostics(harmonized_plasminogen)
plot_diagnostics(harmonized_smoking)
plot_diagnostics(harmonized_stress)


# --- 8d. Significant SNP forest plot (filters out null SNPs) ---

plot_significant_snps <- function(harmonized_dat, title_label) {
  snp_data <- mr_singlesnp(harmonized_dat) %>%
    mutate(
      lower_ci = b - 1.96 * se,
      upper_ci = b + 1.96 * se
    ) %>%
    filter(!(lower_ci < 0 & upper_ci > 0))  # Remove SNPs crossing the null

  ggplot(snp_data, aes(x = b, y = SNP)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(
      title = paste("Significant SNPs:", title_label),
      x     = "Effect Size (log OR)",
      y     = "SNP"
    ) +
    theme(
      plot.title  = element_text(hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
}

plot_significant_snps(harmonized_creatinine,  "Creatinine")
plot_significant_snps(harmonized_mat_smoking, "Maternal Smoking")
plot_significant_snps(harmonized_plasminogen, "Plasminogen")
plot_significant_snps(harmonized_smoking,     "Smoking")
plot_significant_snps(harmonized_stress,      "Absence of Psychosocial Stress")
