# Setting options to suppress scientific notation
options(scipen = 999)

# Setting the working directory to where the data files are located
setwd("C:/path")

# Loading necessary libraries for data manipulation and analysis
library(readr)      # For reading data files
library(vroom)      # For fast data reading
library(tidyr)      # For data tidying
library(tibble)     # For modern data frames
library(dplyr)      # For data manipulation
library(TwoSampleMR) # For Mendelian Randomization (MR) analysis
library(ggplot2)    # For data visualization
library(gt)         # For creating tables
library(meta)       # For meta-analysis
library(ComplexHeatmap) # For creating complex heatmaps
library(reshape2)   # For reshaping data
library(grid)       # For grid graphics

# Define file paths for the data
data_path <- 'C:/path'
data_path1 <- 'C:/path'

# Load outcome and exposure data
perio <- vroom(paste0(data_path, "AcutePerioUKBB.txt")) # Acute Perio outcome data
perio2 <- vroom(paste0(data_path, "AcutePerio.txt")) # Acute Perio outcome data (alternative)
perio3 <- vroom(paste0(data_path, "ChronicPerio.txt")) # Chronic Perio outcome data
exposure3 <- vroom(paste0(data_path, "exposure3.txt")) # Exposure data for creatinine levels
exposure4 <- vroom(paste0(data_path, "exposure4.txt")) # Exposure data for maternal smoking around birth
exposure5 <- vroom(paste0(data_path, "exposure5.txt")) # Exposure data for plasminogen levels
exposure6 <- vroom(paste0(data_path, "exposure6.txt")) # Exposure data for smoking
exposure7 <- vroom(paste0(data_path, "exposure7.txt")) # Exposure data for absence of psychosocial stress

# Prepare outcome data for MR analysis
perio <- na.omit(perio) # Remove rows with missing values
perio <- format_data(perio, type = "outcome",
                     snp_col = "SNP",
                     beta_col = "beta_EUR",
                     se_col = "se_EUR",
                     effect_allele_col = "ref",
                     other_allele_col = "alt",
                     eaf_col = "af_cases_EUR",
                     pval_col = "neglog10_pval_EUR") %>%
  mutate(outcome = 'Acute Perio')
vroom_write(perio, 'AcutePerioUKBB_Adj.txt') # Save formatted data

# Prepare exposure data for MR analysis
exposure <- na.omit(exposure) # Remove rows with missing values
exposure <- format_data(exposure, type = "exposure",
                        snp_col = "variant_id",
                        beta_col = "beta",
                        se_col = "standard_error",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",
                        eaf_col = "effect_allele_frequency",
                        pval_col = "p_value") %>%
  mutate(exposure = 'vitDdeficiency')
vroom_write(exposure, 'vitDdeficiency_ebi_Adj.txt') # Save formatted data

# Set API endpoint for GWAS data
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

# Clump exposure data to identify top hits
tophits <- exposure %>% 
  filter(pval.exposure <= 5e-05) %>% # Filter for significant SNPs
  clump_data(.,
             clump_kb = 10000,
             clump_r2 = 0.001,
             clump_p1 = 5e-05,
             clump_p2 = 1,
             pop = "EUR")
vroom_write(tophits, 'exposure_tophits.txt') # Save clumped data

# Harmonize outcome and exposure data
Harmonized <- harmonise_data(tophits, perio) #done multiple times for each exposure and outcome

# Perform MR analysis with different methods
mr_res3 <- mr(exposure3, method_list = c("mr_egger_regression", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_weighted_median"))
mr_res4 <- mr(filtered_df, method_list = c("mr_egger_regression", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_weighted_median"))
mr_res5 <- mr(exposure5, method_list = c("mr_egger_regression", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_weighted_median"))
mr_res6 <- mr(exposure6, method_list = c("mr_egger_regression", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_weighted_median"))
mr_res7 <- mr(exposure7, method_list = c("mr_egger_regression", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_weighted_median"))

# Generate odds ratios for all MR results
odds_ratios_list3 <- generate_odds_ratios(mr_res3)
odds_ratios_list4 <- generate_odds_ratios(mr_res4)
odds_ratios_list5 <- generate_odds_ratios(mr_res5)
odds_ratios_list6 <- generate_odds_ratios(mr_res6)
odds_ratios_list7 <- generate_odds_ratios(mr_res7)

# Combine results into one data frame
combined_mr_res <- rbind(odds_ratios_list3, odds_ratios_list4, odds_ratios_list5, odds_ratios_list6, odds_ratios_list7)

# Group similar methods and create a column to show exposure name only once
combined_mr_res$group <- as.character(combined_mr_res$exposure)
combined_mr_res$group <- ave(combined_mr_res$group, combined_mr_res$group, FUN = function(x) {
  x[1] <- x[1]  # Keep the first occurrence
  x[-1] <- ""   # Blank out subsequent occurrences
  return(x)
})

# Create a meta-analysis object for all instruments
meta_analysis_all <- metagen(
  TE = combined_mr_res$b,
  seTE = combined_mr_res$se,
  studlab = paste(combined_mr_res$method),  # Use only the method as the label
  data = combined_mr_res,
  sm = "OR" # Specify the effect measure as Odds Ratio
)

# Save the forest plot as a PNG file
png("forest_plot.png", width = 15, height = 7.5, units = "in", res = 300)

# Create and customize the forest plot
forest(meta_analysis_all, 
       studlab = TRUE, 
       comb.fixed = FALSE, 
       comb.random = TRUE, 
       print.tau2 = FALSE, 
       xlab = "Log Odds Ratio for Periodontitis",
       leftcols = c("group", "method", "nsnp", "b", "se", "pval"),
       leftlabs = c("Exposure", "Method", "SNPs", "Effect Size", "SE", "P-value"),
       rightcols = c("or", "or_lci95", "or_uci95"),
       rightlabs = c("OR", "Lower CI", "Upper CI"),
       fontsize = 12, 
       colgap.forest.left = "0.8cm", 
       colgap.forest.right = "0.8cm",
       digits = 2, 
       backtransf = TRUE, 
       col.square = "red", 
       col.square.random = "red", 
       col.study = "black", 
       col.diamond.random = "blue", 
       col.diamond.lines = "blue", 
       col.predict = "purple")

dev.off() # Close the device

# Create a matrix for the heatmap with Odds Ratios
data <- combined_mr_res %>%
  mutate(exposure_outcome = paste(exposure, outcome, sep = "_"))

melted_data <- melt(data, id.vars = c("exposure", "method", "outcome", "pval"), measure.vars = "or")

# Add significance annotations based on p-values
melted_data$significance <- case_when(
  pval < 0.05 ~ "*",
  pval >= 0.05 & pval < 0.06 ~ "+",
  TRUE ~ ""
)

# Calculate the median value for the midpoint of the color scale
median_value <- median(melted_data$value)

# Create and customize the heatmap
heatmap <- ggplot(melted_data, aes(x = exposure, y = method, fill = value)) + 
  geom_tile(color = "white") +
  geom_text(aes(label = paste(significance)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "lightblue", 
                       midpoint = median_value, limit = c(min(melted_data$value), max(melted_data$value)), 
                       space = "Lab", name = "Odds Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.caption = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Heatmap of Exposures with Outcome Periodontitis",
       x = "Exposure", y = "Method", caption = "Significance: (*) p < 0.05, (+) 0.05 <= 0.06")

# Add a side key note
heatmap <- heatmap + annotate("text", x = Inf, y = length(unique(melted_data$method)) + 1, 
                              label = "+ p-values slightly above 0.05 (between 0.05 and 0.06)", 
                              hjust = 1.1, vjust = 1.1, size = 3, color = "black")

# Save the heatmap as a PNG file
ggsave("heatmap.png", plot = heatmap, width = 12, height = 8, bg = "white")

# Perform MR analysis with the harmonized data
mr(Harmonized, method_list = c("mr_egger_regression", "mr_ivw", "mr_ivw_radial", "mr_ivw_mre", "mr_weighted_median"))

# Conduct heterogeneity and pleiotropy tests
b.1test = mr_heterogeneity(Harmonized)
mr_heterogeneity(Harmonized)

b.2test = mr_pleiotropy_test(Harmonized)
mr_pleiotropy_test(Harmonized)

# Perform Leave-One-Out sensitivity analysis
b.3test = mr_leaveoneout(Harmonized)
mr_leaveoneout_plot(b.3test)

# Generate MR scatter plot
sp1 = mr(Harmonized)
mr_scatter_plot(sp1, Harmonized)

# Generate MR forest plot
mr_forest_plot(mr_singlesnp(Harmonized))

# Generate MR funnel plot
mr_funnel_plot(mr_singlesnp(Harmonized))

#check presence of outliers
mr_rucker_cooksdistance(exposure, parameters = default_parameters())
# Remove outliers based on Leave-One-Out plot & Test thresholds for filtering
loo_results <- mr_leaveoneout(exposure)
mr_leaveoneout_plot(loo_results)

# Calculate mean and standard deviation of MR estimates
mean_effect <- mean(loo_results$b)
sd_effect <- sd(loo_results$b)

# Threshold for identifying significant SNPs
threshold_sd <- 2 * sd_effect 

# Identify and filter out SNPs with significant changes
significant_snps <- loo_results[abs(loo_results$b - mean_effect) > threshold_sd, ]
filtered_data <- loo_results[!loo_results$SNP %in% significant_snps$SNP, ]
mr_leaveoneout_plot(filtered_data)

# Print significant SNPs
print("Significant SNPs (Standard Deviation Approach):")
print(significant_snps$SNP)

filtered_df <- exposure %>% filter(!SNP %in% significant_snps$SNP)

# Extract data from mr_singlesnp for forest plot
singlesnp_data <- mr_singlesnp(Harmonized)
singlesnp_data$lower_ci <- singlesnp_data$b - 1.96 * singlesnp_data$se
singlesnp_data$upper_ci <- singlesnp_data$b + 1.96 * singlesnp_data$se

# Define no-effect value
no_effect <- 0

# Filter out non-significant points
significant_data <- singlesnp_data[!(singlesnp_data$lower_ci < no_effect & singlesnp_data$upper_ci > no_effect),]

# Plot the forest plot for significant SNPs
ggplot(significant_data, aes(x = b, y = SNP)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2) +
  geom_vline(xintercept = no_effect, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Relevant SNPs: Exposure name", x = "Effect Size", y = "SNP") +
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.margin = margin(20, 20, 20, 20) 
  )