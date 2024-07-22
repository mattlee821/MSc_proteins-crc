rm(list=ls())
set.seed(821)

# Load libraries
library(data.table)
library(dplyr)
library(TwoSampleMR)
remotes::install_github("mattlee821/functions")
library(gridExtra)
library(ggplot2)
library(patchwork)

# Load and prepare exposure data for CRC
data_exposure_crc <- fread("data/raw/GWAS/data_exposure.txt")
data_exposure_crc <- data_exposure_crc %>%
  filter(id.exposure == 'fernandez-rozadilla_2022_PMID36539618;crc_combined_GECCO_allancestries;CRC')

# Load and prepare outcome data for proteins
data_outcome_protein <- fread("data/raw/GWAS/18878_15_GREM1_GREM1.txt.gz.annotated.gz.exclusions.gz.alleles.gz",
                              header = FALSE, sep = "\t",
                              col.names = c("chr.outcome", "pos.outcome", "SNPID", "SNP", "A1", "A2",
                                            "beta.outcome", "pval.outcome", "min_log10_pval", "se.outcome",
                                            "samplesize.outcome", "ImpMAF", "outcome", "effect_allele.outcome",
                                            "other_allele.outcome", "eaf.outcome"))

# Extract the protein name from the outcome column
extract_protein_name <- function(outcome) {
  sub("^.*/(.*)\\.txt\\.gz\\.unzipped$", "\\1", outcome)
}

# Shorten names for id.exposure, exposure, and outcome
shorten_names <- function(name) {
  gsub("fernandez-rozadilla_2022_PMID36539618;crc_combined_GECCO_allancestries;CRC", "CRC", name)
}

# Filter and prepare outcome data
data_outcome_protein <- data_outcome_protein %>%
  filter(SNP %in% data_exposure_crc$SNP) %>%
  mutate(id.outcome = extract_protein_name(outcome))

# Harmonise data
data_harmonise_crc_protein <- harmonise_data(
  exposure_dat = data_exposure_crc,
  outcome_dat = data_outcome_protein,
  action = 2
)

# Shorten names in data_harmonise_crc_protein
data_harmonise_crc_protein$id.exposure <- shorten_names(data_harmonise_crc_protein$id.exposure)
data_harmonise_crc_protein$exposure <- shorten_names(data_harmonise_crc_protein$exposure)
data_harmonise_crc_protein$outcome <- extract_protein_name(data_harmonise_crc_protein$outcome)
head(data_harmonise_crc_protein)
data_harmonise_crc_protein$outcome <- "GREM1" #shorten
head(data_harmonise_crc_protein)
write.table(data_harmonise_crc_protein, "analysis/002_MR_CRC-proteins/data_harmonise_crc-protein.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Perform MR analysis
data_mr_crc_protein <- mr(data_harmonise_crc_protein)

# Perform MR for each SNP (i.e., Wald ratios)
data_mr_singlesnp <- mr_singlesnp(data_harmonise_crc_protein)

# Test for heterogeneity in MR estimates across individual SNPs
data_mr_heterogeneity <- mr_heterogeneity(data_harmonise_crc_protein, method_list = c("mr_ivw"))

# Check the intercept term in MR Egger regression
data_mr_pleiotropy_test <- mr_pleiotropy_test(data_harmonise_crc_protein)

# Check for outlying SNPs using leave-one-out analysis
data_mr_leaveoneout <- mr_leaveoneout(data_harmonise_crc_protein)

#Table for MR
data_mr_crc_protein
table_mr_crc_protein <- data_mr_crc_protein %>%
  select(exposure, outcome, method, b, se, pval) %>%
  rename(Exposure = exposure, Outcome = outcome, Method = method, Estimate = b, Std_Error = se, P_Value = pval)
print(table_mr_crc_protein)
write.table(table_mr_crc_protein, "analysis/002_MR_CRC-proteins/table_mr_crc-protien_mr-nonexclusion.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#Table for sensitivity 
table_sensitivity <- data.frame(
  Test = c("Heterogeneity Test", "Pleiotropy Test"),
  Statistic = c(data_mr_heterogeneity$Q, data_mr_pleiotropy_test$intercept),
  P_value = c(data_mr_heterogeneity$pval, data_mr_pleiotropy_test$pval)
)
print(table_sensitivity)
write.table(table_sensitivity, "analysis/002_MR_CRC-proteins/table_mr_crc-protien_sensitivity-nonexclusion.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# Generate plots
plot_mr_scatter_plot <- mr_scatter_plot(data_mr_crc_protein, data_harmonise_crc_protein)[[1]]
plot_mr_forest_plot <- mr_forest_plot(data_mr_singlesnp)[[1]]
plot_mr_leaveoneout_plot <- mr_leaveoneout_plot(data_mr_leaveoneout)[[1]]
plot_mr_funnel_plot <- mr_funnel_plot(data_mr_singlesnp)[[1]]
plot_mr_density_plot <- mr_density_plot(data_mr_singlesnp, data_mr_crc_protein)[[1]]

# Mange legends
plot_mr_scatter_plot <- plot_mr_scatter_plot + theme(legend.position = 'bottom')
plot_mr_forest_plot <- plot_mr_forest_plot + theme(legend.position = 'bottom')
plot_mr_leaveoneout_plot <- plot_mr_leaveoneout_plot+ theme(legend.position = 'bottom')
plot_mr_funnel_plot <- plot_mr_funnel_plot+ theme(legend.position = 'bottom')
plot_mr_density_plot <- plot_mr_density_plot+ theme(legend.position = 'bottom')

# Print plots
print(plot_mr_scatter_plot)
print(plot_mr_forest_plot)
print(plot_mr_leaveoneout_plot)
print(plot_mr_funnel_plot)
print(plot_mr_density_plot)

# Combine selected plots
combined_plot <- plot_mr_scatter_plot + plot_mr_forest_plot
print(combined_plot)

### Exclusion of outlier
data_harmonise_exclusion <- data_harmonise_crc_protein %>%
  filter(SNP != "rs2293582")
write.table(data_harmonise_exclusion, "analysis/002_MR_CRC-proteins/data_harmonise_crc-protein_exclusion(rs2293582).txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Perform MR and sensitivity analysis after exclusion
data_mr_exclusion <- mr(data_harmonise_exclusion)
data_mr_singlesnp_exclusion <- mr_singlesnp(data_harmonise_exclusion)
data_mr_heterogeneity_exclusion <- mr_heterogeneity(data_harmonise_exclusion, method_list = c("mr_ivw"))
data_mr_pleiotropy_test_exclusion <- mr_pleiotropy_test(data_harmonise_exclusion)
data_mr_leaveoneout_exclusion <- mr_leaveoneout(data_harmonise_exclusion)

# Table for MR analysis after exclusion
table_mr_exclusion <- data_mr_exclusion %>%
  select(exposure, outcome, method, b, se, pval) %>%
  rename(Exposure = exposure, Outcome = outcome, Method = method, Estimate = b, Std_Error = se, P_Value = pval)
print(table_mr_exclusion)
write.table(table_mr_exclusion, "analysis/002_MR_CRC-proteins/table_mr_crc-protien_mr-exclusion.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#Table for sensitivity after exclusion
table_sensitivity_exclusion <- data.frame(
  Test = c("Heterogeneity Test", "Pleiotropy Test"),
  Statistic = c(data_mr_heterogeneity_exclusion$Q, data_mr_pleiotropy_test_exclusion$intercept),
  P_value = c(data_mr_heterogeneity_exclusion$pval, data_mr_pleiotropy_test_exclusion$pval)
)
print(table_sensitivity_exclusion)
write.table(table_sensitivity_exclusion, "analysis/002_MR_CRC-proteins/table_mr_crc-protien_sensitivity-exclusion.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


# Generate plots after exclusion
plot_mr_scatter_plot_exclusion <- mr_scatter_plot(data_mr_exclusion, data_harmonise_exclusion)[[1]]
plot_mr_forest_plot_exclusion <- mr_forest_plot(data_mr_singlesnp_exclusion)[[1]]
plot_mr_leaveoneout_plot_exclusion <- mr_leaveoneout_plot(data_mr_leaveoneout_exclusion)[[1]]
plot_mr_funnel_plot_exclusion <- mr_funnel_plot(data_mr_singlesnp_exclusion)[[1]]
plot_mr_density_plot_exclusion <- mr_density_plot(data_mr_singlesnp_exclusion, data_mr_exclusion)[[1]]

# Mange legends
plot_mr_scatter_plot_exclusion <- plot_mr_scatter_plot_exclusion + theme(legend.position = 'bottom')
plot_mr_forest_plot_exclusion <- plot_mr_forest_plot_exclusion + theme(legend.position = 'bottom')
plot_mr_leaveoneout_plot_exclusion <- plot_mr_leaveoneout_plot_exclusion+ theme(legend.position = 'bottom')
plot_mr_funnel_plot_exclusion <- plot_mr_funnel_plot_exclusion+ theme(legend.position = 'bottom')
plot_mr_density_plot_exclusion <- plot_mr_density_plot_exclusion+ theme(legend.position = 'bottom')

# Print plots after exclusion
print(plot_mr_scatter_plot_exclusion)
print(plot_mr_forest_plot_exclusion)
print(plot_mr_leaveoneout_plot_exclusion)
print(plot_mr_funnel_plot_exclusion)
print(plot_mr_density_plot_exclusion)

# Combine selected plots
combined_plot_exclusion <- plot_mr_scatter_plot_exclusion + plot_mr_forest_plot_exclusion
print(combined_plot_exclusion)
