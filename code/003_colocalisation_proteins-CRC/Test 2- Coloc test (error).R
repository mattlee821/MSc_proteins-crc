rm(list=ls())
set.seed(821)

# environment ====
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("MRCIEU/genetics.binaRies", force = T)
# remotes::install_github("explodecomputer/plinkbinr", force = F)
# remotes::install_github("chr1swallace/coloc@main", force = T)
# remotes::install_github("sjmgarnier/viridis", force = F)
# remotes::install_github("chr1swallace/coloc@main",build_vignettes=TRUE, force= T)
# remotes::install_github("mattlee821/functions",build_vignettes=TRUE)
# remotes::install_github("explodecomputer/plinkbinr", force= T)
library(genetics.binaRies)
library(plinkbinr)
library(coloc)
library(viridis)
library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(cowplot)

#Load data 
data_exposure <- fread ("analysis/003_colocalisation_proteins-CRC/GREM1-protien-window-coloc.txt")
data_outcome <- fread ("analysis/003_colocalisation_proteins-CRC/GREM1-CRC-coloc.txt")
head(data_exposure)
head(data_outcome)

# Paths to PLINK and reference data
plink_path <- "/Users/faye/Desktop/001_projects/MSc_proteins-crc/code/plink1.9/PLINK"
bfile_path <- "/Users/faye/Desktop/001_projects/MSc_proteins-crc/data/references/1000genomes/phase3/processed/ALL/ALL"

## Harmonize data
data_exposure$id.exposure <- data_exposure$exposure # Add id.exposure columns
data_outcome$id.outcome <- data_outcome$outcome # Add id.outcome columns
data_harmonise <- harmonise_data(
  exposure_dat = data_exposure, 
  outcome_dat = data_outcome,
  action = 2
)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure) # Remove duplicates in Harmonize data
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
data_harmonise_coloc <- split(data_harmonise, data_harmonise$id.exposure) # Split harmonized data by id.exposure
# write.table(data_harmonise, "analysis/003_colocalisation_proteins-CRC/GREM1-data_harmonise_coloc.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t") # Save Harmonize data file
head(data_harmonise)

# Define windows for analysis based on your data
windows <- list(
  "250kb" = subset(data_harmonise, window == "250kb"),
  "500kb" = subset(data_harmonise, window == "500kb"),
  "1mb" = subset(data_harmonise, window == "1mb"),
  "2mb" = subset(data_harmonise, window == "2mb")
)

# Define LD estimation function
ld_matrix_local <- function(snp_ids, bfile_path, plink_path) {
  ld_matrix <- plinkbinr::plink_ld(
    bfile = bfile_path,
    snps = snp_ids,
    plink = plink_path
  )
  return(ld_matrix)
}

# Initialize an empty dataframe for final results
table_coloc <- data.frame()

# Define priors and labels
priors <- list(
  list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
  list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6),
  list(p1 = 1e-4, p2 = 1e-5, p12 = 1e-6),
  list(p1 = 1e-5, p2 = 1e-4, p12 = 1e-6),
  list(p1 = 1e-5, p2 = 1e-5, p12 = 1e-7)
)

label_priors <- c(
  "p1=1e-4;p2=1e-4;p12=1e-5",
  "p1=1e-4;p2=1e-4;p12=1e-6",
  "p1=1e-4;p2=1e-5;p12=1e-6",
  "p1=1e-5;p2=1e-4;p12=1e-6",
  "p1=1e-5;p2=1e-5;p12=1e-7"
)

# Loop over all windows and run colocalization analysis
for (window_name in names(windows)) {
  data_window_harmonised <- windows[[window_name]]
  
  # Make LD matrix
  ld_matrix <- ld_matrix_local(
    snp_ids = data_window_harmonised$SNP,
    bfile_path = bfile_path,
    plink_path = plink_path
  )
  
  # Prepare data for colocalization
  coloc_data_exposure <- list(
    beta = data_window_harmonised$beta.exposure,
    varbeta = data_window_harmonised$se.exposure^2,
    MAF = data_window_harmonised$ImpMAF,
    type = "quant",
    N = data_window_harmonised$samplesize.exposure,
    snp = rownames(ld_matrix),
    LD = ld_matrix,
    position = data_window_harmonised$POS
  )
  
  # For outcome, assuming a default MAF value (e.g., 0.01) since outcome MAF is not available
  default_MAF_outcome <- 0.01
  
  coloc_data_outcome <- list(
    beta = data_window_harmonised$beta.outcome,
    varbeta = data_window_harmonised$se.outcome^2,
    MAF = rep(default_MAF_outcome, length(data_window_harmonised$SNP)), 
    type = "cc",
    N = 6700,
    snp = rownames(ld_matrix),
    LD = ld_matrix,
    position = data_window_harmonised$POS
  )
  
  # Perform colocalization analysis for each prior
  coloc_results <- lapply(priors, function(params) {
    coloc::coloc.abf(
      dataset1 = coloc_data_exposure,
      dataset2 = coloc_data_outcome,
      p1 = params$p1,
      p2 = params$p2,
      p12 = params$p12
    )
  })
  
  # Compile results into a table
  for (j in seq_along(coloc_results)) {
    results_coloc <- coloc_results[[j]]
    results <- data.frame(
      exposure = "GREM1",  
      outcome = "CRC",  
      window = window_name,  # Update with actual window name
      nsnps = results_coloc$summary[["nsnps"]],
      h0 = results_coloc$summary[["PP.H0.abf"]],
      h1 = results_coloc$summary[["PP.H1.abf"]],
      h2 = results_coloc$summary[["PP.H2.abf"]],
      h3 = results_coloc$summary[["PP.H3.abf"]],
      h4 = results_coloc$summary[["PP.H4.abf"]],
      prior_p1 = priors[[j]]$p1,
      prior_p2 = priors[[j]]$p2,
      prior_p12 = priors[[j]]$p12,
      priors_label = label_priors[j],
      finemap_snp_exposure = NA,  # To be edited
      finemap_snp_exposure_PP = NA,   # To be edited
      finemap_snp_outcome = NA,  # To be edited
      finemap_snp_outcome_PP = NA   # To be edited
    )
    table_coloc <- rbind(table_coloc, results)
  }
  
  # Save results table
  write.table(table_coloc, paste0("analysis/results_", window_name, ".txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  # Save sensitivity plot
  pdf(paste0("analysis/sensitivity_plot_", window_name, ".pdf"), height = 10, width = 10)
  coloc_sensitivity(
    obj = coloc_results[[1]],
    rule = "H4 > 0.8", npoints = 100, row = 1, suppress_messages = TRUE,
    trait1_title = "exposure", trait2_title = "CRC",
    dataset1 = NULL, dataset2 = NULL,
    data_check_trait1 = data_exposure, data_check_trait2 = data_outcome
  )
  dev.off()
}
