rm(list = ls())
set.seed(821)

# Load required libraries
library(data.table)
library(dplyr)
library(coloc)
library(cowplot)

# Define paths to PLINK and reference data
plink_path <- "/Users/faye/Desktop/001_projects/MSc_proteins-crc/code/plink1.9/PLINK"
bfile_path <- "/Users/faye/Desktop/001_projects/MSc_proteins-crc/data/references/1000genomes/phase3/processed/ALL/ALL"

# Load data
data_exposure <- fread("analysis/003_colocalisation_proteins-CRC/GREM1-protien-window-coloc.txt")
data_outcome <- fread("analysis/003_colocalisation_proteins-CRC/GREM1-CRC-coloc.txt")
head(data_exposure)
head(data_outcome)

# Define LD estimation function using PLINK
ld_matrix_local <- function(snp_ids, bfile_path, plink_path) {
  ld_matrix <- plinkbinr::plink_ld(
    bfile = bfile_path,
    snps = snp_ids,
    plink = plink_path
  )
  return(ld_matrix)
}

# Define function to calculate standard deviation of Y (sdY)
calculate_sdY <- function(y) {
  # Replace with your method to calculate sdY
  return(sd(y))
}

SNP <- fread("data/raw/GWAS/plasma-proteins.txt") %>%
  filter(exposure == "18878_15_GREM1_GREM1") %>%
  pull(SNP)


#Test - Not seperate window yet
# data check ====
coloc::check_dataset(d = data_exposure, suffix = 1, warn.minp=5e-8)
coloc::check_dataset(d = data_outcome, suffix = 2, warn.minp=5e-8)
## plot dataset
cowplot::plot_grid(
  coloc_plot_dataset(d = data_exposure, label = "exposure"),
  coloc_plot_dataset(d = data_outcome, label = "outcome"),
  coloc_check_alignment(D = data_exposure),
  coloc_check_alignment(D = data_outcome),
  ncol = 2)

# finemap check ====
SNP_causal_exposure <- coloc::finemap.abf(dataset = data_exposure) %>%
  filter(SNP.PP == max(SNP.PP)) %>%
  select(snp, SNP.PP)

SNP_causal_outcome <- coloc::finemap.abf(dataset = data_outcome) %>%
  filter(SNP.PP == max(SNP.PP)) %>%
  select(snp, SNP.PP)


# Define windows for analysis based on your data
windows <- list(
  "250kb" = subset(data_exposure, window == "250kb"),
  "500kb" = subset(data_exposure, window == "500kb"),
  "1mb" = subset(data_exposure, window == "1mb"),
  "2mb" = subset(data_exposure, window == "2mb")
)

# Define priors for coloc analysis
priors <- list(
  list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
  list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6),
  list(p1 = 1e-4, p2 = 1e-5, p12 = 1e-6),
  list(p1 = 1e-5, p2 = 1e-4, p12 = 1e-6),
  list(p1 = 1e-5, p2 = 1e-5, p12 = 1e-7)
)

# Loop over all windows and perform coloc analysis
for (window_name in names(windows)) {
  data_window_exposure <- windows[[window_name]]
  
  # Estimate LD matrix
  ld_matrix <- ld_matrix_local(
    snp_ids = data_window_exposure$SNP,
    bfile_path = bfile_path,
    plink_path = plink_path
  )
  
  # Calculate sdY for exposure data
  sdY_exposure <- calculate_sdY(data_window_exposure$beta.exposure)
  
  # Prepare data for coloc analysis
  coloc_data_exposure <- list(
    beta = data_window_exposure$beta.exposure,
    varbeta = data_window_exposure$se.exposure^2,
    MAF = data_window_exposure$ImpMAF,
    type = "quant",
    N = data_window_exposure$samplesize.exposure,
    snp = rownames(ld_matrix),
    LD = ld_matrix,
    position = data_window_exposure$pos.exposure
  )
  
  # Prepare outcome data for coloc analysis
  coloc_data_outcome <- list(
    beta = data_outcome$beta.outcome,
    varbeta = data_outcome$se.outcome^2,
    MAF = rep(0.01, nrow(data_outcome)), # Assuming default MAF for outcome data
    type = "cc",
    N = 6700,  # Adjust based on actual sample size
    snp = rownames(ld_matrix),
    LD = ld_matrix,
    position = data_outcome$pos.outcome
  )
  
  # Perform coloc analysis for each set of priors
  coloc_results <- lapply(priors, function(params) {
    coloc::coloc.abf(
      dataset1 = coloc_data_exposure,
      dataset2 = coloc_data_outcome,
      p1 = params$p1,
      p2 = params$p2,
      p12 = params$p12
    )
  })
  
  # Output results (replace with your specific requirements)
  for (j in seq_along(coloc_results)) {
    results_coloc <- coloc_results[[j]]
    print(paste("Window:", window_name, "- Prior:", j))
    print(results_coloc)
  }
  
  # Perform data check
  data_check(data_window_exposure)
}
