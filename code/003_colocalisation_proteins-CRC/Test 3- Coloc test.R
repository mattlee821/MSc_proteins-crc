rm(list = ls())
set.seed(821)

# Load required libraries
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("MRCIEU/genetics.binaRies", force = T)
# remotes::install_github("explodecomputer/plinkbinr", force = F)
# remotes::install_github("chr1swallace/coloc@main", force = T)
# remotes::install_github("sjmgarnier/viridis", force = F)
# remotes::install_github("chr1swallace/coloc@main",build_vignettes=TRUE,build = FALSE, force= T)
# remotes::install_github("mattlee821/functions",build_vignettes=TRUE, build = FALSE)
# install.packages("devtools")
# devtools::install_github("boxiangliu/locuscomparer")
library(genetics.binaRies)
library(coloc)
library(viridis)
library(data.table)
library(functions) 
library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(cowplot)
library(plinkbinr)

# Load data
data_exposure <- fread("analysis/003_colocalisation_proteins-CRC/GREM1-protien-window-coloc.txt")
data_outcome <- fread("analysis/003_colocalisation_proteins-CRC/GREM1-CRC-coloc.txt")
head(data_exposure)
head(data_outcome)

#Select SNP
SNP <- fread("data/raw/GWAS/plasma-proteins.txt") %>%
  filter(exposure == "18878_15_GREM1_GREM1") %>%
  pull(SNP)

# Harmonize data ====
data_exposure$id.exposure <- data_exposure$exposure # Add id.exposure columns
data_outcome$id.outcome <- data_outcome$outcome # Add id.outcome columns
data_harmonise <- harmonise_data(
  exposure_dat = data_exposure,
  outcome_dat = data_outcome,
  action = 2
)
data_harmonise$remove_duplicates <- paste0(data_harmonise$SNP, "_", data_harmonise$id.exposure) # Remove duplicates in Harmonize data
data_harmonise <- data_harmonise[!duplicated(data_harmonise$remove_duplicates),]
list_harmonise <- split(data_harmonise, data_harmonise$id.exposure) # Split harmonized data by id.exposure

#LD matrix ====
ld <- readRDS("data/raw/ld.rds")
head(ld)
str(ld)

# format LD matrix and harmonised list ====
i = 1
ld <- ld[which(rownames(ld) %in% list_harmonise[[i]]$SNP), which(colnames(ld) %in% list_harmonise[[i]]$SNP)]
list_harmonise[[i]] <- list_harmonise[[i]][which(list_harmonise[[i]]$SNP %in% rownames(ld)),]
ld <- ld[match(list_harmonise[[i]]$SNP,rownames(ld)),]
ld <- ld[,match(list_harmonise[[i]]$SNP, colnames(ld))]
list_harmonise[[i]] <- list_harmonise[[i]][match(rownames(ld), list_harmonise[[i]]$SNP),]

#MAF ====
data_maf <- fread("data/references/1000genomes/phase3/processed/ALL/stats.stats")
head(data_maf)
data_maf <- data_maf %>%
  filter(Predictor %in% list_harmonise[[i]]$SNP)

# Define function to calculate standard deviation of Y (sdY) ====
calculate_sdY <- function(beta) {
  sdY <- sd(beta, na.rm = TRUE)
  return(sdY)
}
sdY_exposure <- calculate_sdY(list_harmonise[[i]]$beta.exposure) #Use harmonised data, do we need it? Need it in each window?
sdY_outcome <- calculate_sdY(list_harmonise[[i]]$beta.outcome) #Use harmonised data

# Prepare data for coloc analysis
data_coloc_exposure <- list(
  beta = list_harmonise[[i]]$beta.exposure,
  varbeta = list_harmonise[[i]]$se.exposure^2,
  MAF = list_harmonise[[i]]$eaf.exposure,
  type = "quant",
  N = min(list_harmonise[[i]]$samplesize.exposure),
  snp = list_harmonise[[i]]$SNP,
  sdY = sdY_exposure, 
  LD = ld, 
  position = list_harmonise[[i]]$pos.exposure,
  pval = list_harmonise[[i]]$pval.exposure
)

# Prepare outcome data for coloc analysis
data_coloc_outcome <- list(
  beta = list_harmonise[[i]]$beta.outcome,
  varbeta = list_harmonise[[i]]$se.outcome^2,
  MAF = data_maf$MAF, 
  type = "cc",
  N = 254791,  #actual sample size from the paper
  snp = list_harmonise[[i]]$SNP,
  sdY = sdY_outcome,
  LD = ld,
  position = list_harmonise[[i]]$pos.exposure,
  pval = list_harmonise[[i]]$pval.outcome
  )

# VARS ====
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

window <- list(
  "1mb",
  "1mb",
  "1mb",
  "1mb",
  "1mb"
)

# data check ====
# coloc::check_dataset(d = data_coloc_exposure_test, suffix = 1, warn.minp=5e-8)
# coloc::check_dataset(d = data_coloc_outcome_test, suffix = 2, warn.minp=5e-8)
coloc::check_dataset(d = data_coloc_exposure, suffix = 1, warn.minp=5e-8)
coloc::check_dataset(d = data_coloc_outcome, suffix = 2, warn.minp=5e-8)
## plot dataset
# cowplot::plot_grid(
#   coloc_plot_dataset(d = data_coloc_exposure, label = "exposure"),
#   coloc_plot_dataset(d = data_coloc_outcome, label = "outcome"),
#   coloc_check_alignment(D = data_coloc_exposure),
#   coloc_check_alignment(D = data_coloc_outcome),
#   ncol = 2)

# finemap check ====
SNP_causal_exposure <- coloc::finemap.abf(dataset = data_coloc_exposure) %>%
  filter(SNP.PP == max(SNP.PP)) %>%
  select(snp, SNP.PP)

SNP_causal_outcome <- coloc::finemap.abf(dataset = data_coloc_outcome) %>%
  filter(SNP.PP == max(SNP.PP)) %>%
  select(snp, SNP.PP)

# coloc.abf ====
coloc_results <- lapply(priors, function(params) {
  coloc::coloc.abf(
    dataset1 = data_coloc_exposure,
    dataset2 = data_coloc_outcome,
    p1 = params$p1,
    p2 = params$p2,
    p12 = params$p12
  )
})

# results table ====
table_coloc <- data.frame()
# Loop through the coloc_results list
for (j in seq_along(coloc_results)) {
  results_coloc <- coloc_results[[j]]
  results <- data.frame(
    exposure = "exposure",
    exposure_sex = "exposure_sex",
    exposure_study = "exposure_study",
    exposure_data = "exposure_data",
    exposure_population = "exposure_population",
    outcome = "outcome",
    outcome_sex = "outcome_sex",
    outcome_study = "outcome_study",
    outcome_data = "outcome_data",
    outcome_population = "outcome_population",
    SNP = SNP,
    window = window[[j]],
    finemap_snp_exposure = SNP_causal_exposure$snp,
    finemap_snp_exposure_PP = SNP_causal_exposure$SNP.PP,
    finemap_snp_outcome = SNP_causal_outcome$snp,
    finemap_snp_outcome_PP = SNP_causal_outcome$SNP.PP,
    nsnps = results_coloc$summary[["nsnps"]],
    h0 = results_coloc$summary[["PP.H0.abf"]],
    h1 = results_coloc$summary[["PP.H1.abf"]],
    h2 = results_coloc$summary[["PP.H2.abf"]],
    h3 = results_coloc$summary[["PP.H3.abf"]],
    h4 = results_coloc$summary[["PP.H4.abf"]],
    prior_p1 = results_coloc$priors[["p1"]],
    prior_p2 = results_coloc$priors[["p2"]],
    prior_p12 = results_coloc$priors[["p12"]],
    priors_label = label_priors[j]
  )
  # Bind the current results to the accumulated table
  table_coloc <- bind_rows(table_coloc, results)
}

# sensitivity  ====
plot <- coloc_sensitivity(
  obj = coloc_results[[1]],
  rule = "H4 > 0.8", npoints = 100, row = 1, suppress_messages = TRUE,
  trait1_title = "exposure", trait2_title = "outcome",
  dataset1 = NULL, dataset2 = NULL,
  data_check_trait1 = data_coloc_exposure, data_check_trait2 = data_coloc_outcome
)
plot


# 
# #DRAFT
# # data check 
# coloc::check_dataset(d = coloc_data_exposure, suffix = 1, warn.minp=5e-8) #no idea why got an error - lengths of inputs don't match
# coloc::check_dataset(d = coloc_data_outcome, suffix = 2, warn.minp=5e-8) #no idea why got an error - lengths of inputs don't match
# ## plot dataset
# cowplot::plot_grid(
#   coloc_plot_dataset(d = coloc_data_exposure, label = "exposure"),
#   coloc_plot_dataset(d = coloc_data_outcome, label = "outcome"),
#   coloc_check_alignment(D = coloc_data_exposure),
#   coloc_check_alignment(D = coloc_data_outcome),
#   ncol = 2)
# 
# # finemap check ====
# SNP_causal_exposure <- coloc::finemap.abf(dataset = coloc_data_exposure) %>%
#   filter(SNP.PP == max(SNP.PP)) %>%
#   select(snp, SNP.PP)
# 
# SNP_causal_outcome <- coloc::finemap.abf(dataset = coloc_data_outcome) %>%
#   filter(SNP.PP == max(SNP.PP)) %>%
#   select(snp, SNP.PP)
# 
# # Define windows for analysis based on your data - - Not yet using harmonised data 
# windows <- list(
#   "250kb" = subset(data_exposure, window == "250kb"),
#   "500kb" = subset(data_exposure, window == "500kb"),
#   "1mb" = subset(data_exposure, window == "1mb"),
#   "2mb" = subset(data_exposure, window == "2mb")
# )
# 
# # Define priors for coloc analysis
# priors <- list(
#   list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
#   list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6),
#   list(p1 = 1e-4, p2 = 1e-5, p12 = 1e-6),
#   list(p1 = 1e-5, p2 = 1e-4, p12 = 1e-6),
#   list(p1 = 1e-5, p2 = 1e-5, p12 = 1e-7)
# )
# 
# 
# # Loop over all windows and run colocalization analysis - Not yet using harmonised data 
# for (window_name in names(windows)) {
#   data_window <- windows[[window_name]]
#   
#   #sdY calculate
#   sdY_exposure <- estimate_sdY(data_exposure$beta.exposure)
#   sdY_outcome <- estimate_sdY(data_outcome$beta.outcome)
#   
#   #Prepare data for coloc analysis
#     coloc_data_exposure <- list(
#       beta = data_window_exposure$beta.exposure,
#       varbeta = data_window_exposure$se.exposure^2,
#       MAF = data_window_exposure$ImpMAF,
#       type = "quant",
#       N = data_window_exposure$samplesize.exposure,
#       snp = SNP,
#       LD = ld,
#       sdY = sdY_exposure,
#       position = data_window_exposure$pos.exposure
#     )
# 
#     # Prepare outcome data for coloc analysis
#     coloc_data_outcome <- list(
#       beta = data_outcome$beta.outcome,
#       varbeta = data_outcome$se.outcome^2,
#       MAF = rep(0.01, nrow(data_outcome)), # Assuming default MAF for outcome data
#       type = "cc",
#       N = 6700,  # Adjust based on actual sample size
#       snp = SNP,
#       LD = ld,
#       sdY = sdY_outcome,
#       position = data_outcome$pos.outcome
#     )
#   
#   # Perform colocalization analysis for each prior
#   coloc_results <- lapply(priors, function(params) {
#     coloc::coloc.abf(
#       dataset1 = coloc_data_exposure,
#       dataset2 = coloc_data_outcome,
#       p1 = params$p1,
#       p2 = params$p2,
#       p12 = params$p12
#     )
#   })
#   
#   # Compile results into a table
#   for (j in seq_along(coloc_results)) {
#     results_coloc <- coloc_results[[j]]
#     results <- data.frame(
#       exposure = "GREM1",  
#       outcome = "CRC",  
#       window = "250kb",  
#       nsnps = results_coloc$summary[["nsnps"]],
#       h0 = results_coloc$summary[["PP.H0.abf"]],
#       h1 = results_coloc$summary[["PP.H1.abf"]],
#       h2 = results_coloc$summary[["PP.H2.abf"]],
#       h3 = results_coloc$summary[["PP.H3.abf"]],
#       h4 = results_coloc$summary[["PP.H4.abf"]],
#       prior_p1 = priors[[j]]$p1,
#       prior_p2 = priors[[j]]$p2,
#       prior_p12 = priors[[j]]$p12,
#       priors_label = label_priors[j],
#     )
#     table_coloc <- rbind(table_coloc, results)
#   }
#   
#   # Save results table
#   write.table(table_coloc, paste0("analysis/results_", window_name, ".txt"),
#               row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#   
#   # Save sensitivity plot
#   pdf(paste0("analysis/sensitivity_plot_", window_name, ".pdf"), height = 10, width = 10)
#   coloc_sensitivity(
#     obj = coloc_results[[1]],
#     rule = "H4 > 0.8", npoints = 100, row = 1, suppress_messages = TRUE,
#     trait1_title = "exposure", trait2_title = "CRC",
#     dataset1 = NULL, dataset2 = NULL,
#     data_check_trait1 = data_exposure, data_check_trait2 = data_outcome
#   )
#   dev.off()
# }
# 
