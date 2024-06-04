rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(dplyr)
remotes::install_github("mattlee821/functions")
library(functions)
library(TwoSampleMR)

# data_exposure ====
data_exposure <- fread("data/raw/GWAS/plasma-proteins.txt")
data_exposure <- data_exposure %>%
  filter(exposure == "18878_15_GREM1_GREM1")

# data_outcome ====
data_outcome <- TwoSampleMR::read_outcome_data(filename = "data/raw/GWAS/crc_combined_GECCO_allancestries.txt.gz", 
                                               sep = "\t",
                                               snps = data_exposure$SNP, 
                                               phenotype_col = "phenotype", 
                                               id_col = "phenotype", 
                                               snp_col = "SNP", 
                                               beta_col = "BETA", 
                                               se_col = "SE", 
                                               pval_col = "P",
                                               effect_allele_col = "EA", 
                                               other_allele_col = "OA", 
                                               chr_col = "CHR", 
                                               pos_col = "POS",
                                               min_pval = 1e-200,
                                               log_pval = FALSE)

## add missing EAF
data_outcome1 <- functions::missing_EAF(df = data_outcome, 
                                        reference = "/data/GWAS_data/work/references/1000genomes/phase3/ALL/stats.stats",
                                        column_EAF = "eaf.outcome")

# data_harmonise ====
data_harmonise <- harmonise_data(exposure_dat = data_exposure, outcome_dat = data_outcome1, action = 2)

# data_mr ====
data_mr <- mr(data_harmonise)

functions::missing_EAF