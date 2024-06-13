rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(dplyr)
remotes::install_github("mattlee821/functions")
library(TwoSampleMR)

# Load and prepare exposure data for CRC
data_exposure_crc <- TwoSampleMR::read_exposure_data(
  filename = "data/raw/GWAS/crc_combined_GECCO_allancestries.txt.gz",
  sep = "\t",
  snp_col = "SNP", 
  id_col = "phenotype",
  beta_col = "BETA", 
  se_col = "SE", 
  pval_col = "P",
  effect_allele_col = "EA", 
  other_allele_col = "OA", 
  chr_col = "CHR", 
  pos_col = "POS"  )
#> No phenotype name specified, defaulting to 'exposure'.

# Load and prepare outcome data for proteins
data_outcome_protein <- TwoSampleMR::read_outcome_data(
  filename = "data/raw/GWAS/plasma-proteins.txt",
  sep = "\t",
  snps = data_exposure_crc$SNP, # Using the same SNP from CRC data
  id_col = "id.exposure", 
  snp_col = "SNP", 
  beta_col = "beta.exposure", 
  se_col = "se.exposure", 
  pval_col = "pval.exposure",
  effect_allele_col = "effect_allele.exposure", 
  other_allele_col = "other_allele.exposure",
  phenotype_col = "exposure")

#filter for 18878_15_GREM1_GREM1
data_outcome_protein <- data_outcome_protein %>%
  filter(outcome == "18878_15_GREM1_GREM1")

# Harmonise data
data_harmonise_crc_protein <- harmonise_data(
  exposure_dat = data_exposure_crc, 
  outcome_dat = data_outcome_protein, 
  action = 2
)

# Perform MR analysis
data_mr_crc_protein <- mr(data_harmonise_crc_protein)

# Display results
print(data_mr_crc_protein)
