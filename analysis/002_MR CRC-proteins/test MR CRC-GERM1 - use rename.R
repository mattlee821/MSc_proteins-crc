rm(list=ls())
set.seed(821)

# environment ====
library(data.table)
library(dplyr)
remotes::install_github("mattlee821/functions")
library(TwoSampleMR)

#Load exposure data for CRC - Try to use command read_outcome_data but keeps getting errors, so I used rename instead, pending for revision
data_exposure_crc <- fread("data/raw/GWAS/crc_combined_GECCO_allancestries.txt.gz") %>%
  rename(
  SNP = SNP,
  beta.exposure = BETA,
  se.exposure = SE,
  effect_allele.exposure = EA,
  other_allele.exposure = OA,
  chr.exposure = CHR,
  pos.exposure = POS,
  pval.exposure = P,
  exposure = phenotype
) %>%
  mutate(
    id.exposure = "CRC", #there is no id.exposure column which create an error during harmonising data 
  )


# Filter for rs2293582 in the exposure data (No need)
#data_exposure_crc <- data_exposure_crc %>%
#  filter(SNP == "rs2293582")

# Load outcome data for plasma proteins
data_outcome_protein <- fread("data/raw/GWAS/plasma-proteins.txt") %>%
  filter(exposure == "18878_15_GREM1_GREM1") %>%
  rename(
    SNP = SNP,
    beta.outcome = beta.exposure,
    se.outcome = se.exposure,
    effect_allele.outcome = effect_allele.exposure,
    other_allele.outcome = other_allele.exposure,
    pval.outcome = pval.exposure,
    outcome = exposure,
    id.outcome = id.exposure,
  )


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
