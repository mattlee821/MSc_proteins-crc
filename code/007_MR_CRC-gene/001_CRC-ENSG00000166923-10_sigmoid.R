###MR - CRC and Gene Expression (Sigmoid)
rm(list=ls())
set.seed(821)

# Load libraries
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(functions)

#Load and Filter Exposure Data (CRC Data)
data_exposure <- fread("data/raw/GWAS/data_exposure.txt") %>%
  filter(id.exposure == 'fernandez-rozadilla_2022_PMID36539618;crc_combined_GECCO_allancestries;CRC') %>%
  filter(pval.exposure <= 5e-4) %>% #Do we need to filter CRC data based on p-value <= 5e-4? Try - Result doesn't change with/without this p-val (102 obs)
  rename(rsid = SNP,
         pval = pval.exposure)

#clumping
data_exposure <- ld_clump(
  dat= data_exposure,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "data/references/1000genomes/phase3/processed/ALL/ALL"
) %>%
  rename(SNP = rsid,
         pval.exposure = pval)

#load and Prepare Outcome Data (Gene Expression)
data_outcome <- fread("data/raw/GWAS/colon-sigmoid_allpairs_ENSG00000166923.10.txt")
data_outcome_ref <- fread("data/raw/GWAS/reference-chr15.txt")
data_outcome <- data_outcome %>%
  left_join(data_outcome_ref, by = "variant_id") %>%
  rename(
      outcome = gene_id,
      SNP = rs_id_dbSNP151_GRCh38p7,           
      beta.outcome = slope,                   
      se.outcome = slope_se,                  
      effect_allele.outcome = alt,     
      other_allele.outcome = ref,      
      eaf.outcome = maf,  
      pval.outcome = pval_nominal
    ) %>%
    mutate(
      id.outcome = outcome) 

# Filter outcome data
data_outcome <- data_outcome %>%
  filter(SNP %in% data_exposure$SNP)

#Harmonise exposure and outcome data 
data_harmonise <- harmonise_data(
  exposure_dat = data_exposure,
  outcome_dat = data_outcome)
## Shorten names ====
# Shorten names for id.exposure, exposure
shorten_names <- function(name) {
  gsub("fernandez-rozadilla_2022_PMID36539618;crc_combined_GECCO_allancestries;CRC", "CRC", name)
}

# Shorten names
data_harmonise$id.exposure <- shorten_names(data_harmonise$id.exposure)
data_harmonise$exposure <- shorten_names(data_harmonise$exposure)
head(data_harmonise)

#Save harmonise data ====
write.table(data_harmonise, "analysis/007_MR_CRC-gene/001_ENSG00000166923-10_sigmoid_data_harmonise_crc-gene.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#Perform MR Analysis ====
data_mr <- mr(data_harmonise)
print(data_mr)

#Create MR Results Table
table_mr <- data_mr %>%
  select(exposure, outcome, method, b, se, pval) %>%
  rename(Exposure = exposure, Outcome = outcome, Method = method, Estimate = b, Std_Error = se, P_Value = pval)
print(table_mr)
write.table(table_mr, "analysis/007_MR_CRC-gene/001_ENSG00000166923-10_sigmoid_table_mr_crc-gene.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
