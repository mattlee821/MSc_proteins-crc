#006-1_colocalisation_gene-CRC (transverse)
#Prepare outcome data for Colocalization - CRC

rm(list = ls())
set.seed(821)

#load libraries
library(data.table)
library(dplyr)
#remotes::install_github("mattlee821/functions")

data_SNPlist <- fread("analysis/006_colocalisation_gene-CRC/transverse/Exposure_GREM1-gene_transverse_SNPlist.txt" )
head(data_SNPlist)

#load and prepare outcome GWAS data
data_outcome <- fread("data/raw/GWAS/crc_combined_GECCO_allancestries.txt.gz",
                  header = TRUE , sep = "\t",
                  col.names = c("SNP","pval.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome",
                                "other_allele.outcome","beta.outcome", "se.outcome", "I2", "outcome"))
head(data_outcome)

#extracted SNPs from outcome GWAS by SNPs in the windows (using data_SNPlist from exposure data)
data_outcome <- data_outcome %>%
  filter(SNP %in% data_SNPlist$SNP)
head(data_outcome)
write.table(data_outcome, "analysis/006_colocalisation_gene-CRC/transverse/Outcome_GREM1-CRC_transverse_coloc.txt", sep="\t",row.names = FALSE, quote = FALSE, col.names = TRUE)
