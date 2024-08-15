#Prepare outcome data for Colocalization - CRC

rm(list = ls())
set.seed(821)

#load libraries
library(data.table)
library(dplyr)
#remotes::install_github("mattlee821/functions")

data_protein_SNPlist <- fread("analysis/003_colocalisation_proteins-CRC/GREM1-protein_SNPlist.txt" )
head(data_protein_SNPlist)

#load and prepare outcome GWAS data
data_crc_gwas <- fread("data/raw/GWAS/crc_combined_GECCO_allancestries.txt.gz",
                  header = TRUE , sep = "\t",
                  col.names = c("SNP","pval.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome",
                                "other_allele.outcome","beta.outcome", "se.outcome", "I2", "outcome"))
head(data_crc_gwas)

#extracted SNPs from outcome GWAS by SNPs in the windows (using data_protein_SNPlist from exposure data)
data_crc_coloc <- data_crc_gwas %>%
  filter(SNP %in% data_protein_SNPlist$SNP)
head(data_crc_coloc)
write.table(data_crc_coloc, "analysis/003_colocalisation_proteins-CRC/GREM1-CRC_coloc.txt", sep="\t",row.names = FALSE, quote = FALSE, col.names = TRUE)
