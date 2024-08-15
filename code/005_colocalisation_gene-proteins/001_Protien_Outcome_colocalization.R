#005_colocalisation_gene-proteins
#Prepare outcome data for Colocalization - Plasma Protein (GREM1)

rm(list = ls())
set.seed(821)

#load libraries
library(data.table)
library(dplyr)

#Load SNP list from gene expression (exposure) data
data_SNPlist <- fread("analysis/005_colocalisation_gene-proteins/Exposure_GREM1-gene_SNPlist.txt" )
head(data_SNPlist)

#load and prepare outcome 
data_outcome <- fread("data/raw/GWAS/18878_15_GREM1_GREM1.txt.gz.annotated.gz.exclusions.gz.alleles.gz",
                           header = FALSE, sep = "\t",
                           col.names = c("chr.exposure", "pos.outcome", "SNPID", "SNP", "A1", "A2",
                                         "beta.outcome", "pval.outcome", "min_log10_pval", "se.outcome",
                                         "samplesize.outcome", "ImpMAF", "outcome", "effect_allele.outcome",
                                         "other_allele.outcome", "eaf.outcome"))

head(data_outcome)

#extracted SNPs from outcome GWAS by SNPs in the windows (using data_SNPlist from exposure data)
data_outcome <- data_outcome %>%
  filter(SNP %in% data_SNPlist$SNP)

#Shorten exposure names
data_outcome <- data_outcome %>%
  mutate(outcome = sub(".*/(.*)\\.txt\\.gz\\.unzipped", "\\1", outcome))
head(data_outcome)
write.table(data_outcome, "analysis/005_colocalisation_gene-proteins/Outcome_GREM1-protein_coloc.txt", sep="\t",row.names = FALSE, quote = FALSE, col.names = TRUE)
