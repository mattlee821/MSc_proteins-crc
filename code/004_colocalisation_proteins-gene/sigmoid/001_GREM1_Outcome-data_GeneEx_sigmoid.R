#004-2_colocalisation_proteins-gene sigmoid
#Prepare outcome data for Colocalization - Gene expression sigmoid

rm(list = ls())
set.seed(821)

#load libraries
library(data.table)
library(dplyr)

data_SNPlist <- fread("analysis/004_colocalisation_proteins-gene/sigmoid/Exposure_GREM1-protein-SNPlist.txt" )
head(data_SNPlist)

#load and prepare outcome - Gene expression sigmoid
data_outcome <- fread("data/raw/GWAS/colon-sigmoid_allpairs_ENSG00000166923.10.txt")
# data_outcome <- data_outcome %>%    ###QUESTION: DON'T Think we need this line (to cut off by p-value), otherwisw we will have only 4 SNPs/data. Is it correct? (I didn't run it)
#   filter(pval_nominal <= 5e-4)
data_outcome_ref<- fread("data/raw/GWAS/reference-chr15.txt") #load reference data
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
    pval.outcome = pval_nominal,
    chr.outcome = chr, 
    pos.outcome = variant_pos) 
head(data_outcome)

#extracted SNPs from outcome GWAS by SNPs in the windows (using data_SNPlist from exposure data)
data_outcome <- data_outcome %>%
  filter(SNP %in% data_SNPlist$SNP)
head(data_outcome)
write.table(data_outcome, "analysis/004_colocalisation_proteins-gene/sigmoid/Outcome_GREM1-gene_sigmoid_coloc.txt", sep="\t",row.names = FALSE, quote = FALSE, col.names = TRUE)
