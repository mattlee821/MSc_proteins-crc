#Prepare exposure data for Colocalization - Gene expression transverse
rm(list = ls())
set.seed(821)

#load libraries
library(data.table)
library(dplyr)

# Load exposure data
data_exposure <- fread("data/raw/GWAS/colon-transverse_allpairs_ENSG00000166923.10.txt")
# data_exposure <- data_exposure %>%     #?? DON'T Think we need this line, otherwisw we will have only 4 SNPs/data
#   filter(pval_nominal <= 5e-4)
data_exposure_ref<- fread("data/raw/GWAS/reference-chr15.txt") #load reference data
data_exposure <- data_exposure %>%
  left_join(data_exposure_ref, by = "variant_id") %>%
  rename(
    exposure = gene_id,
    SNP = rs_id_dbSNP151_GRCh38p7,           
    beta.exposure = slope,                   
    se.exposure = slope_se,                  
    effect_allele.exposure = alt,     
    other_allele.exposure = ref,      
    eaf.exposure = maf,  
    pval.exposure = pval_nominal,
    chr.exposure = chr, 
    pos.exposure = variant_pos) #NOT SURE how to get sample size for GTX8


#Filter SNP in dataset - use SNP from MR Gene expression-CRC (rs62002705)
data_SNP <- data_exposure %>%
  filter(SNP == "rs62002705") #?? Which SNP should I use - from MR gene expression-CRC? as this this gene expression data for exposure data. Diff SNP make a difference on number of data in each window

#Extract SNP chromosome
data_exposure <- data_exposure %>%
  filter(chr.exposure == data_SNP$chr.exposure)

## Mapping SNP - SNP from MR between gene experession-CRC analysis (Not SNP from protein-CRC)
# #IF use SNP from MR protein-CRC - It will make a difference on number of data in each window
# data_protein <- fread("data/raw/GWAS/plasma-proteins.txt")
# data_protein <- data_protein %>%
#   filter(exposure == "18878_15_GREM1_GREM1")

#Extract SNP position in 250kb, 500kb, 1m and 2m windows
window_250kb <- c(data_SNP$pos.exposure - 125000, data_SNP$pos.exposure + 125000)
data_window_250kb <- data_exposure %>%
  filter(pos.exposure >= window_250kb[1] & pos.exposure <= window_250kb[2])%>% 
  mutate(window="250kb")

window_500kb <- c(data_SNP$pos.exposure - 250000, data_SNP$pos.exposure + 250000)
data_window_500kb <- data_exposure %>%
  filter(pos.exposure >= window_500kb[1] & pos.exposure <= window_500kb[2])%>% 
  mutate(window="500kb")

window_1mb <- c(data_SNP$pos.exposure - 500000, data_SNP$pos.exposure + 500000)
data_window_1mb <- data_exposure %>%
  filter(pos.exposure >= window_1mb[1] & pos.exposure <= window_1mb[2])%>% 
  mutate(window="1mb")

window_2mb <- c(data_SNP$pos.exposure - 1000000, data_SNP$pos.exposure + 1000000)
data_window_2mb <- data_exposure %>%
  filter(pos.exposure >= window_2mb[1] & pos.exposure <= window_2mb[2]) %>% 
  mutate(window="2mb")

#Combine data frames and write the table for combinded windows
data_protein_window <- bind_rows(data_window_250kb, data_window_500kb, data_window_1mb, data_window_2mb)
write.table(data_protein_window, "analysis/005_colocalisation_gene-proteins/Exposure_GREM1-gene_window_coloc.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#Extract only SNP from exposures
data_SNPlist <- data_protein_window %>%
  select(SNP) %>% 
  unique() %>%
  as.data.frame()
write.table(data_SNPlist, "analysis/005_colocalisation_gene-proteins/Exposure_GREM1-gene_SNPlist.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
