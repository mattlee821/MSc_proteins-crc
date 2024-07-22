rm(list = ls())
set.seed(821)

#load libraries
library(data.table)
library(dplyr)
remotes::install_github("mattlee821/functions")

#load and prepare exposure data for proteins
data_protein <- fread("data/raw/GWAS/plasma-proteins.txt")
data_protein <- data_protein %>%
  filter(exposure == "18878_15_GREM1_GREM1")

data_protein_gwas <- fread("data/raw/GWAS/18878_15_GREM1_GREM1.txt.gz.annotated.gz.exclusions.gz.alleles.gz",
                              header = FALSE, sep = "\t",
                              col.names = c("chr.exposure", "pos.exposure", "SNPID", "SNP", "A1", "A2",
                                            "beta.exposure", "pval.exposure", "min_log10_pval", "se.exposure",
                                            "samplesize.exposure", "ImpMAF", "exposure", "effect_allele.exposure",
                                            "other_allele.exposure", "eaf.exposure"))

#Filter SNP in GWAS dataset from plasma protein data
data_protein_SNP <- data_protein_gwas %>%
  filter(SNP == data_protein$SNP)

#Extract SNP chromosome
data_protein_chr <- data_protein_gwas %>%
  filter(chr.exposure == data_protein_SNP$chr.exposure)

#Shorten exposure names
data_protein_chr <- data_protein_chr %>%
  mutate(exposure = sub(".*/(.*)\\.txt\\.gz\\.unzipped", "\\1", exposure))
head(data_protein_chr)

#Extract SNP position in 250kb, 500kb, 1m and 2m windows
window_250kb <- c(data_protein_SNP$pos.exposure - 125000, data_protein_SNP$pos.exposure + 125000)
data_window_250kb <- data_protein_chr %>%
  filter(pos.exposure >= window_250kb[1] & pos.exposure <= window_250kb[2])%>% 
  mutate(window="250kb")

window_500kb <- c(data_protein_SNP$pos.exposure - 250000, data_protein_SNP$pos.exposure + 250000)
data_window_500kb <- data_protein_chr %>%
  filter(pos.exposure >= window_500kb[1] & pos.exposure <= window_500kb[2])%>% 
  mutate(window="500kb")

window_1mb <- c(data_protein_SNP$pos.exposure - 500000, data_protein_SNP$pos.exposure + 500000)
data_window_1mb <- data_protein_chr %>%
  filter(pos.exposure >= window_1mb[1] & pos.exposure <= window_1mb[2])%>% 
  mutate(window="1mb")

window_2mb <- c(data_protein_SNP$pos.exposure - 1000000, data_protein_SNP$pos.exposure + 1000000)
data_window_2mb <- data_protein_chr %>%
  filter(pos.exposure >= window_2mb[1] & pos.exposure <= window_2mb[2]) %>% 
  mutate(window="2mb")

#Combine data frames and write the table for combinded windows
data_protein_window <- bind_rows(data_window_250kb, data_window_500kb, data_window_1mb, data_window_2mb)
write.table(data_protein_window, "analysis/003_colocalisation_proteins-CRC/GREM1-protien_window_coloc.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

#Extract only SNP from exposures
data_protein_SNPlist <- data_protein_window %>%
  select(SNP) %>% 
  unique() %>%
  as.data.frame()
write.table(data_protein_SNPlist, "analysis/003_colocalisation_proteins-CRC/GREM1-protein_SNPlist.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
