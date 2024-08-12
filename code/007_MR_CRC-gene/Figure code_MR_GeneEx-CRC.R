## Figure for MR between gene expression and CRC 
rm(list=ls())
set.seed(821)

# Load required libraries
library(ggforestplot)
library(dplyr)
library(data.table)

# Load and prepare table for MR_gene-crc_sigmoid 
table_mr_gene_crc_sigmoid  <- fread("analysis/006_MR_gene-CRC/001_ENSG00000166923-10_sigmoid_table_mr_gene-crc.txt", header = TRUE, sep = "\t")
table_mr_gene_crc_sigmoid <- table_mr_gene_crc_sigmoid %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 
head(table_mr_gene_crc_sigmoid)

# Load and prepare table for MR_gene-crc_transverse
table_mr_gene_crc_transverse <- fread("analysis/006_MR_gene-CRC/001_ENSG00000166923-10_transverse_table_mr_gene-crc.txt", header = TRUE, sep = "\t")
table_mr_gene_crc_transverse <- table_mr_gene_crc_transverse %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 
head(table_mr_gene_crc_transverse)

# Load and prepare table for MR_crc-gene_sigmoid 
table_mr_crc_gene_sigmoid <- fread("analysis/007_MR_CRC-gene/001_ENSG00000166923-10_sigmoid_table_mr_crc-gene.txt", header = TRUE, sep = "\t")
table_mr_crc_gene_sigmoid <- table_mr_crc_gene_sigmoid %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 
head(table_mr_crc_gene_sigmoid)

# Load and prepare table for MR_crc-gene_transverse
table_mr_crc_gene_transverse <- fread("analysis/007_MR_CRC-gene/001_ENSG00000166923-10_transverse_table_mr_crc-gene.txt", header = TRUE, sep = "\t")
table_mr_crc_gene_transverse <- table_mr_crc_gene_transverse %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 
head(table_mr_crc_gene_transverse)


# Combine data into plot_data
plot_data <- bind_rows(
  table_mr_gene_crc_sigmoid %>%
    mutate(Analysis = "Sigmoid: ENSG00000166923.10 - CRC") %>%
    select(Study = Exposure, Analysis, Method, Estimate, SE = Std_Error, P_Value),
  table_mr_gene_crc_transverse %>%
    mutate(Analysis = "Transverse: ENSG00000166923.10 - CRC") %>%
    select(Study = Exposure, Analysis, Method, Estimate, SE = Std_Error, P_Value), 
  table_mr_crc_gene_sigmoid %>%
    mutate(Analysis = "Sigmoid: CRC - ENSG00000166923.10") %>%
    select(Study = Exposure, Analysis, Method, Estimate, SE = Std_Error, P_Value),
  table_mr_crc_gene_transverse %>%
    mutate(Analysis = "Transverse: CRC - ENSG00000166923.10") %>%
    select(Study = Exposure, Analysis, Method, Estimate, SE = Std_Error, P_Value),
)

# Convert SE and other numeric columns to numeric types
plot_data <- plot_data %>%
  mutate(
    Estimate = as.numeric(Estimate),
    SE = as.numeric(SE),
    P_Value = as.numeric(P_Value)
  )

# Adjust SE to positive values
plot_data$SE <- abs(plot_data$SE)

# Custom order for the Method factor to ensure the correct order in the plot
plot_data$Method <- factor(plot_data$Method, levels = c("Wald ratio", "MR Egger", "Weighted median", "Inverse variance weighted", "Simple mode", "Weighted mode"))

# Plotting using ggforestplot
ggforestplot::forestplot(
  df = plot_data,
  name = Analysis,
  estimate = Estimate,
  se = SE,
  pvalue = P_Value,
  logodds = FALSE,
  psignif = 0.05,
  ci = 0.95,
  colour = Analysis
) + theme(axis.text.y = element_text(size = 10))
