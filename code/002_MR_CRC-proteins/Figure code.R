rm(list=ls())
set.seed(821)

# Load required libraries
library(ggforestplot)
library(dplyr)
library(data.table)

# Load and prepare table_mr_crc_protien_exclusion
table_mr_crc_protien_exclusion <- fread("analysis/002_MR_CRC-proteins/table_mr_crc_protien_exclusion.txt", header = TRUE, sep = "\t")
table_mr_crc_protien_exclusion <- table_mr_crc_protien_exclusion %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 

head(table_mr_crc_protien_exclusion)

# Load and prepare table_mr_protein_crc
table_mr_protein_crc <- fread("analysis/001_MR_proteins-CRC/table_mr_protein-crc.txt", header = TRUE, sep = "\t")
table_mr_protein_crc <- table_mr_protein_crc %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 

# Combine data into plot_data
plot_data <- bind_rows(
  table_mr_protein_crc %>%
    mutate(Analysis = "GREM1-CRC") %>%
    select(Study = Exposure, Analysis, Method, Estimate, SE = Std_Error, P_Value),
  table_mr_crc_protien_exclusion %>%
    mutate(Analysis = "CRC-GREM") %>%
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
  name = Method,
  estimate = Estimate,
  se = SE,
  pvalue = P_Value,
  logodds = FALSE,
  psignif = 0.05,
  ci = 0.95,
  colour = Analysis
) + theme(axis.text.y = element_text(size = 10))
