rm(list=ls())
set.seed(821)

# Load required libraries
library(ggforestplot)
library(dplyr)
library(data.table)
library(ggplot2)


# Load and prepare table_mr_protein_crc
table_mr_protein_crc <- fread("analysis/001_MR_proteins-CRC/table_mr_protein-crc.txt", header = TRUE, sep = "\t")
table_mr_protein_crc <- table_mr_protein_crc %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 

# Load and prepare table for MR_gene-crc_sigmoid 
table_mr_gene_crc_sigmoid  <- fread("analysis/007_MR_gene-CRC/001_ENSG00000166923-10_sigmoid_table_mr_gene-crc.txt", header = TRUE, sep = "\t")
table_mr_gene_crc_sigmoid <- table_mr_gene_crc_sigmoid %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 
head(table_mr_gene_crc_sigmoid)

# Load and prepare table for MR_gene-crc_transverse
table_mr_gene_crc_transverse <- fread("analysis/007_MR_gene-CRC/001_ENSG00000166923-10_transverse_table_mr_gene-crc.txt", header = TRUE, sep = "\t")
table_mr_gene_crc_transverse <- table_mr_gene_crc_transverse %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 
head(table_mr_gene_crc_transverse)

# Combine data into plot_data with an added 'Section' column
plot_data <- bind_rows(
  table_mr_protein_crc %>%
    mutate(
      Section = "Plasma Protein - CRC",
      Analysis = "Plasma GREM1 - CRC"  # No method mentioned
    ) %>%
    select(Study = Exposure, Section, Analysis, Method, Estimate, SE = Std_Error, P_Value),
  
  table_mr_gene_crc_sigmoid %>%
    mutate(
      Section = "Gene Expression - CRC",
      Analysis = "GREM1 gene expression (sigmoid) - CRC"
    ) %>%
    select(Study = Exposure, Section, Analysis, Method, Estimate, SE = Std_Error, P_Value),
  
  table_mr_gene_crc_transverse %>%
    mutate(
      Section = "Gene Expression - CRC",
      Analysis = "GREM1 gene expression (transverse) - CRC"
    ) %>%
    select(Study = Exposure, Section, Analysis, Method, Estimate, SE = Std_Error, P_Value)
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

# Set custom order for the Section factor
plot_data$Section <- factor(plot_data$Section, levels = c("Plasma Protein - CRC", "Gene Expression - CRC"))

# Sort data based on Method to ensure correct order in the plot
plot_data <- plot_data %>%
  arrange(Section)

# Plotting using ggforestplot with facet_wrap(~Section)
ggforestplot::forestplot(
  df = plot_data,
  name = Analysis,
  estimate = Estimate,
  se = SE,
  pvalue = P_Value,
  logodds = FALSE,
  psignif = 0.05,
  ci = 0.95,
  colour = Section  # Use Section for coloring each line
) + 
  facet_wrap(~Section, scales = "free_y", ncol = 1) +  # Arrange by Section in one column
  theme(
    axis.text.y = element_text(size = 10, hjust = 0),  # Align y-axis text (analysis names) to the left
    strip.text = element_text(size = 10),  # Adjust facet label text size
    legend.position = "none",              # Hide legend for clarity
    axis.title.x = element_text(size = 10, hjust = 0.2),  # Adjust x-axis label size and alignment
    plot.margin = margin(10, 10, 10, 10)   # Ensure there is enough margin for the plot
  ) +
  labs(
    x = "Effect estimate (95% CI)",  # Adjusted x-axis label
    y = ""
  )
