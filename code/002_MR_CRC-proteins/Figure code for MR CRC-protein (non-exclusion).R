# Clear the environment
rm(list = ls())
set.seed(821)

# Load required libraries
library(data.table)
library(dplyr)
library(ggforestplot)
library(ggplot2)

# Load and prepare table_mr_crc_protein_nonexclusion
table_mr_crc_protein_nonexclusion <- fread("analysis/002_MR_CRC-proteins/table_mr_crc-protien_mr-nonexclusion.txt", 
                                           header = TRUE, 
                                           sep = "\t")

# Select relevant columns
table_mr_crc_protein_nonexclusion <- table_mr_crc_protein_nonexclusion %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 

# Display the first few rows of the data
head(table_mr_crc_protein_nonexclusion)

# Combine data into plot_data
plot_data <- table_mr_crc_protein_nonexclusion %>%
  mutate(
    Analysis = paste("CRC - Plasma GREM1 (", Method, ")", sep = "")  # Analysis name
  ) %>%
  select(Study = Exposure, Analysis, Method, Estimate, SE = Std_Error, P_Value)

# Convert SE and other numeric columns to numeric types
plot_data <- plot_data %>%
  mutate(
    Estimate = as.numeric(Estimate),
    SE = as.numeric(SE),
    P_Value = as.numeric(P_Value)
  )

# Adjust SE to positive values
plot_data$SE <- abs(plot_data$SE)

# Check unique method names in the dataset for consistency
unique_methods <- unique(plot_data$Method)
print(unique_methods)  # Print unique method names to verify

# Custom order for the Method factor to ensure the correct order in the plot
plot_data$Method <- factor(plot_data$Method, 
                           levels = c("Inverse variance weighted",  # Set IVW first
                                      "MR Egger", 
                                      "Weighted median", 
                                      "Simple mode"))

# Ensure that plot_data is ordered by Method for the forest plot
plot_data <- plot_data %>%
  arrange(Method)

# Plotting using ggforestplot
forest_plot <- ggforestplot::forestplot(
  df = plot_data,
  name = Analysis,  
  estimate = Estimate,
  se = SE,
  pvalue = P_Value,
  logodds = FALSE,
  psignif = 0.05,
  ci = 0.95
) 

# Add the section title above the plot using ggplot2 functionality
forest_plot <- forest_plot +
  ggtitle("CRC - Plasma Protein (including outliers)") +  
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  
    axis.text.y = element_text(size = 10, hjust = 0),
    strip.text = element_text(size = 10),
    legend.position = "none",  
    axis.title.x = element_text(size = 10, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    x = "Effect estimate (95% CI)",  
    y = ""  
  )

# Display the plot
print(forest_plot)
