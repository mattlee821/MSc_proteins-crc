# Clear the environment
rm(list = ls())
set.seed(821)

# Load required libraries
library(data.table)
library(dplyr)
library(ggforestplot)
library(ggplot2)

# Load and prepare table_mr_crc_protein_exclusion
table_mr_crc_protein_exclusion <- fread("analysis/002_MR_CRC-proteins/table_mr_crc-protein_mr-exclusion.txt", header = TRUE, sep = "\t")
table_mr_crc_protein_exclusion <- table_mr_crc_protein_exclusion %>%
  select(Exposure, Outcome, Method, Estimate, Std_Error, P_Value) 

# Display the first few rows of the data
head(table_mr_crc_protein_exclusion)

# Combine data into plot_data
plot_data <- table_mr_crc_protein_exclusion %>%
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

# Create the initial forest plot
forest_plot <- ggforestplot::forestplot(
  df = plot_data,
  name = Analysis,  # Use the Analysis column for the analysis names
  estimate = Estimate,
  se = SE,
  pvalue = P_Value,
  logodds = FALSE,
  psignif = 0.05,
  ci = 0.95
) 

# Customize the plot to change the line color to blue and remove points
forest_plot <- forest_plot +
  ggtitle("CRC - Plasma Protein") +  # Add the section title
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),  # Style the section title
    axis.text.y = element_text(size = 10, hjust = 0),
    strip.text = element_text(size = 10),
    legend.position = "none",  
    axis.title.x = element_text(size = 10, hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    x = "Effect estimate (95% CI)",  
    y = ""  # Remove y-axis label
  )

# Display the plot
print(forest_plot)
