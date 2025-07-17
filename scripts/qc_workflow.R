
#Load packages
library(hakaiApi)
library(plotly)
library(dplyr)
library(lubridate)
library(mgcv)
library(quantreg)
library(splines)
library(purrr)
library(tidyr)
library(here)

# Source functions
source(here("source", "qc_functions.R"))

#Download data using the API

# Initialize the client
client <- Client$new()

#Setting up Query
endpoint <- "/eims/views/output/chlorophyll"
filter <- "site_id=QU39"
chl_url <- paste0("https://hecate.hakai.org/api", endpoint,"?limit=-1&", filter)
data <- client$get(chl_url)


# --------------------- Basic QC Checks--------------------------------------

#Run basic QC checks and balances

#Looks for a variety of errors included negative values, low acid ratios, missing size fractions, duplicated IDs and incorrect calibrations
qc_results <- perform_basic_chlorophyll_qc_checks(data)
# 
# # Access specific results:
negative_values <- qc_results$negative_values
low_acid_ratios <- qc_results$low_acid_ratios
missing_fractions <- qc_results$missing_fractions
replicated_ids <- qc_results$replicated_ids
calibration_summary <- qc_results$calibration_summary
calibration_date_issues <- qc_results$calibration_date_issues

# -----------------------------------

#Run interactive plot to look at QC data incomparison to time series for each depth and filter type

# Method 1: GAM with moving window quantiles (recommended)
result_gam <- run_improved_chlorophyll_qc(data, 
                                          selected_depth = "10", 
                                          selected_filter_type = "20um", 
                                          investigation_year = 2024,
                                          variability_method = "gam_quantiles",
                                          window_days = 14)

# Method 3: Compare with original weekly method
result_weekly <- run_improved_chlorophyll_qc(data,
                                             selected_depth = "5",
                                             selected_filter_type = "Bulk GF/F",
                                             investigation_year = 2024,
                                             variability_method = "weekly")

result_monthly_boxplot <- create_monthly_boxplot_qc( 
  data = data,
  selected_depth = 5,
  selected_filter_type = "GF/F",
  investigation_year = 2023,
  variable = "acid_ratio",
  outlier_method = "IQR"
)

# -----------------------------------------------------

# Run complete QC workflow
sf_qc_results <- run_complete_qc(data)

# Plot specific year with outlier highlighting (default)
plot_bulk_vs_sum_interactive(sf_qc_results, selected_year = 2024,
                             highlight_outliers = TRUE)

# -----------------------

create_monthly_boxplot_qc(data,
                                    selected_depth = 100,
                                    selected_filter_type = "Bulk GF/F",
                                    investigation_year = 2024,
                                    variable = "acid_ratio",
                                    show_outliers = TRUE,
                                    outlier_method = "iqr")




# Required:
#   
# data - Your dataset containing the phytoplankton/chlorophyll data
# 
# Optional parameters (with defaults):
#   
# selected_depth - Water depth in meters (defaults to first available depth)
# selected_filter_type - Filter type (defaults to "GF/F" if available, otherwise first available)
# investigation_year - Year to investigate for outliers (defaults to most recent year)
# variable - Either "chla" or "acid_ratio" (defaults to "chla")
# show_outliers - TRUE/FALSE (defaults to TRUE)
# outlier_method - Either "iqr" or "quantile" (defaults to "iqr")