
#Load packages
library(hakaiApi)
library(plotly)
library(dplyr)
library(stringr)
library(lubridate)
library(mgcv)
library(quantreg)
library(splines)
library(purrr)
library(tidyr)
library(here)

# Source functions
source(here("source", "qc_functions_clean_basic_qc_test.R"))

#Download data using the API

# Initialize the client
client <- Client$new()

#Setting up Query
endpoint <- "/eims/views/output/chlorophyll"
endpoint_hplc <- "/eims/views/output/hplc"
filter <- "site_id=QU39"
chl_url <- paste0("https://hecate.hakai.org/api", endpoint,"?limit=-1&", filter)
hplc_url <- paste0("https://hecate.hakai.org/api", endpoint_hplc,"?limit=-1&", filter)
data <- client$get(chl_url)
data_hplc <- client$get(hplc_url)




# --------------------- Basic QC Checks--------------------------------------

qc_results <- perform_comprehensive_chlorophyll_qc(data)



# View specific issues
qc_negative <- qc_results$negative_values          # Check 1: Negative chlorophyll/phaeo values
qc_acid_ratio <- qc_results$low_acid_ratios         # Check 2: Acid ratios < 1.01
qc_volume <- qc_results$volume_issues           # Check 3: Non-standard volumes (NEW)
qc_log <- qc_results$quality_log_issues      # Check 4: Unusually long quality logs (NEW)
qc_missing_fractions <- qc_results$missing_fractions       # Check 5: Missing size fractions
qc_rep_id <- qc_results$replicated_ids          # Check 6: Replicated hakai_id values
qc_cali_summary <- qc_results$calibration_summary     # Check 7: Calibration consistency


# NEW: View all quality logs in readable format
qc_results$quality_log_summary     # All quality logs with text wrapping
qc_results$quality_log_stats       # Statistics about quality log lengths

#Run basic QC checks and balances

#Need to add station so can apply beyond QU39.
#Would be good to separate dilutions.


# ---------------------- Time-series Checks ----------------------------------

#Run interactive plot to look at QC data incomparison to time series for each depth and filter type

# Method 1: GAM with moving window quantiles (recommended)
result_gam <- timeseries_qc(data, 
                            selected_depth = "5", 
                            selected_filter_type = "Bulk GF/F", 
                            investigation_year = 2024,
                            variability_method = "gam_quantiles",
                            window_days = 14)

# Method 2: Compare with original weekly method
# result_weekly <- timeseries_qc(data,
#                                selected_depth = "5",
#                                selected_filter_type = "Bulk GF/F",
#                                investigation_year = 2024,
#                                variability_method = "weekly")

create_monthly_boxplot_qc( 
  data = data,
  selected_depth = 0,
  selected_filter_type = "GF/F",
  investigation_year = 2024,
  variable = "acid_ratio",
  outlier_method = "IQR"
)

create_monthly_boxplot_qc( 
  data = data,
  selected_depth = 30,
  selected_filter_type = "Bulk GF/F",
  investigation_year = 2024,
  variable = "chla",
  outlier_method = "IQR"
)

# ---------------------- Time-series Checks ----------------------------------

# Run complete QC workflow
sf_qc_results <- run_complete_qc(data)

# Plot specific year with outlier highlighting (default)
plot_bulk_vs_sum_interactive(sf_qc_results, selected_year = 2024,
                             highlight_outliers = TRUE)

qc_metrics <- calculate_qc_metrics(sf_qc_results)

# -----------------------




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

# --------------------------------------------------------------------------

# Compare specific year to baseline
hplc_compare_plot(data, data_hplc, selected_year = 2022)


  

