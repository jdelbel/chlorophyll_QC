library(hakaiApi)

# Initialize the client
client <- Client$new()

#Setting up Query
endpoint <- "/eims/views/output/chlorophyll"
filter <- "site_id=QU39"
chl_url <- paste0("https://hecate.hakai.org/api", endpoint,"?limit=-1&", filter)
data <- client$get(chl_url)


library(plotly)
library(dplyr)
library(lubridate)
library(RColorBrewer)

# Function to create chlorophyll QC plot
create_chlorophyll_qc_plot <- function(data, selected_depth, selected_filter_type, investigation_year, variability_method = "percentile_10_90") {
  
  # Filter data for selected depth and filter type
  depth_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           !is.na(chla),
           !is.na(date)) %>%
    mutate(
      date = as.Date(date),
      chla = as.numeric(chla),
      year = year(date),
      week = week(date),
      chla_flag = as.logical(chla_flag)
    )
  
  # Debug information
  cat("Debug: Filtered data has", nrow(depth_data), "rows\n")
  cat("Debug: Year range:", min(depth_data$year, na.rm = TRUE), "to", max(depth_data$year, na.rm = TRUE), "\n")
  cat("Debug: Selected depth:", selected_depth, "type:", class(selected_depth), "\n")
  cat("Debug: Available depths:", unique(data$line_out_depth), "\n")
  
  # Separate baseline and investigation data
  baseline_data <- depth_data %>%
    filter(year != investigation_year)
  
  investigation_data <- depth_data %>%
    filter(year == investigation_year)
  
  # Calculate weekly statistics for baseline using multiple approaches
  weekly_stats <- baseline_data %>%
    group_by(week) %>%
    summarise(
      mean_chla = mean(chla, na.rm = TRUE),
      median_chla = median(chla, na.rm = TRUE),
      sd_chla = sd(chla, na.rm = TRUE),
      q25 = quantile(chla, 0.25, na.rm = TRUE),
      q75 = quantile(chla, 0.75, na.rm = TRUE),
      q10 = quantile(chla, 0.10, na.rm = TRUE),
      q90 = quantile(chla, 0.90, na.rm = TRUE),
      q05 = quantile(chla, 0.05, na.rm = TRUE),
      q95 = quantile(chla, 0.95, na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      # Option 1: Interquartile range (IQR) - most robust for skewed data
      iqr = q75 - q25,
      iqr_lower = q25 - 1.5 * iqr,
      iqr_upper = q75 + 1.5 * iqr,
      
      # Option 2: 10th-90th percentile range
      p10_90_lower = q10,
      p10_90_upper = q90,
      
      # Option 3: 5th-95th percentile range
      p5_95_lower = q05,
      p5_95_upper = q95,
      
      # Option 4: Median ± 2*MAD (Median Absolute Deviation)
      mad_chla = 1.4826 * median(abs(baseline_data$chla[baseline_data$week == week] - median_chla), na.rm = TRUE),
      mad_lower = pmax(0, median_chla - 2 * mad_chla),
      mad_upper = median_chla + 2 * mad_chla,
      
      # Ensure non-negative lower bounds
      iqr_lower = pmax(0, iqr_lower),
      p10_90_lower = pmax(0, p10_90_lower),
      p5_95_lower = pmax(0, p5_95_lower),
      
      # Set bounds based on method
      upper_bound = case_when(
        variability_method == "sd" ~ mean_chla + sd_chla,
        variability_method == "sd_2" ~ mean_chla + 2 * sd_chla,
        variability_method == "iqr" ~ iqr_upper,
        variability_method == "percentile_10_90" ~ p10_90_upper,
        variability_method == "percentile_5_95" ~ p5_95_upper,
        variability_method == "mad" ~ mad_upper,
        TRUE ~ p10_90_upper
      ),
      lower_bound = case_when(
        variability_method == "sd" ~ pmax(0, mean_chla - sd_chla),
        variability_method == "sd_2" ~ pmax(0, mean_chla - 2 * sd_chla),
        variability_method == "iqr" ~ iqr_lower,
        variability_method == "percentile_10_90" ~ p10_90_lower,
        variability_method == "percentile_5_95" ~ p5_95_lower,
        variability_method == "mad" ~ mad_lower,
        TRUE ~ p10_90_lower
      ),
      method_name = case_when(
        variability_method == "sd" ~ "±1 SD",
        variability_method == "sd_2" ~ "±2 SD",
        variability_method == "iqr" ~ "IQR ± 1.5×IQR",
        variability_method == "percentile_10_90" ~ "10th-90th %ile",
        variability_method == "percentile_5_95" ~ "5th-95th %ile",
        variability_method == "mad" ~ "Median ± 2×MAD",
        TRUE ~ "10th-90th %ile"
      )
    )
  
  # Create the plot
  p <- plot_ly()
  
  # Add variability ribbon (using method-specific name)
  p <- p %>%
    add_ribbons(
      data = weekly_stats,
      x = ~week,
      ymin = ~lower_bound,
      ymax = ~upper_bound,
      fillcolor = "rgba(59, 130, 246, 0.2)",
      line = list(color = "transparent"),
      name = ~paste("Range:", method_name[1]),
      hoverinfo = "text",
      text = ~paste("Week:", week, "<br>Range (", method_name, "):", 
                    round(lower_bound, 3), "-", round(upper_bound, 3), "mg/m³",
                    "<br>Median:", round(median_chla, 3), "mg/m³")
    )
  
  # Add median line (more robust than mean for skewed data)
  p <- p %>%
    add_lines(
      data = weekly_stats,
      x = ~week,
      y = ~median_chla,
      line = list(color = "rgb(59, 130, 256)", width = 2),
      name = "Historical Median",
      hoverinfo = "text",
      text = ~paste("Week:", week, "<br>Median:", round(median_chla, 3), "mg/m³",
                    "<br>Mean:", round(mean_chla, 3), "mg/m³",
                    "<br>n =", n, "samples")
    )
  
  # Add investigation year points
  if (nrow(investigation_data) > 0) {
    # Normal points (not flagged)
    normal_points <- investigation_data %>% filter(!chla_flag | is.na(chla_flag))
    if (nrow(normal_points) > 0) {
      p <- p %>%
        add_markers(
          data = normal_points,
          x = ~week,
          y = ~chla,
          marker = list(color = "rgb(59, 130, 246)", size = 8),
          name = paste(investigation_year, "- Normal"),
          hoverinfo = "text",
          text = ~paste("Date:", date, "<br>Hakai ID:", hakai_id,
                        "<br>Chlorophyll:", round(chla, 3), "mg/m³",
                        "<br>Depth:", line_out_depth, "m")
        )
    }
    
    # Flagged points
    flagged_points <- investigation_data %>% filter(chla_flag == TRUE)
    if (nrow(flagged_points) > 0) {
      p <- p %>%
        add_markers(
          data = flagged_points,
          x = ~week,
          y = ~chla,
          marker = list(color = "rgb(239, 68, 68)", size = 8),
          name = paste(investigation_year, "- Flagged"),
          hoverinfo = "text",
          text = ~paste("Date:", date, "<br>Hakai ID:", hakai_id,
                        "<br>Chlorophyll:", round(chla, 3), "mg/m³",
                        "<br>Depth:", line_out_depth, "m",
                        "<br>⚠ FLAGGED")
        )
    }
  }
  
  # Layout
  p <- p %>%
    layout(
      title = list(
        text = paste("Chlorophyll QC - Depth:", selected_depth, "m -", selected_filter_type, "Filter"),
        font = list(size = 16, family = "Arial, sans-serif")
      ),
      xaxis = list(
        title = "Week of Year",
        showgrid = TRUE,
        gridcolor = "rgba(0,0,0,0.1)"
      ),
      yaxis = list(
        title = "Chlorophyll (mg/m³)",
        showgrid = TRUE,
        gridcolor = "rgba(0,0,0,0.1)"
      ),
      hovermode = "closest",
      legend = list(
        x = 0.01,
        y = 0.99,
        bgcolor = "rgba(255,255,255,0.8)",
        bordercolor = "rgba(0,0,0,0.2)",
        borderwidth = 1
      ),
      plot_bgcolor = "white",
      paper_bgcolor = "white"
    )
  
  return(p)
}

# Function to generate summary statistics (FIXED)
generate_summary_stats <- function(data, selected_depth, selected_filter_type, investigation_year) {
  
  investigation_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           year(as.Date(date)) == investigation_year,
           !is.na(chla)) %>%
    mutate(chla = as.numeric(chla),
           chla_flag = as.logical(chla_flag))
  
  if (nrow(investigation_data) == 0) {
    return(list(
      total_samples = 0,
      flagged_samples = 0,
      mean_chla = NA,
      date_range = "No data"
    ))
  }
  
  stats <- list(
    total_samples = nrow(investigation_data),
    flagged_samples = sum(investigation_data$chla_flag, na.rm = TRUE),
    mean_chla = mean(investigation_data$chla, na.rm = TRUE),
    date_range = paste(min(investigation_data$date), "-", max(investigation_data$date))
  )
  
  return(stats)
}

# Function to get available options
get_available_options <- function(data) {
  depths <- data %>%
    filter(!is.na(line_out_depth)) %>%
    distinct(line_out_depth) %>%
    arrange(as.numeric(line_out_depth)) %>%
    pull(line_out_depth)
  
  years <- data %>%
    filter(!is.na(date)) %>%
    mutate(year = year(as.Date(date))) %>%
    distinct(year) %>%
    arrange(desc(year)) %>%
    pull(year)
  
  filter_types <- data %>%
    filter(!is.na(filter_type)) %>%
    distinct(filter_type) %>%
    pull(filter_type)
  
  return(list(depths = depths, years = years, filter_types = filter_types))
}

# Example usage function (FIXED)
run_chlorophyll_qc <- function(data, selected_depth = NULL, selected_filter_type = NULL, 
                               investigation_year = NULL, variability_method = "percentile_10_90") {
  
  # Get available options
  options <- get_available_options(data)
  
  # Use defaults if not provided
  if (is.null(selected_depth)) {
    selected_depth <- options$depths[1]
  }
  if (is.null(selected_filter_type)) {
    selected_filter_type <- options$filter_types[1]
  }
  if (is.null(investigation_year)) {
    investigation_year <- options$years[1]
  }
  
  # Create the plot
  plot <- create_chlorophyll_qc_plot(data, selected_depth, selected_filter_type, investigation_year, variability_method)
  
  # Generate summary statistics (FIXED - now passing all 4 arguments)
  stats <- generate_summary_stats(data, selected_depth, selected_filter_type, investigation_year)
  
  # Print summary
  cat("=== Chlorophyll QC Summary ===\n")
  cat("Depth:", selected_depth, "m\n")
  cat("Filter Type:", selected_filter_type, "\n")
  cat("Investigation Year:", investigation_year, "\n")
  cat("Variability Method:", variability_method, "\n")
  cat("Total Samples:", stats$total_samples, "\n")
  cat("Flagged Samples:", stats$flagged_samples, "\n")
  cat("Mean Chlorophyll:", round(stats$mean_chla, 3), "mg/m³\n")
  cat("Date Range:", stats$date_range, "\n")
  cat("===============================\n\n")
  
  # Display the plot
  print(plot)
  
  return(list(plot = plot, stats = stats, options = options))
}

# Check available options
options <- get_available_options(data)
options$filter_types  # See available filter types: 'Bulk GF/F', '20um', '3um', 'GF/F'

# Test different variability methods:

# Method 1: 10th-90th percentile (recommended for phytoplankton)
result_p10_90 <- run_chlorophyll_qc(data, 
                                    selected_depth = "5", 
                                    selected_filter_type = "GF/F", 
                                    investigation_year = 2021,
                                    variability_method = "percentile_10_90")

# # Method 2: 5th-95th percentile (even more inclusive)
# result_p5_95 <- run_chlorophyll_qc(data, 
#                                    selected_depth = "5", 
#                                    selected_filter_type = "Bulk GF/F", 
#                                    investigation_year = 2024,
#                                    variability_method = "percentile_5_95")
# 
# # Method 3: IQR + 1.5*IQR (traditional outlier detection)
# result_iqr <- run_chlorophyll_qc(data, 
#                                  selected_depth = "5", 
#                                  selected_filter_type = "Bulk GF/F", 
#                                  investigation_year = 2024,
#                                  variability_method = "iqr")
# 
# # Method 4: Median ± 2*MAD (robust to extreme outliers)
# result_mad <- run_chlorophyll_qc(data, 
#                                  selected_depth = "5", 
#                                  selected_filter_type = "Bulk GF/F", 
#                                  investigation_year = 2024,
#                                  variability_method = "mad")
# 
# # Method 5: Mean ± 2*SD (more inclusive than ±1 SD)
# result_2sd <- run_chlorophyll_qc(data, 
#                                  selected_depth = "5", 
#                                  selected_filter_type = "Bulk GF/F", 
#                                  investigation_year = 2024,
#                                  variability_method = "sd_2")