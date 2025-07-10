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
library(mgcv)
library(quantreg)
library(splines)
library(purrr)
library(tidyr)

# Function to create improved chlorophyll QC plot with seasonal GAM
create_improved_chlorophyll_qc_plot <- function(data, selected_depth, selected_filter_type, investigation_year, 
                                                variability_method = "gam_quantiles", window_days = 14) {
  
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
      day_of_year = yday(date),
      chla_flag = as.logical(chla_flag)
    )
  
  # Debug information
  cat("Debug: Filtered data has", nrow(depth_data), "rows\n")
  cat("Debug: Year range:", min(depth_data$year, na.rm = TRUE), "to", max(depth_data$year, na.rm = TRUE), "\n")
  
  # Separate baseline and investigation data
  baseline_data <- depth_data %>%
    filter(year != investigation_year) %>%
    mutate(
      log_chla = log(chla + 0.01),  # Small offset to handle zeros
      cos1 = cos(2 * pi * day_of_year / 365.25),
      sin1 = sin(2 * pi * day_of_year / 365.25),
      cos2 = cos(4 * pi * day_of_year / 365.25),
      sin2 = sin(4 * pi * day_of_year / 365.25)
    )
  
  investigation_data <- depth_data %>%
    filter(year == investigation_year)
  
  if (nrow(baseline_data) < 20) {
    stop("Insufficient baseline data (< 20 points). Need more historical data.")
  }
  
  # Create prediction grid for smooth curves
  prediction_grid <- tibble(
    day_of_year = 1:365,
    cos1 = cos(2 * pi * day_of_year / 365.25),
    sin1 = sin(2 * pi * day_of_year / 365.25),
    cos2 = cos(4 * pi * day_of_year / 365.25),
    sin2 = sin(4 * pi * day_of_year / 365.25)
  )
  
  # Initialize method_name before conditional blocks
  method_name <- "Unknown Method"
  
  # Method 1: Seasonal GAM approach
  if (variability_method == "gam_quantiles") {
    
    # Fit seasonal GAM for median
    tryCatch({
      gam_model <- gam(log_chla ~ s(day_of_year, bs = "cc", k = 20), 
                       data = baseline_data)
      
      # Predict median values
      prediction_grid$pred_log_median <- predict(gam_model, newdata = prediction_grid)
      prediction_grid$pred_median <- exp(prediction_grid$pred_log_median) - 0.01
      
      # Calculate residuals for variability estimation
      baseline_data$residuals <- residuals(gam_model)
      baseline_data$fitted_log <- fitted(gam_model)
      
      # Moving window quantiles for variability
      calculate_moving_quantiles <- function(target_day, window_days = window_days) {
        # Handle year wrap-around
        if (target_day <= window_days) {
          day_range <- c((365 + target_day - window_days):365, 1:(target_day + window_days))
        } else if (target_day > (365 - window_days)) {
          day_range <- c((target_day - window_days):365, 1:(target_day + window_days - 365))
        } else {
          day_range <- (target_day - window_days):(target_day + window_days)
        }
        
        window_data <- baseline_data %>%
          filter(day_of_year %in% day_range)
        
        if (nrow(window_data) < 5) {
          return(list(q10 = NA, q90 = NA, q05 = NA, q95 = NA))
        }
        
        # Calculate quantiles of residuals
        residual_quantiles <- quantile(window_data$residuals, 
                                       probs = c(0.05, 0.10, 0.90, 0.95), 
                                       na.rm = TRUE)
        
        return(list(
          q05 = residual_quantiles[1],
          q10 = residual_quantiles[2], 
          q90 = residual_quantiles[3],
          q95 = residual_quantiles[4]
        ))
      }
      
      # Apply moving window to each day
      quantile_results <- map_dfr(1:365, ~{
        quants <- calculate_moving_quantiles(.x, window_days)
        tibble(
          day_of_year = .x,
          residual_q05 = quants$q05,
          residual_q10 = quants$q10,
          residual_q90 = quants$q90,
          residual_q95 = quants$q95
        )
      })
      
      # Combine with predictions
      prediction_grid <- prediction_grid %>%
        left_join(quantile_results, by = "day_of_year") %>%
        mutate(
          # Convert back to original scale
          lower_bound_90 = exp(pred_log_median + residual_q10) - 0.01,
          upper_bound_90 = exp(pred_log_median + residual_q90) - 0.01,
          lower_bound_95 = exp(pred_log_median + residual_q05) - 0.01,
          upper_bound_95 = exp(pred_log_median + residual_q95) - 0.01,
          
          # Ensure non-negative bounds
          lower_bound_90 = pmax(0, lower_bound_90),
          upper_bound_90 = pmax(0, upper_bound_90),
          lower_bound_95 = pmax(0, lower_bound_95),
          upper_bound_95 = pmax(0, upper_bound_95)
        )
      
      method_name <- paste0("GAM + Moving Window (±", window_days, " days)")
      
    }, error = function(e) {
      cat("GAM fitting failed:", e$message, "\nFalling back to harmonic regression\n")
      # Fallback to harmonic regression
      harm_model <- lm(log_chla ~ cos1 + sin1 + cos2 + sin2, data = baseline_data)
      prediction_grid$pred_log_median <- predict(harm_model, newdata = prediction_grid)
      prediction_grid$pred_median <- exp(prediction_grid$pred_log_median) - 0.01
      
      # Simple quantiles for fallback
      overall_q10 <- quantile(baseline_data$log_chla, 0.10, na.rm = TRUE)
      overall_q90 <- quantile(baseline_data$log_chla, 0.90, na.rm = TRUE)
      
      prediction_grid <- prediction_grid %>%
        mutate(
          lower_bound_90 = exp(pred_log_median + (overall_q10 - mean(baseline_data$log_chla, na.rm = TRUE))) - 0.01,
          upper_bound_90 = exp(pred_log_median + (overall_q90 - mean(baseline_data$log_chla, na.rm = TRUE))) - 0.01,
          lower_bound_95 = lower_bound_90 * 0.8,
          upper_bound_95 = upper_bound_90 * 1.2,
          lower_bound_90 = pmax(0, lower_bound_90),
          upper_bound_90 = pmax(0, upper_bound_90),
          lower_bound_95 = pmax(0, lower_bound_95),
          upper_bound_95 = pmax(0, upper_bound_95)
        )
      
      method_name <- "Harmonic Regression (GAM fallback)"
    })
    
  } else {
    # Fallback to original weekly method for comparison
    weekly_stats <- baseline_data %>%
      mutate(week = week(date)) %>%
      group_by(week) %>%
      summarise(
        median_chla = median(chla, na.rm = TRUE),
        q10 = quantile(chla, 0.10, na.rm = TRUE),
        q90 = quantile(chla, 0.90, na.rm = TRUE),
        q05 = quantile(chla, 0.05, na.rm = TRUE),
        q95 = quantile(chla, 0.95, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Convert to day of year for plotting
    prediction_grid <- prediction_grid %>%
      mutate(
        week = week(as.Date(paste("2020", day_of_year), format = "%Y %j"))
      ) %>%
      left_join(weekly_stats, by = "week") %>%
      mutate(
        pred_median = median_chla,
        lower_bound_90 = q10,
        upper_bound_90 = q90,
        lower_bound_95 = q05,
        upper_bound_95 = q95
      )
    
    method_name <- "Weekly Averages (Original)"
  }
  
  # Create the plot
  p <- plot_ly()
  
  # Add 95% confidence ribbon (lighter)
  p <- p %>%
    add_ribbons(
      data = prediction_grid,
      x = ~day_of_year,
      ymin = ~lower_bound_95,
      ymax = ~upper_bound_95,
      fillcolor = "rgba(59, 130, 246, 0.1)",
      line = list(color = "transparent"),
      name = "95% Range",
      hoverinfo = "text",
      text = ~paste("Day:", day_of_year, "<br>95% Range:", 
                    round(lower_bound_95, 3), "-", round(upper_bound_95, 3), "mg/m³")
    )
  
  # Add 90% confidence ribbon (darker)
  p <- p %>%
    add_ribbons(
      data = prediction_grid,
      x = ~day_of_year,
      ymin = ~lower_bound_90,
      ymax = ~upper_bound_90,
      fillcolor = "rgba(59, 130, 246, 0.2)",
      line = list(color = "transparent"),
      name = "90% Range",
      hoverinfo = "text",
      text = ~paste("Day:", day_of_year, "<br>90% Range:", 
                    round(lower_bound_90, 3), "-", round(upper_bound_90, 3), "mg/m³",
                    "<br>Median:", round(pred_median, 3), "mg/m³")
    )
  
  # Add median line
  p <- p %>%
    add_lines(
      data = prediction_grid,
      x = ~day_of_year,
      y = ~pred_median,
      line = list(color = "rgb(59, 130, 246)", width = 3),
      name = "Seasonal Median",
      hoverinfo = "text",
      text = ~paste("Day:", day_of_year, "<br>Predicted Median:", round(pred_median, 3), "mg/m³")
    )
  
  # Add investigation year points
  if (nrow(investigation_data) > 0) {
    # Normal points (not flagged)
    normal_points <- investigation_data %>% filter(!chla_flag | is.na(chla_flag))
    if (nrow(normal_points) > 0) {
      p <- p %>%
        add_markers(
          data = normal_points,
          x = ~day_of_year,
          y = ~chla,
          marker = list(color = "rgb(34, 197, 94)", size = 8),
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
          x = ~day_of_year,
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
        text = paste("Improved Chlorophyll QC - Depth:", selected_depth, "m -", selected_filter_type, "Filter<br>", 
                     "<sub>Method:", method_name, "</sub>"),
        font = list(size = 16, family = "Arial, sans-serif")
      ),
      xaxis = list(
        title = "Day of Year",
        showgrid = TRUE,
        gridcolor = "rgba(0,0,0,0.1)",
        range = c(1, 365)
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

# Updated summary function
generate_improved_summary_stats <- function(data, selected_depth, selected_filter_type, investigation_year) {
  
  investigation_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           year(as.Date(date)) == investigation_year,
           !is.na(chla)) %>%
    mutate(chla = as.numeric(chla),
           chla_flag = as.logical(chla_flag))
  
  baseline_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           year(as.Date(date)) != investigation_year,
           !is.na(chla)) %>%
    mutate(chla = as.numeric(chla))
  
  if (nrow(investigation_data) == 0) {
    return(list(
      total_samples = 0,
      flagged_samples = 0,
      mean_chla = NA,
      date_range = "No data",
      baseline_years = 0,
      baseline_samples = 0
    ))
  }
  
  stats <- list(
    total_samples = nrow(investigation_data),
    flagged_samples = sum(investigation_data$chla_flag, na.rm = TRUE),
    mean_chla = mean(investigation_data$chla, na.rm = TRUE),
    median_chla = median(investigation_data$chla, na.rm = TRUE),
    date_range = paste(min(investigation_data$date), "-", max(investigation_data$date)),
    baseline_years = length(unique(year(as.Date(baseline_data$date)))),
    baseline_samples = nrow(baseline_data)
  )
  
  return(stats)
}

# Updated main function
run_improved_chlorophyll_qc <- function(data, selected_depth = NULL, selected_filter_type = NULL, 
                                        investigation_year = NULL, variability_method = "gam_quantiles",
                                        window_days = 14) {
  
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
  plot <- create_improved_chlorophyll_qc_plot(data, selected_depth, selected_filter_type, 
                                              investigation_year, variability_method, window_days)
  
  # Generate summary statistics
  stats <- generate_improved_summary_stats(data, selected_depth, selected_filter_type, investigation_year)
  
  # Print summary
  cat("=== Improved Chlorophyll QC Summary ===\n")
  cat("Depth:", selected_depth, "m\n")
  cat("Filter Type:", selected_filter_type, "\n")
  cat("Investigation Year:", investigation_year, "\n")
  cat("Method:", variability_method, "\n")
  if (variability_method == "gam_quantiles") {
    cat("Window Size: ±", window_days, "days\n")
  }
  cat("Total Samples:", stats$total_samples, "\n")
  cat("Flagged Samples:", stats$flagged_samples, "\n")
  cat("Mean Chlorophyll:", round(stats$mean_chla, 3), "mg/m³\n")
  cat("Median Chlorophyll:", round(stats$median_chla, 3), "mg/m³\n")
  cat("Date Range:", stats$date_range, "\n")
  cat("Baseline Data:", stats$baseline_samples, "samples from", stats$baseline_years, "years\n")
  cat("=====================================\n\n")
  
  # Display the plot
  print(plot)
  
  return(list(plot = plot, stats = stats, options = options))
}

# Function to get available options (reused from original)
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

# ----------------------------------------------------------------------------




# Example usage:

# Method 1: GAM with moving window quantiles (recommended)
result_gam <- run_improved_chlorophyll_qc(data, 
                                          selected_depth = "5", 
                                          selected_filter_type = "3um", 
                                          investigation_year = 2024,
                                          variability_method = "gam_quantiles",
                                          window_days = 14)

# -------------------------------------------------------------------------


# Method 3: Compare with original weekly method
# result_weekly <- run_improved_chlorophyll_qc(data, 
#                                              selected_depth = "5", 
#                                              selected_filter_type = "Bulk GF/F", 
#                                              investigation_year = 2024,
#                                              variability_method = "weekly")
#--------------------------------------------------------------------------

# Function to sum size-fractions and compare to bulk
sum_size_fractions <- function(df) {
  
  # Check if required columns exist
  required_cols <- c("date", "site_id", "line_out_depth", "hakai_id", 
                     "collected", "filter_type", "chla", "chla_flag")
  missing_cols <- setdiff(required_cols, names(df))
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Define size classes, bulk, and bad flags
  size_classes <- c("20um", "3um", "GF/F")
  bulk_class <- "Bulk GF/F"
  bad_flags <- c("SVD")
  
  # Filter out bad quality data and split into size classes
  size_data <- df %>%
    filter(filter_type %in% size_classes) %>%
    filter(!chla_flag %in% bad_flags) %>%  # Remove bad quality data
    group_by(collected, date, site_id, line_out_depth) %>%
    summarise(
      n_size_fractions = n(),
      # Track which size classes are present
      size_classes_present = paste(sort(filter_type), collapse = ", "),
      # Only calculate sum if all 3 size fractions are present
      chla_sum = if_else(n() == 3, sum(chla, na.rm = TRUE), NA_real_),
      # Track if we have complete size fraction set
      complete_set = n() == 3,
      .groups = "drop"
    )
  
  # Get bulk measurements (excluding bad quality)
  bulk_data <- df %>%
    filter(filter_type == bulk_class) %>%
    filter(!chla_flag %in% bad_flags) %>%  # Remove bad quality bulk measurements
    select(collected, date, site_id, line_out_depth, 
           chla_bulk = chla, chla_bulk_flag = chla_flag, hakai_id)
  
  # Join size fraction sums with bulk measurements
  comparison_data <- size_data %>%
    left_join(bulk_data, by = c("collected", "date", "site_id", "line_out_depth"))
  
  # Calculate difference and ratio (only for complete sets with valid bulk)
  comparison_data <- comparison_data %>%
    mutate(
      year = as.numeric(format(as.Date(date), "%Y")),
      chla_diff = if_else(complete_set & !is.na(chla_bulk), 
                          chla_sum - chla_bulk, NA_real_),
      chla_ratio = if_else(complete_set & !is.na(chla_bulk) & chla_bulk > 0, 
                           chla_sum / chla_bulk, NA_real_),
      percent_diff = if_else(complete_set & !is.na(chla_bulk) & chla_bulk > 0, 
                             ((chla_sum - chla_bulk) / chla_bulk) * 100, NA_real_)
    )
  
  return(comparison_data)
}

# Example usage (assuming your data frame is called 'chla_data'):
result <- sum_size_fractions(data)
# 
# # View summary
# summary(result)
# 
# # Check for missing bulk measurements
# missing_bulk <- result %>% filter(is.na(chla_bulk))
# if(nrow(missing_bulk) > 0) {
#   cat("Warning:", nrow(missing_bulk), "size-fraction sets have no corresponding bulk measurement\n")
# }
# 
# # Check for incomplete size-fraction sets (should have 3 size classes)
# incomplete_sets <- result %>% filter(!complete_set)
# if(nrow(incomplete_sets) > 0) {
#   cat("Warning:", nrow(incomplete_sets), "sets have fewer than 3 good-quality size fractions\n")
#   print(incomplete_sets %>% select(collected, site_id, n_size_fractions, size_classes_present))
# }
# 
# # Check how many samples have valid comparisons
# valid_comparisons <- result %>% filter(!is.na(chla_diff))
# cat("Valid bulk vs sum comparisons:", nrow(valid_comparisons), "out of", nrow(result), "total sample sets\n")

# Interactive visualization function for QC
# Interactive visualization function for QC
plot_bulk_vs_sum_interactive <- function(comparison_data, selected_year = NULL) {
  library(plotly)
  library(ggplot2)
  
  # Only plot valid comparisons (complete sets with good bulk data)
  plot_data <- comparison_data %>%
    filter(!is.na(chla_diff) & complete_set)
  
  if(nrow(plot_data) == 0) {
    warning("No valid comparisons available for plotting")
    return(NULL)
  }
  
  # Create year grouping for coloring
  if(!is.null(selected_year)) {
    plot_data <- plot_data %>%
      mutate(year_group = if_else(year == selected_year, 
                                  paste("Year", selected_year), 
                                  "Other years"))
  } else {
    plot_data$year_group <- "All data"
  }
  
  # Calculate regression for full dataset
  full_lm <- lm(chla_sum ~ chla_bulk, data = plot_data)
  
  # Calculate regression for selected year if specified
  year_lm <- NULL
  if(!is.null(selected_year)) {
    year_data <- plot_data %>% filter(year == selected_year)
    if(nrow(year_data) >= 3) {
      year_lm <- lm(chla_sum ~ chla_bulk, data = year_data)
    }
  }
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = chla_bulk, y = chla_sum)) +
    # Add confidence band for full dataset regression
    geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue", alpha = 0.3,
                linetype = "solid", size = 1) +
    # Add 1:1 reference line
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(
      x = "Bulk Chlorophyll-a (mg/m³)",
      y = "Sum of Size Fractions (mg/m³)",
      title = "Interactive Bulk vs Sum of Size Fractions QC Plot"
    ) +
    theme_minimal()
  
  # Add year-specific regression line if available
  if(!is.null(year_lm)) {
    year_data <- plot_data %>% filter(year == selected_year)
    p <- p + geom_smooth(data = year_data, method = "lm", se = FALSE, 
                         color = "darkred", linetype = "solid", size = 1.2)
  }
  
  # Add points with hover information
  p <- p + geom_point(aes(text = paste("Hakai ID:", hakai_id,
                                       "<br>Date:", date,
                                       "<br>Site:", site_id,
                                       "<br>Depth:", line_out_depth, "m",
                                       "<br>Year:", year,
                                       "<br>Bulk:", round(chla_bulk, 3),
                                       "<br>Sum:", round(chla_sum, 3),
                                       "<br>% Diff:", round(percent_diff, 1)),
                          color = if(!is.null(selected_year)) year_group else NULL),
                      alpha = 0.9)
  
  # Add colors if year is selected
  if(!is.null(selected_year)) {
    p <- p + scale_color_manual(values = c("grey50", "darkred")) +
      labs(color = "Year Group")
  } else {
    p <- p + scale_color_manual(values = "steelblue")
  }
  
  # Add regression info to subtitle
  r2_full <- summary(full_lm)$r.squared
  slope_full <- coef(full_lm)[2]
  
  subtitle_text <- paste0("Full dataset: R² = ", round(r2_full, 3), 
                          ", Slope = ", round(slope_full, 3))
  
  if(!is.null(year_lm)) {
    r2_year <- summary(year_lm)$r.squared
    slope_year <- coef(year_lm)[2]
    subtitle_text <- paste0(subtitle_text, " | Year ", selected_year, 
                            ": R² = ", round(r2_year, 3), 
                            ", Slope = ", round(slope_year, 3))
  }
  
  p <- p + labs(subtitle = subtitle_text)
  
  # Convert to interactive plot with better dimensions
  ggplotly(p, tooltip = "text", height = 600, width = 700)
}

# QC metrics function
calculate_qc_metrics <- function(comparison_data) {
  # Filter to valid comparisons
  valid_data <- comparison_data %>%
    filter(!is.na(chla_diff) & complete_set)
  
  if(nrow(valid_data) == 0) {
    warning("No valid data for QC metrics")
    return(NULL)
  }
  
  # Overall metrics
  overall_metrics <- valid_data %>%
    summarise(
      n_samples = n(),
      mean_ratio = mean(chla_ratio, na.rm = TRUE),
      sd_ratio = sd(chla_ratio, na.rm = TRUE),
      median_ratio = median(chla_ratio, na.rm = TRUE),
      mad_ratio = mad(chla_ratio, na.rm = TRUE),
      mean_percent_diff = mean(percent_diff, na.rm = TRUE),
      sd_percent_diff = sd(percent_diff, na.rm = TRUE),
      # Regression metrics
      r_squared = cor(chla_bulk, chla_sum, use = "complete.obs")^2,
      slope = coef(lm(chla_sum ~ chla_bulk, data = .))[2],
      intercept = coef(lm(chla_sum ~ chla_bulk, data = .))[1],
      rmse = sqrt(mean((chla_sum - chla_bulk)^2, na.rm = TRUE))
    ) %>%
    mutate(period = "Full time series")
  
  # Yearly metrics
  yearly_metrics <- valid_data %>%
    group_by(year) %>%
    summarise(
      n_samples = n(),
      mean_ratio = mean(chla_ratio, na.rm = TRUE),
      sd_ratio = sd(chla_ratio, na.rm = TRUE),
      median_ratio = median(chla_ratio, na.rm = TRUE),
      mad_ratio = mad(chla_ratio, na.rm = TRUE),
      mean_percent_diff = mean(percent_diff, na.rm = TRUE),
      sd_percent_diff = sd(percent_diff, na.rm = TRUE),
      # Regression metrics (if enough data)
      r_squared = if_else(n() >= 3, cor(chla_bulk, chla_sum, use = "complete.obs")^2, NA_real_),
      slope = if_else(n() >= 3, coef(lm(chla_sum ~ chla_bulk))[2], NA_real_),
      intercept = if_else(n() >= 3, coef(lm(chla_sum ~ chla_bulk))[1], NA_real_),
      rmse = sqrt(mean((chla_sum - chla_bulk)^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(period = paste("Year", year))
  
  # Combine results
  all_metrics <- bind_rows(overall_metrics, yearly_metrics)
  
  return(all_metrics)
}

# QC comparison function
compare_year_qc <- function(comparison_data, selected_year) {
  metrics <- calculate_qc_metrics(comparison_data)
  
  if(is.null(metrics)) return(NULL)
  
  overall <- metrics %>% filter(period == "Full time series")
  year_data <- metrics %>% filter(period == paste("Year", selected_year))
  
  if(nrow(year_data) == 0) {
    cat("No data available for year", selected_year, "\n")
    return(NULL)
  }
  
  cat("=== QC COMPARISON FOR YEAR", selected_year, "===\n\n")
  
  cat("Sample sizes:\n")
  cat("  Full time series:", overall$n_samples, "samples\n")
  cat("  Year", selected_year, ":", year_data$n_samples, "samples\n\n")
  
  cat("Bulk vs Sum Ratio (ideally ~1.0):\n")
  cat("  Full series - Mean:", round(overall$mean_ratio, 3), "±", round(overall$sd_ratio, 3), "\n")
  cat("  Year", selected_year, "- Mean:", round(year_data$mean_ratio, 3), "±", round(year_data$sd_ratio, 3), "\n")
  cat("  Difference from overall:", round(year_data$mean_ratio - overall$mean_ratio, 3), "\n\n")
  
  cat("Percent Difference:\n")
  cat("  Full series - Mean:", round(overall$mean_percent_diff, 1), "%\n")
  cat("  Year", selected_year, "- Mean:", round(year_data$mean_percent_diff, 1), "%\n\n")
  
  cat("Regression Metrics:\n")
  cat("  Full series - R²:", round(overall$r_squared, 3), "Slope:", round(overall$slope, 3), "\n")
  if(!is.na(year_data$r_squared)) {
    cat("  Year", selected_year, "- R²:", round(year_data$r_squared, 3), "Slope:", round(year_data$slope, 3), "\n")
  }
  
  return(metrics)
}

# Example plotting:
# plot_bulk_vs_sum(result)

# NEW USAGE EXAMPLES:
# 
# # 1. Interactive plot with hover information
plot_bulk_vs_sum_interactive(result)
# 
# 2. Interactive plot highlighting a specific year
plot_bulk_vs_sum_interactive(result, selected_year = 2024)
# 
# # 3. Calculate QC metrics for all years
qc_metrics <- calculate_qc_metrics(result)

test <- compare_year_qc(result, 2022)









