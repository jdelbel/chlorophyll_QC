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
# result <- sum_size_fractions(data)
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

# -----------------------------------------------------------------------------


# Function to perform basic chlorophyll QC checks
perform_basic_chlorophyll_qc_checks <- function(data) {
  
  # Convert date to Date class if it isn't already
  data <- data %>%
    mutate(date = as.Date(date))
  
  # Check 1: Records with negative chlorophyll and/or phaeopigments
  cat("=== CHECK 1: NEGATIVE VALUES ===\n")
  negative_values <- data %>%
    filter(chla < 0 | phaeo < 0 | (!is.na(chla) & !is.na(phaeo) & (chla < 0 | phaeo < 0))) %>%
    select(hakai_id, site_id, line_out_depth, date, chla, phaeo, 
           any_of(c("chla_flag", "phaeo_flag"))) %>%
    mutate(
      chla = round(chla, 2),
      phaeo = round(phaeo, 2)
    ) %>%
    arrange(date, site_id, line_out_depth)
  
  if (nrow(negative_values) > 0) {
    cat("Found", nrow(negative_values), "records with negative chlorophyll and/or phaeopigments:\n")
    print(negative_values)
    
    # Create issues directory if it doesn't exist
    if (!dir.exists(here("issues"))) {
      dir.create(here("issues"), recursive = TRUE)
    }
    
    # Save to CSV
    csv_filename <- here("issues", paste0("negative_values_", Sys.Date(), ".csv"))
    write.csv(negative_values, csv_filename, row.names = FALSE)
    cat("Negative values saved to:", csv_filename, "\n")
  } else {
    cat("No records with negative chlorophyll or phaeopigments found.\n")
  }
  cat("\n")
  
  # Check 4: Missing size fractions for each collection time
  cat("=== CHECK 4: MISSING SIZE FRACTIONS ===\n")
  
  if ("collected" %in% names(data) && "size_fraction" %in% names(data)) {
    # Define expected size fractions (excluding Bulk GF/F)
    expected_fractions <- c("GF/F", "3um", "20um")
    
    # Check for missing size fractions by collection time
    missing_fractions <- data %>%
      filter(!is.na(collected), !is.na(size_fraction)) %>%
      select(collected, site_id, line_out_depth, size_fraction, 
             any_of(c("hakai_id", "date"))) %>%
      # Group by collection parameters (excluding Bulk GF/F)
      filter(!grepl("Bulk", size_fraction, ignore.case = TRUE)) %>%
      group_by(collected, site_id, line_out_depth) %>%
      summarise(
        available_fractions = paste(sort(unique(size_fraction)), collapse = ", "),
        n_fractions = n_distinct(size_fraction),
        missing_fractions = paste(setdiff(expected_fractions, unique(size_fraction)), collapse = ", "),
        first_date = min(date, na.rm = TRUE),
        sample_hakai_ids = paste(unique(hakai_id), collapse = ", "),
        .groups = 'drop'
      ) %>%
      filter(n_fractions < 3 | missing_fractions != "") %>%
      arrange(collected, site_id, line_out_depth)
    
    if (nrow(missing_fractions) > 0) {
      cat("Found", nrow(missing_fractions), "collection times missing size fractions:\n")
      print(missing_fractions)
      
      # Create issues directory if it doesn't exist
      if (!dir.exists(here("issues"))) {
        dir.create(here("issues"), recursive = TRUE)
      }
      
      # Save to CSV
      csv_filename <- here("issues", paste0("missing_size_fractions_", Sys.Date(), ".csv"))
      write.csv(missing_fractions, csv_filename, row.names = FALSE)
      cat("Missing size fractions saved to:", csv_filename, "\n")
    } else {
      cat("All collection times have complete size fraction sets (GF/F, 3um, 20um).\n")
    }
  } else {
    missing_cols <- setdiff(c("collected", "size_fraction"), names(data))
    cat("WARNING: Cannot check size fractions. Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    missing_fractions <- data.frame()
  }
  cat("\n")
  
  # Check 5: Replicated hakai_id
  cat("=== CHECK 5: REPLICATED HAKAI_ID ===\n")
  
  if ("hakai_id" %in% names(data)) {
    replicated_ids <- data %>%
      count(hakai_id, name = "n_records") %>%
      filter(n_records > 1) %>%
      left_join(
        data %>%
          select(hakai_id, site_id, line_out_depth, date, 
                 any_of(c("collected", "size_fraction", "chla", "phaeo"))) %>%
          mutate(
            chla = round(chla, 2),
            phaeo = round(phaeo, 2)
          ),
        by = "hakai_id"
      ) %>%
      arrange(hakai_id, date, site_id, line_out_depth)
    
    if (nrow(replicated_ids) > 0) {
      cat("Found", length(unique(replicated_ids$hakai_id)), "replicated hakai_id values in", nrow(replicated_ids), "total records:\n")
      print(replicated_ids)
      
      # Create issues directory if it doesn't exist
      if (!dir.exists(here("issues"))) {
        dir.create(here("issues"), recursive = TRUE)
      }
      
      # Save to CSV
      csv_filename <- here("issues", paste0("replicated_hakai_ids_", Sys.Date(), ".csv"))
      write.csv(replicated_ids, csv_filename, row.names = FALSE)
      cat("Replicated hakai_id records saved to:", csv_filename, "\n")
    } else {
      cat("No replicated hakai_id values found.\n")
    }
  } else {
    cat("WARNING: Cannot check for replicated hakai_id. Missing 'hakai_id' column.\n")
    replicated_ids <- data.frame()
  }
  
  # Check 2: Records with acid ratios < 1.01
  cat("=== CHECK 2: LOW ACID RATIOS ===\n")
  
  # Calculate acid ratio if not already present
  if (!"acid_ratio" %in% names(data)) {
    if (all(c("before_acid", "after_acid") %in% names(data))) {
      data <- data %>%
        mutate(acid_ratio = before_acid / after_acid)
      cat("Calculated acid_ratio from before_acid/after_acid\n")
    } else {
      cat("ERROR: Cannot find acid ratio data. Need either 'acid_ratio' column or 'before_acid' and 'after_acid' columns.\n")
      cat("Available columns:", paste(names(data), collapse = ", "), "\n\n")
      low_acid_ratios <- data.frame()
    }
  }
  
  if ("acid_ratio" %in% names(data) || all(c("before_acid", "after_acid") %in% names(data))) {
    low_acid_ratios <- data %>%
      filter(!is.na(acid_ratio) & acid_ratio < 1.01) %>%
      select(hakai_id, site_id, line_out_depth, date, acid_ratio, 
             any_of(c("before_acid", "after_acid", "chla", "phaeo", "chla_flag", "phaeo_flag"))) %>%
      mutate(
        chla = round(chla, 2),
        phaeo = round(phaeo, 2)
      ) %>%
      arrange(date, site_id, line_out_depth)
    
    if (nrow(low_acid_ratios) > 0) {
      cat("Found", nrow(low_acid_ratios), "records with acid ratios < 1.01:\n")
      print(low_acid_ratios)
      
      # Create issues directory if it doesn't exist
      if (!dir.exists(here("issues"))) {
        dir.create(here("issues"), recursive = TRUE)
      }
      
      # Save to CSV
      csv_filename <- here("issues", paste0("low_acid_ratios_", Sys.Date(), ".csv"))
      write.csv(low_acid_ratios, csv_filename, row.names = FALSE)
      cat("Low acid ratios saved to:", csv_filename, "\n")
    } else {
      cat("No records with acid ratios < 1.01 found.\n")
    }
  }
  cat("\n")
  
  # Check 3: Calibration consistency check
  cat("=== CHECK 3: CALIBRATION CONSISTENCY ===\n")
  
  # Define calibration-related columns to check
  calibration_cols <- c("flurometer_serial_no", "calibration", "acid_ratio_correction_factor", 
                        "acid_coefficient", "calibration_slope")
  
  # Check which calibration columns are available
  available_cal_cols <- calibration_cols[calibration_cols %in% names(data)]
  missing_cal_cols <- calibration_cols[!calibration_cols %in% names(data)]
  
  if (length(missing_cal_cols) > 0) {
    cat("WARNING: Missing calibration columns:", paste(missing_cal_cols, collapse = ", "), "\n")
  }
  
  if (length(available_cal_cols) > 0) {
    cat("Checking calibration consistency for columns:", paste(available_cal_cols, collapse = ", "), "\n\n")
    
    # Check if we have flurometer_serial_no for grouping
    if ("flurometer_serial_no" %in% available_cal_cols) {
      
      # Create calibration summary by fluorometer serial number - include analyzed dates
      calibration_summary <- data %>%
        filter(!is.na(date)) %>%
        select(date, all_of(available_cal_cols), any_of("analyzed")) %>%
        # Remove rows where all calibration values are NA
        filter(if_all(all_of(available_cal_cols), ~ !is.na(.))) %>%
        arrange(date) %>%
        # Group by fluorometer serial number first, then by calibration parameters
        group_by(flurometer_serial_no) %>%
        nest() %>%
        mutate(
          calibration_periods = map(data, ~ {
            summary_data <- .x %>%
              group_by(across(all_of(available_cal_cols[available_cal_cols != "flurometer_serial_no"]))) %>%
              summarise(
                start_date = min(date, na.rm = TRUE),
                end_date = max(date, na.rm = TRUE),
                n_records = n(),
                date_range = paste(min(date, na.rm = TRUE), "to", max(date, na.rm = TRUE)),
                .groups = 'drop'
              )
            
            # Add analyzed date columns if analyzed column exists
            if ("analyzed" %in% names(.x)) {
              summary_data <- .x %>%
                group_by(across(all_of(available_cal_cols[available_cal_cols != "flurometer_serial_no"]))) %>%
                summarise(
                  start_date = min(date, na.rm = TRUE),
                  end_date = max(date, na.rm = TRUE),
                  start_analyzed = min(as.Date(analyzed), na.rm = TRUE),
                  end_analyzed = max(as.Date(analyzed), na.rm = TRUE),
                  n_records = n(),
                  date_range = paste(min(date, na.rm = TRUE), "to", max(date, na.rm = TRUE)),
                  analyzed_range = paste(min(as.Date(analyzed), na.rm = TRUE), "to", max(as.Date(analyzed), na.rm = TRUE)),
                  .groups = 'drop'
                )
            }
            
            summary_data %>% arrange(if ("start_analyzed" %in% names(.)) start_analyzed else start_date)
          })
        ) %>%
        select(-data) %>%
        unnest(calibration_periods)
      
      if (nrow(calibration_summary) > 0) {
        cat("Calibration periods by fluorometer serial number:\n")
        
        # Print summary for each fluorometer
        for (serial_num in unique(calibration_summary$flurometer_serial_no)) {
          cat("\n--- Fluorometer Serial Number:", serial_num, "---\n")
          serial_data <- calibration_summary %>% filter(flurometer_serial_no == serial_num)
          print(serial_data %>% select(-flurometer_serial_no))
        }
        
        # Check for calibration date issues if calibration column exists
        calibration_date_issues <- NULL
        if ("calibration" %in% available_cal_cols && "analyzed" %in% names(data)) {
          cat("\n=== CALIBRATION DATE VALIDATION ===\n")
          
          # Use the already computed calibration_summary for date validation
          calibration_date_issues <- calibration_summary %>%
            mutate(
              calibration_date = as.Date(calibration, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y")),
              date_issue = !is.na(calibration_date) & start_analyzed < calibration_date
            ) %>%
            filter(date_issue == TRUE) %>%
            select(flurometer_serial_no, calibration, calibration_date, start_analyzed, end_analyzed,
                   start_date, end_date, n_records, everything())
          
          if (nrow(calibration_date_issues) > 0) {
            cat("WARNING: Found", nrow(calibration_date_issues), "calibration periods where analyzed date is before calibration date:\n")
            print(calibration_date_issues)
            
            # Create issues directory if it doesn't exist
            if (!dir.exists(here("issues"))) {
              dir.create(here("issues"), recursive = TRUE)
            }
            
            # Save to CSV
            csv_filename <- here("issues", paste0("calibration_date_issues_", Sys.Date(), ".csv"))
            write.csv(calibration_date_issues, csv_filename, row.names = FALSE)
            cat("Calibration date issues saved to:", csv_filename, "\n")
            
          } else {
            cat("No calibration date issues found (analyzed dates are after calibration dates).\n")
          }
          
        } else if ("calibration" %in% available_cal_cols && !"analyzed" %in% names(data)) {
          cat("\n=== CALIBRATION DATE VALIDATION ===\n")
          cat("WARNING: 'analyzed' column not found. Cannot perform calibration date validation.\n")
          cat("Available date columns:", paste(names(data)[grepl("date|time|analyzed", names(data), ignore.case = TRUE)], collapse = ", "), "\n")
        }
        
      } else {
        cat("No complete calibration records found in the data.\n")
      }
      
    } else {
      cat("Cannot perform fluorometer-specific analysis: flurometer_serial_no column not found.\n")
      calibration_summary <- data.frame()
    }
    
  } else {
    cat("No calibration columns found in the data.\n")
    calibration_summary <- data.frame()
  }
  
  # Summary statistics
  cat("\n=== SUMMARY ===\n")
  cat("Total records processed:", nrow(data), "\n")
  cat("Records with negative values:", nrow(negative_values), "\n")
  if (exists("low_acid_ratios")) {
    cat("Records with low acid ratios:", nrow(low_acid_ratios), "\n")
  }
  if (exists("missing_fractions")) {
    cat("Collection times missing size fractions:", nrow(missing_fractions), "\n")
  }
  if (exists("replicated_ids")) {
    cat("Replicated hakai_id values:", length(unique(replicated_ids$hakai_id)), "\n")
  }
  if (exists("calibration_summary")) {
    cat("Distinct calibration combinations:", nrow(calibration_summary), "\n")
  }
  
  # Return results as a list
  results <- list(
    negative_values = negative_values,
    summary_stats = list(
      total_records = nrow(data),
      negative_records = nrow(negative_values)
    )
  )
  
  if (exists("low_acid_ratios")) {
    results$low_acid_ratios <- low_acid_ratios
    results$summary_stats$low_acid_ratio_records <- nrow(low_acid_ratios)
  }
  
  if (exists("missing_fractions")) {
    results$missing_fractions <- missing_fractions
    results$summary_stats$missing_fraction_collections <- nrow(missing_fractions)
  }
  
  if (exists("replicated_ids")) {
    results$replicated_ids <- replicated_ids
    results$summary_stats$replicated_hakai_ids <- length(unique(replicated_ids$hakai_id))
  }
  
  if (exists("calibration_summary")) {
    results$calibration_summary <- calibration_summary
    results$summary_stats$calibration_combinations <- nrow(calibration_summary)
  }
  
  if (exists("calibration_date_issues") && !is.null(calibration_date_issues)) {
    results$calibration_date_issues <- calibration_date_issues
    results$summary_stats$calibration_date_issues <- nrow(calibration_date_issues)
  }
  
  return(results)
}

#---------------------------------------------------------------------------


# Function to create interactive monthly box plots for outlier detection
create_monthly_boxplot_qc <- function(data, selected_depth = NULL, selected_filter_type = NULL, investigation_year = NULL, 
                                      variable = "chla", show_outliers = TRUE, outlier_method = "iqr") {
  
  # Variable configuration
  var_config <- list(
    chla = list(
      column = "chla",
      title = "Chlorophyll",
      units = "mg/m³",
      color = "rgb(34, 197, 94)",
      flag_column = "chla_flag"
    ),
    acid_ratio = list(
      column = "acid_ratio",
      title = "Acid Ratio",
      units = "(before/after)",
      color = "rgb(168, 85, 247)",
      flag_column = "acid_ratio_flag"
    )
  )
  
  if (!variable %in% names(var_config)) {
    stop("Variable must be either 'chla' or 'acid_ratio'")
  }
  
  config <- var_config[[variable]]
  
  # Set defaults if parameters are NULL
  if (is.null(selected_depth)) {
    available_depths <- sort(unique(data$line_out_depth), na.last = TRUE)
    selected_depth <- available_depths[1]
    cat("Using default depth:", selected_depth, "m\n")
  }
  
  if (is.null(selected_filter_type)) {
    available_filters <- unique(data$filter_type)
    if ("GF/F" %in% available_filters) {
      selected_filter_type <- "GF/F"
    } else {
      selected_filter_type <- available_filters[1]
    }
    cat("Using default filter type:", selected_filter_type, "\n")
  }
  
  if (is.null(investigation_year)) {
    available_years <- sort(unique(year(as.Date(data$date))), decreasing = TRUE)
    investigation_year <- available_years[1]
    cat("Using default investigation year:", investigation_year, "\n")
  }
  
  # Calculate acid_ratio if it doesn't exist
  if (variable == "acid_ratio" && !"acid_ratio" %in% names(data)) {
    data <- data %>%
      mutate(acid_ratio = ifelse(after_acid > 0, before_acid / after_acid, NA))
  }
  
  # Filter data for selected depth and filter type
  filtered_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           !is.na(!!sym(config$column)),
           !is.na(date)) %>%
    mutate(
      date = as.Date(date),
      value = as.numeric(!!sym(config$column)),
      year = year(date),
      month = month(date),
      month_name = month.name[month],
      flag = if(config$flag_column %in% names(data)) as.logical(!!sym(config$flag_column)) else FALSE
    ) %>%
    filter(!is.na(value))
  
  # Debug information
  cat("Debug: Filtered data has", nrow(filtered_data), "rows\n")
  cat("Debug: Year range:", min(filtered_data$year, na.rm = TRUE), "to", max(filtered_data$year, na.rm = TRUE), "\n")
  
  # Separate baseline and investigation data
  baseline_data <- filtered_data %>%
    filter(year != investigation_year)
  
  # Remove SVD flagged data from baseline - always use chla_flag for SVD filtering
  if ("chla_flag" %in% names(data)) {
    baseline_data <- baseline_data %>%
      filter(is.na(chla_flag) | chla_flag != "SVD")
    cat("Excluded SVD flagged data from baseline using chla_flag\n")
  }
  
  investigation_data <- filtered_data %>%
    filter(year == investigation_year)
  
  if (nrow(baseline_data) < 10) {
    stop("Insufficient baseline data (< 10 points). Need more historical data.")
  }
  
  # Calculate monthly statistics for baseline data
  monthly_stats <- baseline_data %>%
    group_by(month, month_name) %>%
    summarise(
      median_val = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      iqr = q75 - q25,
      lower_fence = q25 - 1.5 * iqr,
      upper_fence = q75 + 1.5 * iqr,
      whisker_lower = pmax(lower_fence, min(value, na.rm = TRUE)),
      whisker_upper = pmin(upper_fence, max(value, na.rm = TRUE)),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    arrange(month)
  
  # Identify outliers in investigation data
  investigation_with_outliers <- investigation_data %>%
    left_join(monthly_stats, by = c("month", "month_name")) %>%
    mutate(
      is_outlier = case_when(
        outlier_method == "iqr" ~ value < lower_fence | value > upper_fence,
        outlier_method == "quantile" ~ value < quantile(baseline_data$value, 0.05, na.rm = TRUE) | 
          value > quantile(baseline_data$value, 0.95, na.rm = TRUE),
        TRUE ~ FALSE
      ),
      outlier_type = case_when(
        value < lower_fence ~ "Low",
        value > upper_fence ~ "High",
        TRUE ~ "Normal"
      )
    )
  
  # Create the plot
  p <- plot_ly()
  
  # Add baseline box plots for each month
  for (i in 1:nrow(monthly_stats)) {
    month_data <- baseline_data %>%
      filter(month == monthly_stats$month[i])
    
    if (nrow(month_data) > 0) {
      # Fix: Create text vector with same length as month_data
      hover_text <- paste("Month:", monthly_stats$month_name[i],
                          "<br>Baseline Data",
                          "<br>", config$title, ":", round(month_data$value, 3), config$units,
                          "<br>Date:", month_data$date,
                          "<br>Hakai ID:", month_data$hakai_id)
      
      p <- p %>%
        add_boxplot(
          y = month_data$value,
          x = rep(monthly_stats$month_name[i], nrow(month_data)),
          name = paste("Baseline", monthly_stats$month_name[i]),
          marker = list(color = "rgba(59, 130, 246, 0.7)"),
          line = list(color = "rgb(59, 130, 246)"),
          showlegend = i == 1,
          legendgroup = "baseline",
          hoverinfo = "text",
          text = hover_text
        )
    }
  }
  
  # Add investigation year points
  if (nrow(investigation_with_outliers) > 0) {
    # Normal points (not outliers and not flagged)
    normal_points <- investigation_with_outliers %>% 
      filter(!is_outlier & (!flag | is.na(flag)))
    
    if (nrow(normal_points) > 0) {
      p <- p %>%
        add_markers(
          data = normal_points,
          x = ~month_name,
          y = ~value,
          marker = list(color = config$color, size = 10, symbol = "circle"),
          name = paste(investigation_year, "- Normal"),
          hoverinfo = "text",
          text = ~paste("Date:", date, "<br>Hakai ID:", hakai_id,
                        "<br>", config$title, ":", round(value, 3), config$units,
                        "<br>Month:", month_name,
                        "<br>Depth:", selected_depth, "m")
        )
    }
    
    # Outlier points
    outlier_points <- investigation_with_outliers %>% 
      filter(is_outlier)
    
    if (nrow(outlier_points) > 0) {
      p <- p %>%
        add_markers(
          data = outlier_points,
          x = ~month_name,
          y = ~value,
          marker = list(color = "rgb(239, 68, 68)", size = 12, symbol = "diamond"),
          name = paste(investigation_year, "- Outliers"),
          hoverinfo = "text",
          text = ~paste("Date:", date, "<br>Hakai ID:", hakai_id,
                        "<br>", config$title, ":", round(value, 3), config$units,
                        "<br>Month:", month_name,
                        "<br>Outlier Type:", outlier_type,
                        "<br>Depth:", selected_depth, "m",
                        "<br>⚠ OUTLIER")
        )
    }
    
    # Flagged points (already flagged in original data)
    flagged_points <- investigation_with_outliers %>% 
      filter(flag == TRUE)
    
    if (nrow(flagged_points) > 0) {
      p <- p %>%
        add_markers(
          data = flagged_points,
          x = ~month_name,
          y = ~value,
          marker = list(color = "rgb(220, 38, 127)", size = 10, symbol = "x"),
          name = paste(investigation_year, "- Pre-flagged"),
          hoverinfo = "text",
          text = ~paste("Date:", date, "<br>Hakai ID:", hakai_id,
                        "<br>", config$title, ":", round(value, 3), config$units,
                        "<br>Month:", month_name,
                        "<br>Depth:", selected_depth, "m",
                        "<br>⚠ ALREADY FLAGGED")
        )
    }
  }
  
  # Layout
  p <- p %>%
    layout(
      title = list(
        text = paste("Monthly", config$title, "Box Plot QC - Depth:", selected_depth, "m -", selected_filter_type, "Filter<br>", 
                     "<sub>Outlier Method:", toupper(outlier_method), "</sub>"),
        font = list(size = 16, family = "Arial, sans-serif")
      ),
      xaxis = list(
        title = "Month",
        showgrid = TRUE,
        gridcolor = "rgba(0,0,0,0.1)",
        categoryorder = "array",
        categoryarray = month.name
      ),
      yaxis = list(
        title = paste(config$title, "(", config$units, ")"),
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
  
  return(list(plot = p, outlier_data = investigation_with_outliers, monthly_stats = monthly_stats))
}


