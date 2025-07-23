# Streamlined chlorophyll QC function with seasonal GAM
create_chlorophyll_qc_plot <- function(data, selected_depth, selected_filter_type, 
                                       investigation_year, variability_method = "gam_quantiles", 
                                       window_days = 14) {
  
  # Filter and prepare data
  depth_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           !is.na(chla), !is.na(date),
           chla_flag != "SVD" | is.na(chla_flag)) %>%
    mutate(
      date = as.Date(date),
      chla = as.numeric(chla),
      year = year(date),
      day_of_year = yday(date),
      chla_flag = as.logical(chla_flag)
    )
  
  if (nrow(depth_data) < 20) {
    stop("Insufficient data (< 20 points). Need more historical data.")
  }
  
  # Split data
  baseline_data <- depth_data %>%
    filter(year != investigation_year) %>%
    mutate(
      log_chla = log(chla + 0.01),
      cos1 = cos(2 * pi * day_of_year / 365.25),
      sin1 = sin(2 * pi * day_of_year / 365.25),
      cos2 = cos(4 * pi * day_of_year / 365.25),
      sin2 = sin(4 * pi * day_of_year / 365.25)
    )
  
  investigation_data <- depth_data %>% filter(year == investigation_year)
  
  # Prediction grid
  prediction_grid <- tibble(
    day_of_year = 1:365,
    cos1 = cos(2 * pi * day_of_year / 365.25),
    sin1 = sin(2 * pi * day_of_year / 365.25),
    cos2 = cos(4 * pi * day_of_year / 365.25),
    sin2 = sin(4 * pi * day_of_year / 365.25)
  )
  
  # Calculate seasonal patterns and bounds
  if (variability_method == "gam_quantiles") {
    prediction_grid <- fit_gam_model(baseline_data, prediction_grid, window_days)
    method_name <- paste0("GAM + Moving Window (±", window_days, " days)")
  } else {
    prediction_grid <- fit_weekly_model(baseline_data, prediction_grid)
    method_name <- "Weekly Averages"
  }
  
  # Create plot
  create_qc_plotly(prediction_grid, investigation_data, investigation_year, 
                   selected_depth, selected_filter_type, method_name)
}

# GAM fitting helper function
fit_gam_model <- function(baseline_data, prediction_grid, window_days) {
  
  # Fit GAM with fallback to harmonic regression
  tryCatch({
    gam_model <- gam(log_chla ~ s(day_of_year, bs = "cc", k = 20), data = baseline_data)
    prediction_grid$pred_log_median <- predict(gam_model, newdata = prediction_grid)
    baseline_data$residuals <- residuals(gam_model)
    
    # Calculate moving window quantiles
    quantile_results <- map_dfr(1:365, ~calculate_moving_quantiles(.x, baseline_data, window_days))
    
    prediction_grid %>%
      left_join(quantile_results, by = "day_of_year") %>%
      mutate(
        pred_median = exp(pred_log_median) - 0.01,
        lower_bound_90 = pmax(0, exp(pred_log_median + residual_q10) - 0.01),
        upper_bound_90 = pmax(0, exp(pred_log_median + residual_q90) - 0.01),
        lower_bound_95 = pmax(0, exp(pred_log_median + residual_q05) - 0.01),
        upper_bound_95 = pmax(0, exp(pred_log_median + residual_q95) - 0.01)
      )
    
  }, error = function(e) {
    # Fallback to harmonic regression
    harm_model <- lm(log_chla ~ cos1 + sin1 + cos2 + sin2, data = baseline_data)
    pred_log <- predict(harm_model, newdata = prediction_grid)
    
    # Simple quantile adjustment
    q_adjust <- quantile(baseline_data$log_chla, c(0.05, 0.10, 0.90, 0.95), na.rm = TRUE) - 
      mean(baseline_data$log_chla, na.rm = TRUE)
    
    prediction_grid %>%
      mutate(
        pred_median = exp(pred_log) - 0.01,
        lower_bound_90 = pmax(0, exp(pred_log + q_adjust[2]) - 0.01),
        upper_bound_90 = pmax(0, exp(pred_log + q_adjust[3]) - 0.01),
        lower_bound_95 = pmax(0, exp(pred_log + q_adjust[1]) - 0.01),
        upper_bound_95 = pmax(0, exp(pred_log + q_adjust[4]) - 0.01)
      )
  })
}

# Weekly model helper function
fit_weekly_model <- function(baseline_data, prediction_grid) {
  weekly_stats <- baseline_data %>%
    mutate(week = week(as.Date(paste("2020", day_of_year), format = "%Y %j"))) %>%
    group_by(week) %>%
    summarise(
      pred_median = median(exp(log_chla) - 0.01, na.rm = TRUE),
      lower_bound_90 = quantile(exp(log_chla) - 0.01, 0.10, na.rm = TRUE),
      upper_bound_90 = quantile(exp(log_chla) - 0.01, 0.90, na.rm = TRUE),
      lower_bound_95 = quantile(exp(log_chla) - 0.01, 0.05, na.rm = TRUE),
      upper_bound_95 = quantile(exp(log_chla) - 0.01, 0.95, na.rm = TRUE),
      .groups = 'drop'
    )
  
  prediction_grid %>%
    mutate(week = week(as.Date(paste("2020", day_of_year), format = "%Y %j"))) %>%
    left_join(weekly_stats, by = "week")
}

# Moving window quantile calculation
calculate_moving_quantiles <- function(target_day, baseline_data, window_days) {
  # Handle year wrap-around
  day_range <- if (target_day <= window_days) {
    c((365 + target_day - window_days):365, 1:(target_day + window_days))
  } else if (target_day > (365 - window_days)) {
    c((target_day - window_days):365, 1:(target_day + window_days - 365))
  } else {
    (target_day - window_days):(target_day + window_days)
  }
  
  window_data <- baseline_data %>% filter(day_of_year %in% day_range)
  
  if (nrow(window_data) < 5) {
    return(tibble(day_of_year = target_day, 
                  residual_q05 = NA, residual_q10 = NA, 
                  residual_q90 = NA, residual_q95 = NA))
  }
  
  quantiles <- quantile(window_data$residuals, c(0.05, 0.10, 0.90, 0.95), na.rm = TRUE)
  tibble(
    day_of_year = target_day,
    residual_q05 = quantiles[1],
    residual_q10 = quantiles[2],
    residual_q90 = quantiles[3],
    residual_q95 = quantiles[4]
  )
}

# Plotly creation helper
create_qc_plotly <- function(prediction_grid, investigation_data, investigation_year, 
                             selected_depth, selected_filter_type, method_name) {
  
  p <- plot_ly() %>%
    # 95% confidence ribbon
    add_ribbons(
      data = prediction_grid, x = ~day_of_year,
      ymin = ~lower_bound_95, ymax = ~upper_bound_95,
      fillcolor = "rgba(59, 130, 246, 0.1)", line = list(color = "transparent"),
      name = "95% Range", hovertemplate = "Day: %{x}<br>95% Range: %{ymin:.3f}-%{ymax:.3f} mg/m³"
    ) %>%
    # 90% confidence ribbon
    add_ribbons(
      data = prediction_grid, x = ~day_of_year,
      ymin = ~lower_bound_90, ymax = ~upper_bound_90,
      fillcolor = "rgba(59, 130, 246, 0.2)", line = list(color = "transparent"),
      name = "90% Range", hovertemplate = "Day: %{x}<br>90% Range: %{ymin:.3f}-%{ymax:.3f} mg/m³"
    ) %>%
    # Median line
    add_lines(
      data = prediction_grid, x = ~day_of_year, y = ~pred_median,
      line = list(color = "rgb(59, 130, 246)", width = 3),
      name = "Seasonal Median", hovertemplate = "Day: %{x}<br>Median: %{y:.3f} mg/m³"
    )
  
  # Add investigation year data
  if (nrow(investigation_data) > 0) {
    # Normal points
    normal_points <- investigation_data %>% filter(!chla_flag | is.na(chla_flag))
    if (nrow(normal_points) > 0) {
      p <- p %>% add_markers(
        data = normal_points, x = ~day_of_year, y = ~chla,
        marker = list(color = "rgb(34, 197, 94)", size = 8),
        name = paste(investigation_year, "- Normal"),
        hovertemplate = "Date: %{customdata}<br>Chlorophyll: %{y:.3f} mg/m³<br>Depth: %{text} m",
        customdata = ~date, text = ~line_out_depth
      )
    }
    
    # Flagged points
    flagged_points <- investigation_data %>% filter(chla_flag == TRUE)
    if (nrow(flagged_points) > 0) {
      p <- p %>% add_markers(
        data = flagged_points, x = ~day_of_year, y = ~chla,
        marker = list(color = "rgb(239, 68, 68)", size = 8),
        name = paste(investigation_year, "- Flagged"),
        hovertemplate = "Date: %{customdata}<br>Chlorophyll: %{y:.3f} mg/m³<br>Depth: %{text} m<br>⚠ FLAGGED",
        customdata = ~date, text = ~line_out_depth
      )
    }
  }
  
  # Layout
  p %>% layout(
    title = paste("Chlorophyll QC - Depth:", selected_depth, "m -", selected_filter_type, 
                  "Filter<br><sub>Method:", method_name, "</sub>"),
    xaxis = list(title = "Day of Year", showgrid = TRUE, range = c(1, 365)),
    yaxis = list(title = "Chlorophyll (mg/m³)", showgrid = TRUE),
    hovermode = "closest",
    legend = list(x = 0.01, y = 0.99, bgcolor = "rgba(255,255,255,0.8)"),
    plot_bgcolor = "white"
  )
}

# Streamlined summary function
generate_summary_stats <- function(data, selected_depth, selected_filter_type, investigation_year) {
  investigation_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           year(as.Date(date)) == investigation_year,
           !is.na(chla),
           chla_flag != "SVD" | is.na(chla_flag)) %>%
    mutate(chla = as.numeric(chla), chla_flag = as.logical(chla_flag))
  
  baseline_data <- data %>%
    filter(line_out_depth == selected_depth,
           filter_type == selected_filter_type,
           year(as.Date(date)) != investigation_year,
           !is.na(chla),
           chla_flag != "SVD" | is.na(chla_flag))
  
  if (nrow(investigation_data) == 0) {
    return(list(total_samples = 0, flagged_samples = 0, mean_chla = NA, 
                date_range = "No data", baseline_years = 0, baseline_samples = 0))
  }
  
  list(
    total_samples = nrow(investigation_data),
    flagged_samples = sum(investigation_data$chla_flag, na.rm = TRUE),
    mean_chla = mean(investigation_data$chla, na.rm = TRUE),
    median_chla = median(investigation_data$chla, na.rm = TRUE),
    date_range = paste(range(investigation_data$date), collapse = " - "),
    baseline_years = length(unique(year(as.Date(baseline_data$date)))),
    baseline_samples = nrow(baseline_data)
  )
}

# Main function with cleaner output
timeseries_qc <- function(data, selected_depth = NULL, selected_filter_type = NULL, 
                          investigation_year = NULL, variability_method = "gam_quantiles",
                          window_days = 14) {
  
  # Get defaults if not provided
  options <- get_available_options(data)
  selected_depth <- selected_depth %||% options$depths[1]
  selected_filter_type <- selected_filter_type %||% options$filter_types[1]
  investigation_year <- investigation_year %||% options$years[1]
  
  # Create plot and stats
  plot <- create_chlorophyll_qc_plot(data, selected_depth, selected_filter_type, 
                                     investigation_year, variability_method, window_days)
  stats <- generate_summary_stats(data, selected_depth, selected_filter_type, investigation_year)
  
  # Print concise summary
  cat("Time Series QC Summary:\n",
      "Depth: ", selected_depth, "m | Filter: ", selected_filter_type, 
      " | Year: ", investigation_year, "\n",
      "Samples: ", stats$total_samples, " (", stats$flagged_samples, " flagged) | ",
      "Mean: ", round(stats$mean_chla, 3), " mg/m³\n",
      "Baseline: ", stats$baseline_samples, " samples from ", stats$baseline_years, " years\n",
      sep = "")
  
  print(plot)
  list(plot = plot, stats = stats, options = options)
}

# Helper function for available options
get_available_options <- function(data) {
  list(
    depths = data %>% filter(!is.na(line_out_depth)) %>% 
      distinct(line_out_depth) %>% arrange(as.numeric(line_out_depth)) %>% pull(),
    years = data %>% filter(!is.na(date)) %>% 
      mutate(year = year(as.Date(date))) %>% distinct(year) %>% 
      arrange(desc(year)) %>% pull(),
    filter_types = data %>% filter(!is.na(filter_type)) %>% 
      distinct(filter_type) %>% pull()
  )
}

# Null coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x

# ----------------------------------------------------------------------------

# Enhanced function to sum size-fractions and compare to bulk with acid ratio QC
sum_size_fractions <- function(df) {
  
  # Check if required columns exist
  required_cols <- c("date", "site_id", "line_out_depth", "hakai_id", 
                     "collected", "filter_type", "chla", "chla_flag",
                     "before_acid", "after_acid")
  missing_cols <- setdiff(required_cols, names(df))
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Create issues folder if it doesn't exist
  if(!dir.exists("issues")) {
    dir.create("issues")
  }
  
  # Define size classes, bulk, and bad flags
  size_classes <- c("20um", "3um", "GF/F")
  bulk_class <- "Bulk GF/F"
  bad_flags <- c("SVD")
  
  # Calculate acid ratios for all data
  df <- df %>%
    mutate(
      acid_ratio = before_acid / after_acid,
      acid_ratio_flag = case_when(
        is.na(acid_ratio) ~ "Missing",
        acid_ratio < 1.1 ~ "Low",
        acid_ratio > 2.0 ~ "High",
        TRUE ~ "OK"
      )
    )
  
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
      # Check for any problematic acid ratios in size fractions
      any_bad_acid_ratio = any(acid_ratio_flag %in% c("Low", "High")),
      bad_acid_details = paste(filter_type[acid_ratio_flag %in% c("Low", "High")], 
                               collapse = ", "),
      .groups = "drop"
    )
  
  # Get bulk measurements (excluding bad quality) with acid ratio info
  bulk_data <- df %>%
    filter(filter_type == bulk_class) %>%
    filter(!chla_flag %in% bad_flags) %>%  # Remove bad quality bulk measurements
    select(collected, date, site_id, line_out_depth, 
           chla_bulk = chla, chla_bulk_flag = chla_flag, hakai_id,
           bulk_acid_ratio = acid_ratio, bulk_acid_ratio_flag = acid_ratio_flag)
  
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
                             ((chla_sum - chla_bulk) / chla_bulk) * 100, NA_real_),
      # Combined acid ratio flag (bad if ANY filter type has bad acid ratio)
      combined_acid_flag = case_when(
        any_bad_acid_ratio | bulk_acid_ratio_flag %in% c("Low", "High") ~ "Bad",
        TRUE ~ "OK"
      )
    )
  
  return(comparison_data)
}

# Function to detect outliers using residuals
detect_outliers <- function(comparison_data, residual_threshold = 2.5) {
  
  # Filter to valid comparisons
  valid_data <- comparison_data %>%
    filter(!is.na(chla_diff) & complete_set)
  
  if(nrow(valid_data) < 10) {
    warning("Insufficient data for outlier detection (need at least 10 points)")
    return(comparison_data %>% mutate(outlier_flag = "Insufficient_data",
                                      high_confidence_outlier = FALSE))
  }
  
  # Fit regression model
  model <- lm(chla_sum ~ chla_bulk, data = valid_data)
  
  # Calculate residuals and standardized residuals
  valid_data <- valid_data %>%
    mutate(
      predicted = predict(model, newdata = .),
      residual = chla_sum - predicted,
      standardized_residual = residual / sd(residual, na.rm = TRUE)
    )
  
  # Flag outliers based on standardized residuals
  valid_data <- valid_data %>%
    mutate(
      outlier_flag = case_when(
        abs(standardized_residual) > residual_threshold ~ "Outlier",
        TRUE ~ "Normal"
      ),
      # High confidence outlier: outlier AND bad acid ratio
      high_confidence_outlier = outlier_flag == "Outlier" & combined_acid_flag == "Bad"
    )
  
  # Join back with original data
  result <- comparison_data %>%
    left_join(
      valid_data %>% select(collected, date, site_id, line_out_depth, 
                            predicted, residual, standardized_residual, 
                            outlier_flag, high_confidence_outlier),
      by = c("collected", "date", "site_id", "line_out_depth")
    ) %>%
    mutate(
      outlier_flag = if_else(is.na(outlier_flag), "No_comparison", outlier_flag),
      high_confidence_outlier = if_else(is.na(high_confidence_outlier), FALSE, high_confidence_outlier)
    )
  
  return(result)
}

# Function to export problematic samples to issues folder
export_problem_samples <- function(comparison_data, original_df) {
  
  # Create timestamp for filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Function to get size fraction hakai_ids for problematic samples
  get_size_fraction_ids <- function(problem_data) {
    if(nrow(problem_data) == 0) return(problem_data)
    
    # Get the size fraction hakai_ids for each problematic sample
    size_fraction_ids <- problem_data %>%
      left_join(
        original_df %>%
          filter(filter_type %in% c("20um", "3um", "GF/F")) %>%
          group_by(collected, date, site_id, line_out_depth) %>%
          summarise(
            size_fraction_hakai_ids = paste(hakai_id, collapse = "; "),
            .groups = "drop"
          ),
        by = c("collected", "date", "site_id", "line_out_depth")
      ) %>%
      # Keep both bulk and size fraction IDs for reference
      rename(bulk_hakai_id = hakai_id) %>%
      select(collected, date, site_id, line_out_depth, year, 
             bulk_hakai_id, size_fraction_hakai_ids, everything())
    
    return(size_fraction_ids)
  }
  
  # High confidence outliers
  high_conf_outliers <- comparison_data %>%
    filter(high_confidence_outlier == TRUE) %>%
    arrange(date, site_id, line_out_depth) %>%
    get_size_fraction_ids()
  
  if(nrow(high_conf_outliers) > 0) {
    write.csv(high_conf_outliers, 
              file = paste0("issues/sf_qc_high_confidence_outliers_", timestamp, ".csv"),
              row.names = FALSE)
    cat("Exported", nrow(high_conf_outliers), "high confidence outliers to issues folder\n")
  }
  
  # All outliers (for reference)
  all_outliers <- comparison_data %>%
    filter(outlier_flag == "Outlier") %>%
    arrange(date, site_id, line_out_depth) %>%
    get_size_fraction_ids()
  
  if(nrow(all_outliers) > 0) {
    write.csv(all_outliers, 
              file = paste0("issues/sf_qc_all_outliers_", timestamp, ".csv"),
              row.names = FALSE)
    cat("Exported", nrow(all_outliers), "total outliers to issues folder\n")
  }
  
  # Samples with bad acid ratios (regardless of outlier status)
  bad_acid_samples <- comparison_data %>%
    filter(combined_acid_flag == "Bad") %>%
    arrange(date, site_id, line_out_depth) %>%
    get_size_fraction_ids()
  
  if(nrow(bad_acid_samples) > 0) {
    write.csv(bad_acid_samples, 
              file = paste0("issues/sf_qc_bad_acid_ratios_", timestamp, ".csv"),
              row.names = FALSE)
    cat("Exported", nrow(bad_acid_samples), "samples with bad acid ratios to issues folder\n")
  }
  
  # Summary report
  summary_report <- data.frame(
    Category = c("Total samples", "Valid comparisons", "Outliers", 
                 "Bad acid ratios", "High confidence outliers"),
    Count = c(nrow(comparison_data),
              nrow(comparison_data %>% filter(outlier_flag %in% c("Outlier", "Normal"))),
              nrow(all_outliers),
              nrow(bad_acid_samples),
              nrow(high_conf_outliers)),
    Percentage = c(100,
                   round(nrow(comparison_data %>% filter(outlier_flag %in% c("Outlier", "Normal"))) / nrow(comparison_data) * 100, 1),
                   round(nrow(all_outliers) / nrow(comparison_data) * 100, 1),
                   round(nrow(bad_acid_samples) / nrow(comparison_data) * 100, 1),
                   round(nrow(high_conf_outliers) / nrow(comparison_data) * 100, 1))
  )
  
  write.csv(summary_report, 
            file = paste0("issues/sf_qc_summary_", timestamp, ".csv"),
            row.names = FALSE)
  
  return(summary_report)
}

# Enhanced interactive plotting function with outlier highlighting
plot_bulk_vs_sum_interactive <- function(comparison_data, selected_year = NULL, highlight_outliers = TRUE) {
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
  
  # Create point colors and shapes based on year AND outlier status
  if(!is.null(selected_year)) {
    # When year is selected, prioritize year grouping with outlier overlay
    plot_data <- plot_data %>%
      mutate(
        point_color = year_group,
        point_shape = case_when(
          high_confidence_outlier == TRUE ~ "High Confidence Outlier",
          outlier_flag == "Outlier" ~ "Outlier", 
          combined_acid_flag == "Bad" ~ "Bad Acid Ratio",
          TRUE ~ "Normal"
        ),
        point_size = case_when(
          high_confidence_outlier == TRUE ~ 4,
          outlier_flag == "Outlier" ~ 3,
          combined_acid_flag == "Bad" ~ 2.5,
          TRUE ~ 2
        )
      )
  } else {
    # When no year selected, use outlier status for coloring
    plot_data <- plot_data %>%
      mutate(
        point_color = case_when(
          high_confidence_outlier == TRUE ~ "High Confidence Outlier",
          outlier_flag == "Outlier" ~ "Outlier",
          combined_acid_flag == "Bad" ~ "Bad Acid Ratio",
          TRUE ~ "Normal"
        ),
        point_shape = "Normal",
        point_size = 2
      )
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
      title = "Interactive Bulk vs Sum QC Plot with Outlier Detection"
    ) +
    theme_minimal()
  
  # Add year-specific regression line if available
  if(!is.null(year_lm)) {
    year_data <- plot_data %>% filter(year == selected_year)
    p <- p + geom_smooth(data = year_data, method = "lm", se = FALSE, 
                         color = "darkred", linetype = "solid", size = 1.2)
  }
  
  # Add points with hover information
  hover_text <- paste("Hakai ID:", plot_data$hakai_id,
                      "<br>Date:", plot_data$date,
                      "<br>Site:", plot_data$site_id,
                      "<br>Depth:", plot_data$line_out_depth, "m",
                      "<br>Year:", plot_data$year,
                      "<br>Bulk:", round(plot_data$chla_bulk, 3),
                      "<br>Sum:", round(plot_data$chla_sum, 3),
                      "<br>% Diff:", round(plot_data$percent_diff, 1),
                      "<br>Outlier:", plot_data$outlier_flag,
                      "<br>Acid Flag:", plot_data$combined_acid_flag)
  
  # Add points with appropriate styling
  if(!is.null(selected_year)) {
    # Create named vector for colors dynamically
    year_label <- paste("Year", selected_year)
    color_values <- c("grey50", "darkred")
    names(color_values) <- c("Other years", year_label)
    
    p <- p + geom_point(aes(color = point_color, 
                            shape = point_shape, 
                            size = point_size,
                            text = hover_text), 
                        alpha = 0.8) +
      scale_color_manual(values = color_values) +
      scale_shape_manual(values = c("Normal" = 16, 
                                    "Bad Acid Ratio" = 17,
                                    "Outlier" = 4, 
                                    "High Confidence Outlier" = 8)) +
      scale_size_identity() +
      labs(color = "Year Group", shape = "QC Status")
  } else {
    # Outlier-based coloring
    p <- p + geom_point(aes(color = point_color, text = hover_text), 
                        alpha = 0.8, size = 2) +
      scale_color_manual(values = c("Normal" = "steelblue", 
                                    "Bad Acid Ratio" = "orange",
                                    "Outlier" = "red", 
                                    "High Confidence Outlier" = "darkred")) +
      labs(color = "Sample Status")
  }
  
  # Add regression info to subtitle
  r2_full <- summary(full_lm)$r.squared
  slope_full <- coef(full_lm)[2]
  
  outlier_count <- sum(plot_data$outlier_flag == "Outlier", na.rm = TRUE)
  high_conf_count <- sum(plot_data$high_confidence_outlier == TRUE, na.rm = TRUE)
  
  subtitle_text <- paste0("Full dataset: R² = ", round(r2_full, 3), 
                          ", Slope = ", round(slope_full, 3))
  
  if(!is.null(year_lm)) {
    r2_year <- summary(year_lm)$r.squared
    slope_year <- coef(year_lm)[2]
    subtitle_text <- paste0(subtitle_text, " | Year ", selected_year, 
                            ": R² = ", round(r2_year, 3), 
                            ", Slope = ", round(slope_year, 3))
  }
  
  subtitle_text <- paste0(subtitle_text, " | Outliers: ", outlier_count,
                          " | High Confidence: ", high_conf_count)
  
  p <- p + labs(subtitle = subtitle_text)
  
  # Convert to interactive plot
  ggplotly(p, tooltip = "text", height = 600, width = 800)
}

# Complete workflow function
run_complete_qc <- function(df, residual_threshold = 2.5, export_issues = TRUE) {
  
  cat("=== RUNNING COMPLETE QC WORKFLOW ===\n\n")
  
  # Step 1: Calculate size fraction sums and basic QC
  cat("Step 1: Calculating size fraction sums and acid ratios...\n")
  comparison_data <- sum_size_fractions(df)
  
  # Step 2: Detect outliers
  cat("Step 2: Detecting outliers...\n")
  comparison_data <- detect_outliers(comparison_data, residual_threshold)
  
  # Step 3: Export problem samples
  if(export_issues) {
    cat("Step 3: Exporting problem samples...\n")
    summary_report <- export_problem_samples(comparison_data, df)
    print(summary_report)
  }
  
  # Step 4: Print summary
  cat("\n=== QC SUMMARY ===\n")
  cat("Total samples processed:", nrow(comparison_data), "\n")
  cat("Valid comparisons:", sum(!is.na(comparison_data$chla_diff) & comparison_data$complete_set), "\n")
  cat("Outliers detected:", sum(comparison_data$outlier_flag == "Outlier", na.rm = TRUE), "\n")
  cat("High confidence outliers:", sum(comparison_data$high_confidence_outlier == TRUE, na.rm = TRUE), "\n")
  cat("Samples with bad acid ratios:", sum(comparison_data$combined_acid_flag == "Bad", na.rm = TRUE), "\n")
  
  return(comparison_data)
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

# ----------------------------------------------------------------------------

#Code for HPLC comparisons - not sure if this is useful.

# HPLC Chlorophyll Comparison Function for QC Workflow
# Author: Biological Oceanography QC Tool
# Purpose: Compare Fluorophotometric vs HPLC chlorophyll measurements

#' Compare chlorophyll-a measurements between Fluorophotometry and HPLC
#' 
#' @param data Dataframe with Fluorophotometric chlorophyll data
#' @param data_hplc Dataframe with HPLC chlorophyll data  
#' @param selected_year Numeric year to highlight (optional)
#' @param interactive Logical, create interactive plot (default TRUE)
#' @param show_annual_stats Logical, display annual regression statistics (default TRUE)
#' @return ggplot or plotly object depending on interactive parameter

hplc_compare_plot <- function(data, data_hplc, selected_year = NULL, 
                              interactive = TRUE, show_annual_stats = TRUE) {
  
  # Load required libraries
  require(dplyr)
  require(ggplot2)
  require(plotly)
  require(lubridate)
  
  # Data processing
  message("Processing data...")
  
  # Filter and join datasets
  data_filtered <- data %>%
    filter(filter_type == "Bulk GF/F")
  
  combined_data <- data_filtered %>%
    inner_join(data_hplc, 
               by = c("collected", "site_id", "line_out_depth"),
               suffix = c("_chla", "_hplc")) %>%
    # Remove samples with SVD flags (but keep NA flags)
    filter((is.na(chla_flag) | chla_flag != "SVD") & 
             (is.na(all_chl_a_flag) | all_chl_a_flag != "SVD")) %>%
    # Extract year from date
    mutate(
      year_chla = year(date_chla),
      year_hplc = year(date_hplc),
      year = year_chla
    ) %>%
    # Remove rows with missing chlorophyll values
    filter(!is.na(chla) & !is.na(all_chl_a))
  
  if(nrow(combined_data) == 0) {
    stop("No data remaining after filtering. Check your datasets and filter criteria.")
  }
  
  # Calculate baseline regression (all data)
  baseline_model <- lm(all_chl_a ~ chla, data = combined_data)
  baseline_summary <- summary(baseline_model)
  
  # Create prediction intervals for baseline
  chla_range <- seq(min(combined_data$chla, na.rm = TRUE), 
                    max(combined_data$chla, na.rm = TRUE), 
                    length.out = 100)
  baseline_pred <- predict(baseline_model, 
                           newdata = data.frame(chla = chla_range), 
                           interval = "confidence", 
                           level = 0.95)
  
  baseline_bands <- data.frame(
    chla = chla_range,
    fit = baseline_pred[, "fit"],
    lwr = baseline_pred[, "lwr"],
    upr = baseline_pred[, "upr"]
  )
  
  # Annual regression if year is selected
  annual_model <- NULL
  annual_summary <- NULL
  annual_bands <- NULL
  year_data <- NULL
  
  if (!is.null(selected_year)) {
    year_data <- combined_data %>% filter(year == selected_year)
    
    if (nrow(year_data) >= 3) {  # Need at least 3 points for regression
      annual_model <- lm(all_chl_a ~ chla, data = year_data)
      annual_summary <- summary(annual_model)
      
      # Create prediction intervals for annual regression
      if (nrow(year_data) > 0) {
        annual_pred <- predict(annual_model, 
                               newdata = data.frame(chla = chla_range), 
                               interval = "confidence", 
                               level = 0.95)
        annual_bands <- data.frame(
          chla = chla_range,
          fit = annual_pred[, "fit"],
          lwr = annual_pred[, "lwr"],
          upr = annual_pred[, "upr"]
        )
      }
    }
  }
  
  # Determine axis limits
  max_val <- max(c(combined_data$chla, combined_data$all_chl_a), na.rm = TRUE)
  min_val <- min(c(combined_data$chla, combined_data$all_chl_a), na.rm = TRUE)
  
  # Create plot
  p <- ggplot() +
    # Baseline confidence bands
    geom_ribbon(data = baseline_bands, 
                aes(x = chla, ymin = lwr, ymax = upr), 
                alpha = 0.2, fill = "blue") +
    
    # Annual confidence bands (if available)
    {if (!is.null(annual_bands)) {
      geom_ribbon(data = annual_bands, 
                  aes(x = chla, ymin = lwr, ymax = upr), 
                  alpha = 0.3, fill = "darkgreen")
    }} +
    
    # 1:1 line
    geom_abline(intercept = 0, slope = 1, 
                color = "red", linetype = "dashed", size = 1) +
    
    # Baseline regression line
    geom_line(data = baseline_bands, 
              aes(x = chla, y = fit), 
              color = "blue", size = 1.2) +
    
    # Annual regression line (if available)
    {if (!is.null(annual_bands)) {
      geom_line(data = annual_bands, 
                aes(x = chla, y = fit), 
                color = "darkgreen", size = 1.2)
    }} +
    
    # All data points (baseline)
    geom_point(data = combined_data, 
               aes(x = chla, y = all_chl_a, 
                   text = paste("Date (chla):", date_chla, 
                                "<br>Date (HPLC):", date_hplc,
                                "<br>Site:", site_id,
                                "<br>Depth:", line_out_depth,
                                "<br>Chl-a:", round(chla, 3),
                                "<br>All Chl-a:", round(all_chl_a, 3),
                                "<br>Year:", year)), 
               alpha = 0.4, color = "gray60", size = 1.5) +
    
    # Selected year points (if specified)
    {if (!is.null(year_data) && nrow(year_data) > 0) {
      geom_point(data = year_data, 
                 aes(x = chla, y = all_chl_a,
                     text = paste("Date (chla):", date_chla, 
                                  "<br>Date (HPLC):", date_hplc,
                                  "<br>Site:", site_id,
                                  "<br>Depth:", line_out_depth,
                                  "<br>Chl-a:", round(chla, 3),
                                  "<br>All Chl-a:", round(all_chl_a, 3),
                                  "<br>Year:", year)), 
                 color = "darkgreen", size = 2.5, alpha = 0.9)
    }} +
    
    # Labels and formatting
    labs(
      x = "Chlorophyll-a (μg/L) - Fluorometric",
      y = "All Chlorophyll-a (μg/L) - HPLC",
      title = ifelse(!is.null(selected_year), 
                     paste("Chlorophyll Comparison: Year", selected_year, "vs Baseline"),
                     "Chlorophyll Comparison: All Years"),
      caption = paste("Red dashed: 1:1 line | Blue: Baseline regression ± 95% CI",
                      ifelse(!is.null(selected_year), "| Green: Annual regression ± 95% CI", ""))
    ) +
    
    coord_fixed(ratio = 1, xlim = c(min_val, max_val * 1.05), 
                ylim = c(min_val, max_val * 1.05)) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
  
  # Print statistics
  available_years <- sort(unique(combined_data$year))
  
  cat("\n=== HPLC vs Fluorophotometric Chlorophyll Comparison ===\n")
  cat("Total samples after filtering:", nrow(combined_data), "\n")
  cat("Available years:", paste(available_years, collapse = ", "), "\n\n")
  
  cat("BASELINE REGRESSION (All Data):\n")
  cat("Equation: HPLC = ", round(coef(baseline_model)[1], 4), " + ", 
      round(coef(baseline_model)[2], 4), " × Fluoro\n", sep="")
  cat("R²:", round(baseline_summary$r.squared, 4), "\n")
  cat("RMSE:", round(sqrt(mean(baseline_model$residuals^2)), 4), "\n")
  cat("Correlation:", round(cor(combined_data$chla, combined_data$all_chl_a), 4), "\n")
  cat("p-value:", format.pval(baseline_summary$coefficients[2,4]), "\n\n")
  
  if (!is.null(selected_year) && !is.null(annual_model) && show_annual_stats) {
    cat("ANNUAL REGRESSION (", selected_year, "):\n", sep="")
    cat("Samples:", nrow(year_data), "\n")
    cat("Equation: HPLC = ", round(coef(annual_model)[1], 4), " + ", 
        round(coef(annual_model)[2], 4), " × Fluoro\n", sep="")
    cat("R²:", round(annual_summary$r.squared, 4), "\n")
    cat("RMSE:", round(sqrt(mean(annual_model$residuals^2)), 4), "\n")
    cat("Correlation:", round(cor(year_data$chla, year_data$all_chl_a), 4), "\n")
    cat("p-value:", format.pval(annual_summary$coefficients[2,4]), "\n")
    
    # Compare slopes
    slope_diff <- coef(annual_model)[2] - coef(baseline_model)[2]
    cat("Slope difference from baseline:", round(slope_diff, 4), "\n\n")
  }
  
  # Return plot
  if (interactive) {
    return(ggplotly(p, tooltip = "text") %>%
             layout(title = list(text = p$labels$title, font = list(size = 16))))
  } else {
    return(p)
  }
}

# Convenience function for quick year comparisons
quick_year_compare <- function(data, data_hplc, year) {
  return(hplc_compare_plot(data, data_hplc, selected_year = year, 
                           interactive = TRUE, show_annual_stats = TRUE))
}
