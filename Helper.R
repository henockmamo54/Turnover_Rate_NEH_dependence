library(dplyr)
library(ggplot2)
library(tidyr) # For data manipulation functions like unnest() and nest()

#' Adjust Px(t) values based on NEH dependence.
#'
#' This function performs a linear regression of Px(t) against the Number of
#' Exchangeable Hydrogens (NEH) for each time point. It then calculates an
#' adjusted Px(t) value by subtracting the NEH-dependent component.
#'
#' @param pxt_df A data frame containing at least `Pxt`, `NEH`, and `Time` columns.
#' @return A data frame with the same columns as `pxt_df`, plus `slope`,
#'   `intercept`, and `Pxt_adj` (the adjusted Px(t) value).
#' @importFrom dplyr group_by reframe left_join mutate
#' @importFrom stats lm coef
adjust_pxt <- function(pxt_df) {
  # For each time point, fit a linear model Pxt ~ NEH and extract parameters.
  # `reframe` is used as the modern replacement for `do`.
  model_params <- pxt_df |>
    group_by(Time) |>
    reframe({
      fit <- lm(Pxt ~ NEH, data = cur_data())
      tibble::tibble(
        slope     = round(coef(fit)[2], 4),
        intercept = round(coef(fit)[1], 4)
      )
    })

  # Join model parameters back to the original data and compute adjusted Px(t).
  pxt_data_joined <- pxt_df |>
    left_join(model_params, by = "Time") |>
    mutate(Pxt_adj = Pxt - NEH * slope)

  return(pxt_data_joined)
}


#' Adjust Px(t) values for each protein separately.
#'
#' This function is similar to `adjust_pxt`, but it performs the linear
#' regression and adjustment for each protein (identified by `SourceFile`)
#' at each time point.
#'
#' @param pxt_df A data frame containing `Pxt`, `NEH`, `Time`, and `SourceFile`.
#' @return A data frame with added `slope`, `intercept`, and `Pxt_adj` columns.
#' @importFrom dplyr group_by reframe left_join mutate
#' @importFrom stats lm coef
adjust_pxt_all <- function(pxt_df) {
  # Group by both Time and SourceFile to get protein-specific model parameters.
  model_params <- pxt_df |>
    group_by(Time, SourceFile) |>
    reframe({
      fit <- lm(Pxt ~ NEH, data = cur_data())
      tibble::tibble(
        slope     = round(coef(fit)[2], 4),
        intercept = round(coef(fit)[1], 4)
      )
    })

  # Join parameters and calculate the adjusted Px(t).
  pxt_data_joined <- pxt_df |>
    left_join(model_params, by = c("Time", "SourceFile")) |>
    mutate(Pxt_adj = Pxt - NEH * slope)

  return(pxt_data_joined)
}

 
make_pxt_facet <- function(df, y_var, y_lab, tag,
                           out_dir = ".", w_mm = 180, h_mm = 120,
                           units = "mm") {
  time_labels <- c(
    `0` ="No labelling", `1` ="1-day",  `2` ="2-day",  `3` ="3-day",
    `4` ="4-day",        `5` ="5-day",  `6` ="6-day",
    `14`="13-day",       `21`="21-day"
  )
  df <- dplyr::mutate(df, Time = factor(Time, levels = names(time_labels)))

  # Calculate regression statistics for each time point for plot annotations.
  reg_stats <- df |>
    dplyr::group_by(Time) |>
    dplyr::reframe({
      # Construct the formula dynamically to correctly use the y_var column per group
      fit <- lm(as.formula(paste(y_var, "~ NEH")), data = cur_data())
      summary_fit <- summary(fit)
      tibble::tibble(
        slope     = coef(fit)[2],
        intercept = coef(fit)[1],
        ci_l      = confint(fit, level = 0.95)[2, 1],
        ci_u      = confint(fit, level = 0.95)[2, 2],
        mean_y    = mean(.data[[y_var]], na.rm = TRUE),
        sd_y      = sd(.data[[y_var]], na.rm = TRUE),
        pval      = summary_fit$coefficients[2, 4],
        r_sq      = summary_fit$r.squared
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      slope_chr      = formatC(abs(slope), digits = 4, format = "f"),
      int_chr        = formatC(intercept,  digits = 4, format = "f", flag = "+"),
      mean_chr       = formatC(mean_y,     digits = 3, format = "f"),
      sd_chr         = formatC(sd_y,       digits = 3, format = "f"),
      ci_l_chr       = formatC(ci_l,       digits = 5, format = "f"),
      ci_u_chr       = formatC(ci_u,       digits = 5, format = "f"),
      r2_chr         = formatC(r_sq,       digits = 2, format = "f"),
      # Create strings for plotmath expressions.
      eqn_str = sprintf("p[x](t) == %s * ' * ' * N[EH]* %s", slope_chr, int_chr),
      stats_r2_str   = sprintf("R^2 == %s", r2_chr),
      stats_mean_str = sprintf("bar(p[x](t)) == %s * \"\\u00B1\" * %s", mean_chr, sd_chr),
      # Note: '~~' creates a space in plotmath expressions.
      ci_str         = sprintf("\"95%% CI of slope:\"~~'['*%s*','~~%s*']'", ci_l_chr, ci_u_chr),
      pval_str       = ifelse(pval < 0.001, "\"p < 0.001\"", sprintf("p == %.3f", pval))
    )
  
  facet_theme <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text       = ggplot2::element_text(size = 8, face = "bold"),
      axis.title       = ggplot2::element_text(size = 10),
      axis.text        = ggplot2::element_text(size = 9)
    )
  
  # NOTE: The hjust/vjust values for geom_text are manually tuned for the
  # specific plot dimensions and may need adjustment if dimensions change.
  p <- ggplot2::ggplot(df, ggplot2::aes(NEH, .data[[y_var]])) +
    ggplot2::geom_point(size = 1.0, alpha = 0.80, colour = "black") +
    ggplot2::geom_smooth(method = "lm",
                         linewidth = 1, colour = "#E15759", se = FALSE) +
    # Annotations for top facets (Times 0-6 days)
    ggplot2::geom_text(data = dplyr::filter(reg_stats, !(Time %in% c("14", "21"))),
                       ggplot2::aes(x = -Inf, y = Inf, label = eqn_str),
                       hjust = -0.07, vjust = 1.55, parse = TRUE, size = 1.7) +
    ggplot2::geom_text(data = dplyr::filter(reg_stats, !(Time %in% c("14", "21"))),
                       ggplot2::aes(x = -Inf, y = Inf, label = stats_mean_str),
                       hjust = -0.08, vjust = 2.75, parse = TRUE, size = 1.7) +
    ggplot2::geom_text(data = dplyr::filter(reg_stats, !(Time %in% c("14", "21"))),
                       ggplot2::aes(x = -Inf, y = Inf, label = ci_str),
                       hjust = -0.04, vjust = 5.15, parse = TRUE, size = 1.7) +
    ggplot2::geom_text(data = dplyr::filter(reg_stats, !(Time %in% c("14", "21"))),
                       ggplot2::aes(x = -Inf, y = Inf, label = pval_str),
                       hjust = -0.13, vjust = 6.30, parse = TRUE, size = 1.7) +
  # Annotations for bottom facets (Times 14 & 21 days), positioned differently.
  ggplot2::geom_text(data = dplyr::filter(reg_stats, Time %in% c("14", "21")),
                     ggplot2::aes(x = -Inf, y = -Inf, label = eqn_str),
                     hjust = -0.06, vjust = -6.85, parse = TRUE, size = 1.7) +
  
    ggplot2::geom_text(data = dplyr::filter(reg_stats, Time %in% c("14", "21")),
                       ggplot2::aes(x = -Inf, y = -Inf, label = stats_mean_str),
                       hjust = -0.09, vjust = -4.35, parse = TRUE, size = 1.7) +
    ggplot2::geom_text(data = dplyr::filter(reg_stats, Time %in% c("14", "21")),
                       ggplot2::aes(x = -Inf, y = -Inf, label = ci_str),
                       hjust = -0.05, vjust = -3.85, parse = TRUE, size = 1.7) +
    ggplot2::geom_text(data = dplyr::filter(reg_stats, Time %in% c("14", "21")),
                       ggplot2::aes(x = -Inf, y = -Inf, label = pval_str),
                       hjust = -0.19, vjust = -2.6, parse = TRUE, size = 1.7) +
    ggplot2::facet_wrap(~ Time, nrow = 3, dir = "v",
                        labeller = ggplot2::as_labeller(time_labels)) +
    ggplot2::labs(
      x = expression(
        "Number of exchangeable hydrogens (" *
          N[plain("EH")] *
          ")"),
      y = y_lab
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(5)) +
    facet_theme
  
  # Save the plot in both PDF and PNG formats.
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE) # Ensure output directory exists
  ggplot2::ggsave(file.path(out_dir, paste0(tag, ".pdf")),
                  p, device = cairo_pdf,
                  width = w_mm, height = h_mm, units = units)
  ggplot2::ggsave(file.path(out_dir, paste0(tag, ".png")),
                  p, device = ragg::agg_png, dpi = 1200,
                  width = w_mm, height = h_mm, units = units)
  invisible(p)
}


#' Robustly read a CSV file with potentially duplicated or messy headers.
#'
#' This function reads a CSV file, skipping a specified number of initial lines.
#' It is designed to handle cases where column headers might be duplicated by
#' making them unique using `make.unique`. It also attempts to convert columns
#' to numeric types where appropriate.
#'
#' @param file_path The path to the CSV file.
#' @param skip The number of lines to skip at the beginning of the file.
#' @return A data frame with the contents of the CSV file. If an error occurs,
#'   an empty data frame is returned and a message is printed.
#' @note For very large files, consider using `readr::read_csv` or
#'   `data.table::fread` for better performance.
#' @importFrom utils read.csv
#' @importFrom stats na.omit
read_csv_robust <- function(file_path, skip = 0) {
  tryCatch({
    # Read the file content
    lines <- readLines(file_path)
    
    # Skip specified number of lines
    data_lines <- lines[(skip + 1):length(lines)]
    
    # Read the header
    header <- unlist(strsplit(data_lines[1], ","))
    header <- make.unique(trimws(header))  # Ensure unique column names
    
    # Read the data
    data <- read.csv(text = paste(data_lines[-1], collapse = "\n"), 
                     header = FALSE, 
                     stringsAsFactors = FALSE, 
                     check.names = FALSE)
    
    # Assign column names
    colnames(data) <- header
    
    # Convert numeric columns efficiently
    data[] <- lapply(data, function(x) {
      num_x <- suppressWarnings(as.numeric(x))
      if(all(!is.na(num_x))) num_x else x
    })
    
    return(data)
  }, error = function(e) {
    message(paste("Error reading", basename(file_path), ":", e$message))
    return(data.frame())
  })
}


#' Robustly read a CSV file with a custom routine for unique headers.
#'
#' This function reads a CSV file, similar to `read_csv_robust`, but uses a
#' custom implementation to create unique column names by appending suffixes
#' like ".0", ".1", etc., to duplicates.
#'
#' @param file_path The path to the CSV file.
#' @param skip The number of lines to skip at the beginning of the file.
#' @return A data frame. Returns an empty data frame on error.
#' @note This is a custom implementation for making headers unique. The standard
#'   `make.unique()` function (used in `read_csv_robust`) is generally sufficient.
#' @importFrom utils read.csv
#' @importFrom stats na.omit
read_csv_robust2 <- function(file_path, skip = 0) {
  tryCatch({
    # Read file content and process headers
    lines <- readLines(file_path)
    data_lines <- lines[(skip + 1):length(lines)]
    
    # Process column headers
    header <- trimws(unlist(strsplit(data_lines[1], ",")))
    freq <- table(header)
    new_header <- character(length(header))
    counters <- list()
    
    for (i in seq_along(header)) {
      col_name <- header[i]
      if (freq[col_name] > 1) {
        if (!col_name %in% names(counters)) {
          counters[[col_name]] <- 0
        } else {
          counters[[col_name]] <- counters[[col_name]] + 1
        }
        new_header[i] <- paste0(col_name, ".", counters[[col_name]])
      } else {
        new_header[i] <- col_name
      }
    }
    
    # Read data with processed headers
    data <- read.csv(
      text = paste(data_lines[-1], collapse = "\n"),
      header = FALSE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    colnames(data) <- new_header
    
    # Convert numeric columns
    data[] <- lapply(data, function(x) {
      num_x <- suppressWarnings(as.numeric(x))
      if (all(!is.na(num_x))) num_x else x
    })
    
    return(data)
  }, error = function(e) {
    message(paste("Error reading", basename(file_path), ":", e$message))
    return(data.frame())
  })
}


#' Calculate the Relative Isotopic Abundance (RIA) from experimental data.
#'
#' This function calculates the RIA for a given peptide at a specific time point.
#' It expects isotopomer intensity columns named in a "I0.t", "I1.t", ... pattern.
#' It attempts to find data for two consecutive time indices (`time_index` and
#' `time_index + 1`), calculates RIA for each if possible, and returns their average.
#' If only one is available, it returns that one.
#'
#' @param time_index The base time index (integer) to look for (e.g., 0 for T0).
#'   The function will also look for `time_index + 1`.
#' @param index The row index in `merged_data` for the peptide of interest.
#' @param merged_data A data frame containing isotopomer intensity columns.
#' @return The calculated RIA as a numeric value, or `NA` if it cannot be computed.
#' @note This function assumes the M0 isotopomer (I0) is used for RIA calculation.
#'   It also writes to a global `warning_ria_log_file` which must be defined elsewhere.
get_ria <- function(time_index, index, merged_data) {
  # Dynamically construct column names for two consecutive time points
  i_cols1 <- paste0("I", 0:5, ".", time_index)
  i_cols2 <- paste0("I", 0:5, ".", time_index+1)

  cols1_exist <- all(i_cols1 %in% colnames(merged_data))
  cols2_exist <- all(i_cols2 %in% colnames(merged_data))

  if (!cols1_exist && !cols2_exist) {
    msg <- paste("Missing all RIA columns for time_index", time_index, "and", time_index + 1, "at row index", index)
    warning(msg)
    # This assumes warning_ria_log_file is a global variable.
    try(write(msg, file = warning_ria_log_file, append = TRUE), silent = TRUE)
    return(NA_real_)
  }

  # Helper to calculate RIA from a vector of isotopomer intensities
  calculate_ria_single <- function(values) {
    values <- values[!is.na(values)]
    # Ensure there are values and their sum is positive to avoid division by zero
    if (length(values) > 0 && sum(values) > 0) {
      return(values[1] / sum(values)) # RIA = I0 / sum(I0..In)
    }
    return(NA_real_)
  }

  i1 <- NA_real_
  if (cols1_exist) {
    i_values1 <- as.numeric(merged_data[index, i_cols1])
    i1 <- calculate_ria_single(i_values1)
  }

  i2 <- NA_real_
  if (cols2_exist) {
    i_values2 <- as.numeric(merged_data[index, i_cols2])
    i2 <- calculate_ria_single(i_values2)
  }

  # Average the valid RIAs. `mean` with `na.rm=T` handles all cases.
  mean_ria <- mean(c(i1, i2), na.rm = TRUE)

  # mean() of an all-NA vector is NaN. Return NA instead.
  return(ifelse(is.nan(mean_ria), NA_real_, mean_ria))
}


#' Compute the fraction of newly synthesized protein, Px(t).
#'
#' @param time_index The time index passed to `get_ria`.
#' @param index The row index in `merged_data`.
#' @param neh The Number of Exchangeable Hydrogens for the peptide.
#' @param merged_data The main data frame containing `M0` and intensity columns.
#' @return The calculated Px(t) value, or `NA`.
#' @note This function depends on a global variable `ph` (body water enrichment)
#'   which must be defined in the calling environment.
get_pxt <- function(time_index, index, neh, merged_data) {
  i0_0_teo <- (merged_data$M0[index] / 100)
  i0_t_exp <- get_ria(time_index, index, merged_data)

  # Ensure all components are valid numbers before calculation.
  if (!is.na(i0_0_teo) && !is.na(i0_t_exp) && i0_0_teo != 0 && !is.na(neh) && neh > 0) {
    # The `ph` variable is expected to be defined in the global environment.
    return((1-ph)*(1 - ((i0_t_exp / i0_0_teo)^(1 / neh))))
  } else {
    return(NA_real_)
  }
}


#' Fit an exponential decay curve to estimate protein turnover rate (k).
#'
#' This function uses non-linear least squares (nls) to fit the model:
#' `y ~ cf + (1 - cf) * exp(-k * x)`
#' where `y` is the adjusted Px(t), `x` is time, and `cf` is the 
#' the ratio of I0_asymptote to the monoisotopic RIA
#'
#' @param cf the ratio of I0_asymptote to the monoisotopic RIA
#' @param x A numeric vector of time points.
#' @param y A numeric vector of corresponding (adjusted) Px(t) values.
#' @return A list containing the estimated rate `k`, the Residual Sum of
#'   Squares `RSS`, and the R-squared `rsq_k`. Returns `NA` for these values
#'   if the fit fails.
#' @importFrom stats nls predict coef resid
d2ome_curve_fit_adj_pxt <- function(cf, x, y) {
  # Input validation
  if (length(x) != length(y) || length(x) < 2 || any(is.na(x)) || any(is.na(y))) {
    warning("Invalid input: x and y must be same length, at least length 2, and contain no NA.")
    return(list(k = NA_real_, RSS = NA_real_, rsq_k = NA_real_))
  }
  
  # Handle cases where y doesn't vary
  if (length(unique(y)) < 2) {
    warning("Y values are constant — curve fitting is not possible.")
    return(list(k = NA_real_, RSS = NA_real_, rsq_k = NA_real_))
  }
  
  # Attempt curve fitting with error handling
  fit_result <- tryCatch({
    model <- nls(
      y ~ cf + (1 - cf) * exp(-k * x),
      start = list(k = 0.01) # A reasonable starting value for k
    )

    # Compute R-squared for the non-linear model
    y_pred <- predict(model)
    y_obs  <- y
    ss_res <- sum((y_obs - y_pred)^2)
    ss_tot <- sum((y_obs - mean(y_obs))^2)
    rsq    <- if (ss_tot > 0) 1 - (ss_res / ss_tot) else 1

    k_est <- coef(model)["k"]
    rss <- sum(resid(model)^2)
    return(list(k = k_est, RSS = rss, rsq_k = rsq))
  }, error = function(e) {
    warning(paste("Curve fitting failed:", e$message))
    return(list(k = NA_real_, RSS = NA_real_, rsq_k = NA_real_))
  })

  return(fit_result)
}
