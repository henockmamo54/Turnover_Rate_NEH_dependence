
library(dplyr)


##### Compute slope and intercept for each time point via linear regression and adjust Pₓ(t)######
adjust_pxt<-function(pxt_df)
{

  model_params <- pxt_df %>%
  group_by(Time) %>%
  do({
    fit <- lm(Pxt ~ NEH, data = .)
    data.frame(
      slope     = round(coef(fit)[2], 4),
      intercept = round(coef(fit)[1], 4)
    )
  })


  pxt_data_joined <- pxt_df %>%
  left_join(model_params, by = "Time") %>%
  mutate(Pxt_adj = Pxt - NEH * slope)


  return( pxt_data_joined)
}


##### Compute slope and intercept for each time point via linear regression and adjust Pₓ(t) for each protein separately.######
 
adjust_pxt_all<-function(pxt_df)
{

  model_params <- pxt_df %>%
  group_by(Time,SourceFile) %>%
  do({
    fit <- lm(Pxt ~ NEH, data = .)
    data.frame(
      slope     = round(coef(fit)[2], 4),
      intercept = round(coef(fit)[1], 4)
    )
  })


  pxt_data_joined <- pxt_df %>%
  left_join(model_params, by = c("Time","SourceFile")) %>%
  mutate(Pxt_adj = Pxt - NEH * slope)


# }
  return( pxt_data_joined)
}

############### Create Px(t) plot based on 'pxt_data_adj' created by 'compute_turnover_rate.R' ##############
make_pxt_facet <- function(df, y_var, y_lab, tag,
                           out_dir = ".", w_mm = 180, h_mm = 120,
                           units = "mm") {
  time_labels <- c(
    `0` ="No labelling", `1` ="1-day",  `2` ="2-day",  `3` ="3-day",
    `4` ="4-day",        `5` ="5-day",  `6` ="6-day",
    `14`="13-day",       `21`="21-day"
  )
  df <- dplyr::mutate(df, Time = factor(Time, levels = names(time_labels)))
  
  reg_stats <- df |>
    dplyr::group_by(Time) |>
    dplyr::reframe({
      fit <- lm(get(y_var) ~ NEH)
      summary_fit <- summary(fit)
      tibble::tibble(
        slope     = coef(fit)[2],
        intercept = coef(fit)[1],
        ci_l      = confint(fit, level = 0.95)[2, 1],
        ci_u      = confint(fit, level = 0.95)[2, 2],
        mean_y    = mean(get(y_var)),
        sd_y      = sd(get(y_var)),
        pval      = summary_fit$coefficients[2, 4]
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
      # eqn_str        = sprintf("p[x](t) == %s * N[EH]* %s", slope_chr, int_chr),
      eqn_str = sprintf("p[x](t) == %s * ' * ' * N[EH]* %s", slope_chr, int_chr),
      
      stats_r2_str   = sprintf("R^2 == %s", r2_chr),
      stats_mean_str = sprintf("bar(p[x](t)) == %s * \"\\u00B1\" * %s", mean_chr, sd_chr),
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
  
  p <- ggplot2::ggplot(df, ggplot2::aes(NEH, .data[[y_var]])) +
    # ggplot2::geom_point(size = 1.8, alpha = 0.80, colour = "black") +
    ggplot2::geom_point(size = 1.0, alpha = 0.80, colour = "black") +
    ggplot2::geom_smooth(method = "lm",
                         linewidth = 1, colour = "#E15759", se = FALSE) +
    # Top-row facets (Times 0-6)
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
    ##################################
  # Bottom-row facets (Times 14 & 21)
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
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(out_dir, paste0(tag, ".pdf")),
                  p, device = cairo_pdf,
                  width = w_mm, height = h_mm, units = units)
  ggplot2::ggsave(file.path(out_dir, paste0(tag, ".png")),
                  p, device = ragg::agg_png, dpi = 1200,
                  width = w_mm, height = h_mm, units = units)
  invisible(p)
}

################################


######## This script reads the dataset files, merges them, and stores them in a DataFrame called 'merged'.############

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
###########################################
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


##############################
#######################################################################################


get_ria <- function(time_index, index, merged_data,imers=0) {
  # Construct column names dynamically
  i_cols1 <- paste0("I", 0:5, ".", time_index)
  i_cols2 <- paste0("I", 0:5, ".", time_index+1)
  # Check if all expected columns exist in merged_data
  if ((!all(i_cols1 %in% colnames(merged_data))) &&(!all(i_cols2 %in% colnames(merged_data)))) {
    msg<-paste("Missing RIA columns for time_index", time_index, "at index", index)
    warning(msg)
    write(msg, file = warning_ria_log_file , append = TRUE)
    return(NA)  # Return NA if expected columns are missing
  }
  
  # Extract values, suppress warnings about NAs introduced by coercion
  i_values1 <- as.numeric(merged_data[index, i_cols1])
  i_values2 <- as.numeric(merged_data[index, i_cols2])
  # Remove NAs and calculate RIA
  i_values1 <- i_values1[!is.na(i_values1)]
  i_values2 <- i_values2[!is.na(i_values2)]
  
  i_values_avg<-0
  f<-0
  f1<-0
  f2<-0
  if ((sum(i_values1, na.rm = TRUE)>0)&& length(i_values1) > 0 && !is.na(i_values1[1])) {
    i1<-(i_values1[1] / sum(i_values1, na.rm = TRUE))
    
    #return(ifelse(imers==0,(i_values[1] / sum(i_values, na.rm = TRUE)),(i_values[imers+1] / sum(i_values, na.rm = TRUE))))
    #return(i_values[1] / sum(i_values, na.rm = TRUE))
  } else {
    f<-f+1
    f1<-1
    #return(NA)  # Return NA if all values are NA
  }
  if ((sum(i_values2, na.rm = TRUE)>0)&&length(i_values2) > 0 && !is.na(i_values2[1])) {
    i2<-(i_values2[1] / sum(i_values2, na.rm = TRUE))
    
    #return(ifelse(imers==0,(i_values[1] / sum(i_values, na.rm = TRUE)),(i_values[imers+1] / sum(i_values, na.rm = TRUE))))
    #return(i_values[1] / sum(i_values, na.rm = TRUE))
  } else {
    f<-f+1
    f2<-1
    #return(NA)  # Return NA if all values are NA
  }
  if (f == 2) {
    return(NA)
  } else if (f == 0) {
    return((i1 + i2) / 2)
  } else if (f1 == 0) {
    return(i1)
  } else if (f2 == 0) {
    return(i2)
  }
}


#########################################################
## compute Px(t) #####
get_pxt <- function(time_index, index, neh, merged_data) {
  i0_0_teo <- (merged_data$M0[index] / 100)
  # i0_0_teo <-compute_M0(merged_data$Peptide[index])
  i0_t_exp <- get_ria(time_index, index, merged_data,0)
  #return((1 - (i0_t_exp / i0_0_teo)^(1 / neh)))
  if(!is.nan(i0_0_teo) &&!is.nan(i0_t_exp) && !(i0_0_teo==0)){
    return((1-ph)*(1 - ((i0_t_exp / i0_0_teo)^(1 / neh))))
  } else {
    return(NA)  # Return NA if m0 or i0(t) is NA
  }
}
######################################################### 
## Fitting function to estimate the turnover rate (k) ####
d2ome_curve_fit_adj_pxt <- function(cf, x, y,r) {
  # Input validation
  if (length(x) != length(y) || length(x) < 2 || any(is.na(x)) || any(is.na(y))) {
    warning("Invalid input: x and y must be same length, at least length 2, and contain no NA.")
    return(list(k = NA, RSS = NA))
  }
  
  # Handle cases where y doesn't vary
  if (length(unique(y)) < 2) {
    warning("Y values are constant — curve fitting is not possible.")
    return(list(k = NA, RSS = NA))
  }
  
  # Attempt curve fitting with error handling
  fit_result <- tryCatch({
    
    model <- nls(
      y ~ cf + (1 - cf) * exp(-k * x),
      start     = list(k = 0.00001)      # Starting value for k
     
    )
    ######### compute Rsqured for the model ###############
    y_pred <- predict(model)
    y_obs  <- y
    ss_res <- sum((y_obs - y_pred)^2)
    ss_tot <- sum((y_obs - mean(y_obs))^2)
    rsq    <- 1 - (ss_res / ss_tot)
    ############################################################
    k_est <- coef(model)["k"]
    rss <- sum(resid(model)^2)
    return(list(k = k_est, RSS = rss,rsq_k=rsq))
    
  }, error = function(e) {
    warning(paste("Curve fitting failed:", e$message))
    return(list(k = NA, RSS = NA,rsq_k=NA))
  })
  
  #return(fit_result)
}
############ End Helper Functions ##########
###########################################################################
