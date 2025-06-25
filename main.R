# --- LIBRARIES ---
library(dplyr)
library(ggplot2)
library(purrr) # For functional programming tools like map()
library(tidyr) # For data manipulation functions like unnest() and nest()
library(readr) # For writing CSV files efficiently

# Input and Output Paths
folder_path <- "SampleData"  # Folder containing the input CSV files
log_dir     <- "log"         # Directory to store log files
output_dir  <- "results"     # Directory to store output files

# Log file names
warning_ria_log_file <- file.path(log_dir, "warning_ria_error.txt")
error_log_file <- file.path(log_dir, "warning_curve_fit_errors.txt")

# Output file name
final_results_file <- file.path(output_dir, "final_turnover_rates.csv")
protein_summary_file <- file.path(output_dir, "protein_turnover_summary.csv")
protein_summary_traditional_filter_file <- file.path(output_dir, "protein_summary_traditional_filter.csv")

# Analysis Parameters
rsquared_threshold <- 0.80
exp_time <- c(0, 1, 2, 3, 4, 5, 6, 14, 21) # Experimental time points
pw <- 0.0399  # Body water enrichment constant
ph <- 1.5574E-4 # A constant used in Px(t) calculation

# --- SOURCE HELPER FUNCTIONS ---
source("Helper.R")

# --- SETUP ---
# Create the log directory if it doesn't exist
if (!dir.exists(log_dir)) {
  dir.create(log_dir)
}


# --- 1. DATA LOADING AND PREPARATION ---
message("Starting data loading and preparation...")

# Find all ".Quant.csv" files in the specified folder.
quant_files <- list.files(folder_path, pattern = "\\.Quant\\.csv$", full.names = TRUE)
# Name the vector with the base file names to get meaningful IDs from map_dfr
names(quant_files) <- sub("\\.Quant\\.csv$", "", basename(quant_files))

# Use map_dfr to read, process, and row-bind all files in a single step.
# This is more efficient and readable than a for-loop with a list.
all_merged_data <- map_dfr(quant_files, function(quant_file) {
  rate_file <- sub("\\.Quant\\.csv$", ".RateConst.csv", quant_file)

  if (!file.exists(rate_file)) {
    message(paste("Skipping", basename(quant_file), "- corresponding RateConst file not found"))
    return(NULL) # map_dfr will ignore NULL returns
  }

  data_quant <- read_csv_robust2(quant_file, skip = 3)
  data_rate <- read_csv_robust2(rate_file, skip = 0)

  if (nrow(data_quant) == 0 || nrow(data_rate) == 0) {
    message(paste("Skipping", basename(quant_file), "- empty data detected"))
    return(NULL)
  }

  # Standardize column names for a consistent merge key
  data_rate <- rename(data_rate, Peptide = any_of("Peptides"))

  inner_join(data_quant, data_rate, by = c("Peptide", "Charge"))
}, .id = "SourceFile")

# Filter the merged data based on the R-squared threshold from the initial fit
merged_filtered <- all_merged_data %>%
  filter(as.numeric(Rsquared) >= rsquared_threshold)

message(sprintf("Loaded and merged data for %d peptides from %d files.",
                nrow(merged_filtered), length(quant_files)))


# --- 2. CALCULATE Px(t) FOR EACH PEPTIDE AT EACH TIME POINT ---
message("Calculating Px(t) for each peptide...")

# We build a list of data frames and bind them all at once.
pxt_list <- vector("list", nrow(merged_filtered))

for (i in seq_len(nrow(merged_filtered))) {
  peptide_row <- merged_filtered[i, ]

  # For each peptide, create a tibble of its Px(t) over the time course
  pxt_results_for_peptide <- tibble(Time = exp_time) %>%
    mutate(
      # The time_index formula seems specific to the dataset structure.
      time_index = 2 * (row_number() - 1),
      Pxt = map_dbl(time_index, ~get_pxt(
        time_index = .x,
        index = i, # Pass the row index in the original dataframe
        neh = peptide_row$NEH,
        merged_data = merged_filtered
      ))
    ) %>%
    filter(!is.na(Pxt) & is.finite(Pxt)) # Filter out failed calculations

  # If there are valid Px(t) values, add peptide metadata
  if (nrow(pxt_results_for_peptide) > 0) {
    pxt_list[[i]] <- pxt_results_for_peptide %>%
      mutate(
        SourceFile = peptide_row$SourceFile,
        peptide_name = peptide_row$Peptide,
        peptide_charge = peptide_row$Charge,
        peptide_len = nchar(peptide_row$Peptide),
        NEH = peptide_row$NEH,
        m0 = peptide_row$M0 / 100,
        k_ratio_csv = peptide_row$RateConstants,
        Rsquared_csv = peptide_row$Rsquared
      )
  }
}

# Combine the list of data frames into a single, tidy data frame
pxt_data <- bind_rows(pxt_list)


# --- 3. ADJUST Px(t) AND FIT TURNOVER RATE (k) ---
message("Adjusting Px(t) and fitting turnover rates...")

# Adjust Px(t) values by correcting for NEH-dependent effects.
pxt_data_adj <- adjust_pxt_all(pxt_data)

# Define a helper function to safely run the curve fit on a data subset
safe_d2ome_fit <- function(data, neh, pw, ph) {
  if (nrow(data) < 2) { # Fitting requires at least 2 points
    return(tibble(k = NA_real_, RSS = NA_real_, Rsquared_k = NA_real_))
  }
  # Calculate coefficients based on formulas from the source literature
  I0_asymp_I0_0 <- (1 - (pw / (neh * (1 - ph))))^neh # this value representens the ratio of I0_asymptote to the monoisotopic RIA
  B_values_adjust9 <- (1 - (data$Pxt_adj / (neh * (1 - ph))))^neh

  fit_result <- d2ome_curve_fit_adj_pxt(cf = I0_asymp_I0_0, x = data$Time, y = B_values_adjust9)

  tibble(k = fit_result$k, RSS = fit_result$RSS, Rsquared_k = fit_result$rsq_k)
}

# Use a split-apply-combine approach to fit each peptide's time course.
# This is a robust and idiomatic 'tidyverse' alternative to a for-loop.
fit_results_df <- pxt_data_adj %>%
  group_by(SourceFile, peptide_name, peptide_charge, NEH, k_ratio_csv, Rsquared_csv) %>%
  nest() %>% # Nest time-course data into a list-column
  mutate(fit_result = map(data, ~{
    tryCatch({
      safe_d2ome_fit(.x, neh = first(NEH), pw = pw, ph = ph)
    }, error = function(e) {
      msg <- sprintf("Error fitting peptide %s (charge %s) from file %s: %s", first(peptide_name), first(peptide_charge), first(SourceFile), e$message)
      warning(msg)
      write(msg, file = error_log_file, append = TRUE)
      tibble(k = NA_real_, RSS = NA_real_, Rsquared_k = NA_real_) # Return NA on error
    })
  })) %>%
  unnest(fit_result) %>% # Unnest results back into columns
  ungroup() %>%
  select(-data, -RSS) # Remove intermediate and unwanted columns

# Filter out peptides where the fitting failed
final_results <- fit_results_df %>%
  filter(!is.na(k)) %>%
  # Rename columns to the desired final format
  rename(
    Protein = SourceFile,
    Peptide = peptide_name,
    Charge = peptide_charge,
    k_old_methold = k_ratio_csv,
    Rsquared_old_methold = Rsquared_csv,
    k_new_methold = k,
    Rsquared_new_methold = Rsquared_k
  )

# --- 4. SAVE RESULTS ---
message(sprintf("Saving final results for %d peptides to '%s'...", nrow(final_results), final_results_file))

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Write the final results data frame to a CSV file
write_csv(final_results, final_results_file)

# --- 5. GENERATE AND SAVE PROTEIN-LEVEL SUMMARY ---
message("Generating protein-level summary...")

# Define the R-squared threshold for including peptides in the protein summary
protein_summary_r2_threshold <- 0.8

protein_summary <- final_results %>%
  group_by(Protein) %>%
  summarise(
    # Old method: filter by old R-squared and summarize old k
    peptides_used_old = sum(Rsquared_old_methold > protein_summary_r2_threshold, na.rm = TRUE),
    median_k_old = median(k_old_methold[Rsquared_old_methold > protein_summary_r2_threshold], na.rm = TRUE),
    sd_k_old = sd(k_old_methold[Rsquared_old_methold > protein_summary_r2_threshold], na.rm = TRUE),

    # New method: filter by new R-squared and summarize new k
    peptides_used_new = sum(Rsquared_new_methold > protein_summary_r2_threshold, na.rm = TRUE),
    median_k_new = median(k_new_methold[Rsquared_new_methold > protein_summary_r2_threshold], na.rm = TRUE),
    sd_k_new = sd(k_new_methold[Rsquared_new_methold > protein_summary_r2_threshold], na.rm = TRUE)
  ) %>%
  # Replace NaN with NA for cases where sd() is calculated on a single peptide
  mutate(across(where(is.numeric), ~if_else(is.nan(.), NA_real_, .))) %>%
  # Remove proteins that have no quantifiable peptides passing the threshold in either method
  filter(peptides_used_old > 0 | peptides_used_new > 0)

# Save the protein-level summary to a new CSV file
write_csv(protein_summary, protein_summary_file)
message(sprintf("Protein-level summary saved to '%s'.", protein_summary_file))

# --- 6. GENERATE PROTEIN SUMMARY USING ONLY TRADITIONAL R-SQUARED FILTER ---
message("Generating protein-level summary using only the traditional R-squared filter...")

protein_summary_traditional <- final_results %>%
  group_by(Protein) %>%
  summarise(
    # Use only peptides that pass the old R-squared threshold for all calculations
    peptides_used = sum(Rsquared_old_methold > protein_summary_r2_threshold, na.rm = TRUE),
    mean_k_old = mean(k_old_methold[Rsquared_old_methold > protein_summary_r2_threshold], na.rm = TRUE),
    sd_k_old = sd(k_old_methold[Rsquared_old_methold > protein_summary_r2_threshold], na.rm = TRUE),
    mean_k_new = mean(k_new_methold[Rsquared_old_methold > protein_summary_r2_threshold], na.rm = TRUE),
    sd_k_new = sd(k_new_methold[Rsquared_old_methold > protein_summary_r2_threshold], na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), ~if_else(is.nan(.), NA_real_, .))) %>%
  filter(peptides_used > 0)

# Save the new summary to a separate CSV file
write_csv(protein_summary_traditional, protein_summary_traditional_filter_file)
message(sprintf("Protein summary with traditional filter saved to '%s'.", protein_summary_traditional_filter_file))

message("Script finished.")
