# Turnover Rates and Numbers of Exchangeable Hydrogens in Deuterated Water Labeled Samples 

This R script automates the computational approach to eliminate the dependence of the turnover rates on the number of exchangeable hydrogens. It reads quantitative peptide data and corresponding rate constant information, merges them, calculates deuterium uptake over time (Pxt), and then performs non-linear curve fitting to determine the kinetic exchange rates (k) for individual peptides. 

## Features

*   Automated read and merging of d2ome quantification outputs(`.Quant.csv` and `.RateConst.csv` files.)
*   Filtering of data based on R-squared thresholds.
*   Calculation and adjustment of Pxt values .
*   Non-linear curve fitting using a exponential decay model to determine peptide exchange rates.
*   Robust error handling and logging for failed curve fits.
*   Generation of a comprehensive output dataframe with calculated turnover rates.

## Prerequisites

*   **R**: A recent version of R (e.g., 4.x.x) is recommended.
*   **R Packages**: The following packages are required. You can install them using:
    ```R
    install.packages(c("dplyr", "minpack.lm", "nlme", "MASS", "ggplot2", "readr"))
    ```
    *   `dplyr` for data manipulation.
    *   `minpack.lm` for non-linear least squares fitting.
    *   `nlme` for linear and non-linear mixed effects models (though not explicitly used in the provided snippet, it's loaded).
    *   `MASS` for various statistical functions (loaded but not explicitly used in the snippet).
    *   `ggplot2` for data visualization (loaded but not explicitly used for plotting in the snippet).
    *   `readr` (implied by `read_csv_robust2` in `Helper.R`) for efficient CSV reading.

## Setup and Usage

1.  **Place Files**: Ensure the main script (`main.R`) and the `Helper.R` file are in the same directory.

2.  **Prepare Your Data**:
    *   Organize your `.Quant.csv` and `.RateConst.csv` files within a designated input folder.
    *   For every `.Quant.csv` file, there must be a corresponding `.RateConst.csv` file with the same base name (e.g., `Protein1.Quant.csv` and `Protein1.RateConst.csv`).

3.  **Configure the Script**:
    *   Open the main R script in an R editor (e.g., RStudio).
    *   **Crucially, set the `folder_path` variable** to the directory containing your input data files:
        ```R
        folder_path <- "path/to/your/data/folder" # <-- IMPORTANT: Set this path!
        ```
    *   Adjust other parameters as needed for your specific experiment:
        *   `rsquared_threshold`: Minimum R-squared value for initial data filtering (default: 0.80).
        *   `exp_time`: A numeric vector defining your experimental time points (e.g., `c(0, 1, 2, 3, 4, 5, 6, 14, 21)`).
        *   `ph=1.5574E-4`, `pw`: Parameters related to water and hydrogen concentrations/properties (defaults provided are for a specific dataset, adjust if necessary).
    *   (Optional) Modify the paths for `warning_ria_log_file` and `error_log_file` if you want logs saved elsewhere.

4.  **Run the Script**:
    Execute the script from your R environment:
    ```R
    source("main.R") 
    ```
    The script will process all `.Quant.csv` and corresponding `.RateConst.csv` files found in the `folder_path`.

5.  **Review Output**:
    The primary output will be stored in the `results_adj_df` dataframe within your R session. You can then save this dataframe to a CSV file or further analyze it:
    ```R
    write.csv(results_adj_df, "results_adj_df.csv", row.names = FALSE)
    ```

## Input Data Format

The script expects two types of CSV files for each protein:

*   **`.Quant.csv`**: Contains quantitative peptide data.
    *   The script uses `read_csv_robust2(..., skip = 3)`, implying the first 3 rows are headers/metadata to be skipped.
    *   It should contain columns like `"Peptide"`, `"Charge"`, and time-point specific deuterium uptake values (e.g., `X0`, `X1`, `X2`, etc., corresponding to your `exp_time` points).

*   **`.RateConst.csv`**: Contains rate constant information.
    *   The script uses `read_csv_robust2(..., skip = 0)`.
    *   It should contain columns like `"Peptides"` (which is renamed to `"Peptide"`), `"Charge"`, `"RateConstants"`, `"Rsquared"`, `"RootMeanRSS"`, `"NDP"`, `"M0"`, `"NEH"`.

### `Helper.R`

This external file is crucial and is expected to contain the following R functions:

*   `read_csv_robust2(file_path, skip)`: A robust CSV reading function.
*   `get_pxt(time_index, index, neh, merged_data)`: Function to extract Pxt values from the merged data.
*   `adjust_pxt_all(pxt_data)`: Function to perform adjustments on the Pxt values.
*   `ElementalComposition(peptide_sequence)`: Function to calculate the number of exchangeable hydrogens (NH) from a peptide sequence.
*   `d2ome_curve_fit_adj_pxt(cof, times, B_values, maxrate)`: Function to perform the non-linear curve fitting for turnover rate determination.

## Output Data (`results_adj_df`)

The `results_adj_df` dataframe contains the calculated turnover parameters for each processed peptide. Key columns include:

*   `SourceFile`: Original input file name.
*   `peptide_name`: Name of the peptide.
*   `k`: The calculated turnover exchange rate constant (turnover rate) from the curve fit.
*   `RSS`: Residual Sum of Squares from the curve fit.
*   `neh`: Number of exchangeable hydrogens for the peptide.
*   `rate`: The original rate constant from the input `.RateConst.csv` file.
*   `Rsquared_k`: R-squared value for the turnover curve fit.
*   `peptide_charge`: Charge state of the peptide.
*   `RootMeanRSS`, `NDP`, `RSS_csv`, `Rsquared_csv`: Metrics from the initial `.RateConst.csv` file.

## Error Handling and Logging

The script includes `tryCatch` blocks to gracefully handle errors during file processing and curve fitting. Errors are logged to the specified `error_log_file` (default: `/log/warning_curve_fit_errors.txt`), helping in debugging and identifying problematic data points.
