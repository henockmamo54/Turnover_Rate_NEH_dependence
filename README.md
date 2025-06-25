# Analysis of the dependence of protein turnover rates on the number of exchangeable hydrogens in deuterated water-labeled samples.

This project provides an R-based workflow for calculating protein turnover rates from mass spectrometry data, eliminating the dependence of the turnover rates on the Number of Exchangeable Hydrogens (NEH).

## Overview

The analysis pipeline performs the following key steps:

1.  **Data Ingestion**: Reads and merges `.Quant.csv` and `.RateConst.csv` files from a specified data directory.
2.  **px(t) Calculation**: Calculates the fraction of newly synthesized protein, px(t), for each peptide at multiple time points.
3.  **NEH-based px(t) Adjustment**: Performs a linear regression of px(t) against NEH for each time point and calculates an adjusted px(t) value to correct for systematic biases.
4.  **Turnover Rate Fitting**: Fits the adjusted px(t) data to a non-linear exponential decay model to estimate the protein turnover rate constant (`k`).
5.  **Output Generation**: Saves the final calculated turnover rates to a clean, analysis-ready CSV file.

The core of this method is the `adjust_pxt_all` function, which corrects for the observed dependency of px(t) on NEH, aiming to provide a more accurate estimation of turnover rates.

## Features

*   Automated read and merging of d2ome quantification outputs(`.Quant.csv` and `.RateConst.csv` files.)
*   Filtering of data based on R-squared thresholds.
*   Calculation and adjustment of Pxt values .
*   Non-linear curve fitting using a exponential decay model to determine peptide exchange rates.
*   Robust error handling and logging for failed curve fits.
*   Generation of a comprehensive output dataframe with calculated turnover rates.
  
## Project Structure

```
Turnover_Rate_NEH_dependence/
├── SampleData/
│   ├── sample1.Quant.csv
│   ├── sample1.RateConst.csv
│   └── ... (other data files)
├── results/
│   └── final_turnover_rates.csv  (Output file)
├── log/
│   ├── warning_ria_error.txt
│   └── warning_curve_fit_errors.txt
├── main.R                # Main script to run the entire analysis
├── Helper.R              # Contains all helper functions for the analysis
└── README.md             # This file
```

## Prerequisites

Before running the analysis, you need to have R installed on your system, along with the following R packages:

-   `dplyr`
-   `ggplot2`
-   `purrr`
-   `tidyr`
-   `readr`

You can install all the required packages by running the following command in your R console:

```r
install.packages(c("dplyr", "ggplot2", "purrr", "tidyr", "readr"))
```

## Setup and Usage

1.  **Clone the Repository**:
    ```bash
    git clone <repository-url>
    cd Turnover_Rate_NEH_dependence
    ```

2.  **Prepare Your Data**:
    -   Place your raw data files (`.Quant.csv` and `.RateConst.csv`) inside the `SampleData` directory.
    -   The script expects paired files, for example, `my_sample.Quant.csv` and `my_sample.RateConst.csv`.

3.  **Configure the Analysis (Optional)**:
    -   Open the `main.R` script.
    -   In the `--- CONFIGURATION ---` section at the top, you can adjust parameters such as:
        -   `folder_path`: The directory containing your input data.
        -   `output_dir`: The directory where results will be saved.
        -   `rsquared_threshold`: The minimum R-squared value from the initial rate constant calculation to include a peptide in the analysis.
        -   `exp_time`: A vector of the experimental time points.

4.  **Run the Analysis**:
    -   Open R or RStudio.
    -   Set your working directory to the root of the project folder.
    -   Source the main script to execute the entire workflow:
    ```r
    source("main.R")
    ```
    -   The script will print progress messages to the console and notify you upon completion.

## Output

The primary output is the `final_turnover_rates.csv` file located in the `results` directory. This file contains the calculated turnover rates for each peptide and includes the following columns:

-   `Protein`: The name of the source protein/sample file.
-   `Peptide`: The peptide sequence.
-   `Charge`: The charge state of the peptide.
-   `NEH`: The Number of Exchangeable Hydrogens.
-   `k_old_methold`: The turnover rate constant from the original `.RateConst.csv` file.
-   `Rsquared_old_methold`: The R-squared value from the original fit.
-   `k_new_methold`: The newly calculated turnover rate constant after NEH adjustment.
-   `Rsquared_new_methold`: The R-squared value for the new turnover rate fit.

Log files containing any warnings or errors encountered during the run are saved in the `log` directory.


### `Helper.R`

This external file is crucial and is expected to contain the following R functions:

*   `read_csv_robust2(file_path, skip)`: A robust CSV reading function.
*   `get_pxt(time_index, index, neh, merged_data)`: Function to extract Pxt values from the merged data.
*   `adjust_pxt_all(pxt_data)`: Function to perform adjustments on the Pxt values.
*   `d2ome_curve_fit_adj_pxt(cof, times, B_values, maxrate)`: Function to perform the non-linear curve fitting for turnover rate determination.

