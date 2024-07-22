# Protein Residuals Script README

**Author:** Michela Carlotta Massi  
**Date:** 2023-12-04

## Overview

This R script is designed for generating residuals from linear regression models for association studies involving protein data (i.e. pQTL studies). The script utilizes various R packages, including `dplyr`, `data.table`, `tidyverse`, and `RJSONIO`, to preprocess and analyze protein datasets.

## Prerequisites

Before running the script, ensure you have the necessary R packages installed. You can install them using the following commands:

```R
install.packages(c("dplyr", "data.table", "tidyverse", "RJSONIO"))
```

## Experiment Configuration

The script relies on a JSON configuration file (`residuals_config.json`) that contains experiment-specific parameters. In this file, you can specify the paths to protein data, covariates data, genetic principal components (or supply one path for a dataset containing all information), and define the output path and filename for the results. The `conf_file_path` variable in the script points to this configuration file.

**NOTES** 
- Under `model_definition` the user should define the linear regression model that to extract residuals from. **In the definition, covariates and principal components should be called with the names of the columns in the datasets provided as input** (see Data Requirements).
Currently, NAs in any of the variables included in the model definition will result in `NA` for that subject in the residuals' dataframe.

- The key `id_cols` in the configuration file requires a list containing the names of all the ID columns that will be concatenated to the output of the script. Therefore, the naming should reflect columns available in at least one of the dataframes supplied to the script.

- The config file is nested to allow for extracting residuals from different models, or using different datasets. The structure of the configuration file should remain nested irrespectively of the number of experiments to run. In case of a single experiment, the output will be saved as defined in `results_filename`, otherwise, the script will attach to the output name of each experiment the name defined at the first level of the nested config file.

- The variable `add_suffix` in the configuration file defines whether the columns of the residuals will be saved as "protein_name_res" (if `add_suffix = 1`) or "protein_name" (otherwise).

```R
# Example of residuals_config.json
{
  "main": {
        "n_datasets" : 3,
        "proteins_data_path": "/path/to/protein_data.txt",
        "gen_PC_path": "/path/to/genPCs.txt",
        "covariates_data_path" : "/path/to/covariates.txt",
        "unique_dataset_path" : " ", # change in case n_datasets = 1
        "n_start_prot": None, # change in case n_datasets = 1
        "n_end_prot": None, # change in case n_datasets = 1
        "results_path": "/path/to/results/",
        "results_filename": "Resulting_residuals",
        "add_suffix": TRUE,
        "model_definition": "~  age + sex + PC1+PC2+PC3+PC4+PC5"
  },
  "additional_experiment": {
    // ... additional configurations ...
  }
}
```

## Data Requirements

The script supports two data configurations:

1. **Multiple Datasets (`n_datasets` = 3):**
 - **Proteins Data**: A .txt data file containing protein data. It should be formatted as a data frame **containing only proteins** with appropriate column names (i.e. names of proteins to be analyzed). The location and name of the file should be defined in `residuals_config.json` as `proteins_data_path`.

- **Covariates Data**: A .txt data file containing covariates data. It should be formatted as a data frame with the column names that need to match those in the model definition. The location and name of the file should be defined in `residuals_config.json` as `covariates_data_path`. **To ensure the linear model function in the script treats categorical covariates as such, it is important to use string-format categories (i.e. ["Jan", "Feb",...]. Numbers should be used for binary variables only (i.e. [0,1]).**

- **Genetic Principal Components**: A .txt file containing genetic principal components. The location and name of the file should be defined in `residuals_config.json` as `gen_PCs_path`.

**NOTE**: the script does not perform any ID matching. **The 3 files supplied in this case should contain the same samples in the same order**.
  
2. **Single Dataset (n_datasets = 1):**
   - Unique Dataset: Unique dataset (.txt file) containing all information needed to fit the linear model, including proteins, covariates, and genetic principal components.. This case requires the range of indices identifying the proteins in the dataframe. Provide it using the`n_start_prot` and `n_end_prot` parameters in the configuration file. 

## Functions

### inormal

The `inormal` function performs the inverse rank-normal transformation on a given vector.

### lm_reg_residuals

The `lm_reg_residuals` function conducts linear regression on a specified response variable (`Y`) using a provided dataset and model definition. The function returns the residuals of the regression.

## Execution

1. **Load Experiment Configurations**: Set the path to the experiment configuration file (`residuals_config.json`).

```R
conf_file_path = "residuals_config.json"
configs <- fromJSON(conf_file_path)
```

2. **Iterate through Experiments**: The script iterates through each experiment specified in the configuration file, performs the regression analysis, and saves the residuals as a text file.

3. **Data Transformation and Regression Analysis**: Protein data is inverse rank-normal transformed, combined with covariates and genetic principal components, and used as input for linear regression models. Residuals are obtained for each protein individually.

4. **Output Files**: Residuals are saved as text files with the specified filename in the results directory. If there are multiple configurations, the filename is appended with the experiment name. 

## Running the Script

1. Ensure R and the required packages are installed.
2. Create a valid experiment configuration file (`residuals_config.json`).
3. Modify the `conf_file_path` variable to point to the configuration file.
4. Execute the script.
