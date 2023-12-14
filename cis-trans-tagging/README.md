# pQTL Cis/Trans Tagger 

## Overview
This Python script is designed for the specific task of processing pQTL signals to determine if each SNP is in cis or trans position relative to target proteins. It does so by processing data from GWAS (Genome-Wide Association Studies) files, comparing these with a mapping file that contains information about target proteins and their locations.

## Requirements
- **Python 3.x**: The script is written in Python and requires Python 3.
- **Pandas**: Used for data manipulation and analysis.
- **mpmath**: A Python library for arbitrary-precision arithmetic, used for precise numerical operations.
- **JSON Configuration File**: A user-defined file that specifies various parameters and file paths needed for the script to run.

## Configuration File
A JSON configuration file is required to run the script. This file should contain the following parameters:
- **`map_file_path`**: Path to the CSV file that maps proteins to their genomic locations.
- **`gwas_files_dir`**: Directory containing the SNP files from GWAS studies.
- **`id_column_name`**: The column in the mapping file that uniquely identifies each protein.
- **Column Names**: Specific column names in the SNP and mapping files for processing (e.g., `snp_file_chrom_col`, `pval_column`).
- **`region_buffer_value`**: A numerical value that extends the genomic region of interest around each protein (default +/. 500K bases).
- **`pvalue_threshold`**: A threshold for filtering pQTL signals based on their p-value.
- **`columns_to_drop`**: A list of columns to be excluded from the final output.
- **`save_path`**: Destination directory for saving output and log files.
- **`output_filename`**: The filename for the final output in pickle format.

## Script Workflow

1. **Parsing Config File**: The script starts by reading the JSON configuration file to get all necessary parameters.
2. **Logging Setup**: It sets up logging to capture detailed information about the script’s execution and any warnings. These logs are saved in a log file named `cis_trans_tagging_log.txt` in the specified `save_path`.
3. **Processing Each SeqID**: For each unique SeqID (representing a protein or a group of proteins) in the mapping file:
   - Identifies the genomic location of the target protein.
   - Finds all corresponding SNP files.
   - If more than 2 SNP files are found for a SeqID, logs a warning and processes only the first file.
4. **SNP File Processing**: For each SNP file:
   - Filters SNPs based on the p-value threshold.
   - Tags each SNP as 'cis' or 'trans' relative to the target protein's location.
   - Calculates distances from the genomic region of interest (to flag pQTL signals that are trans, but "proximal" to the protein coding region).
6. **Combining Results**: All processed data for each SeqID are combined into a single DataFrame.
7. **Saving Output**: The final combined DataFrame is saved as a pickle file in the specified location with the filename provided in the configuration file. This file format is required to safely store p-values with decimals beyond floats memory allocation.

**NOTES** P-values in the GWAS results should be expressed in scientific notation (e.g. '1.2345e-678'). 

## Running the Script
To run the script, use the following command in the terminal:
```
python cis_trans_tagger.py <path_to_config_file>
```

