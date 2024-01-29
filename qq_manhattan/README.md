# Bulk Stratified QQ-plots and Manhattan Plots

## Overview
This Python script is designed for processing protein Genome-Wide Association Study (GWAS) summary statistics, focusing on generating (Stratified) QQ-plots and Manhattan plots using the GWASlab package.
The script is designed to parallelize plots' generation with an externally defined number of parallel jobs.

**Note**: The script is currently tailored to handle summary statistics produced by the **Regenie** algorithm. Future versions of the script will expand compatibility to include data from other software.

## Configuration
The script is configured through a JSON file. Each key in the JSON file represents a parameter used in the script.

### Configuration File (`plots_config.json`):
The configuration file should have the following structure:

```json
{
    "gwas_files_list_dir": "/path/to/summary_stats_list.txt",
    "save_path": "/path/to/output/plots/",
    "scaled": true,
    "build": "19",
    "skip": 1,
    "cut_dist": 2,
    "sig_level": 5e-8,
    "maf_bins": [[0, 0.01], [0.01, 0.05], [0.05, 1]],
    "stratified": true,
    "n_processes": 4
}
```

#### Parameters:
- **gwas_files_list_dir**: Path to a text file listing the GWAS summary statistics file directories.
- **save_path**: Directory to save the generated plots.
- **scaled**: Boolean to determine if `MLOG10P` is used for plotting.
- **build**: Genome build version (e.g., "19" for GRCh37).
- **skip**: Value to skip variants with `MLOG10P` less than this for faster plotting.
- **cut_dist**: Distance used to calculate the `cut` value for QQ plots.
- **sig_level**: Significance level for the QQ plot.
- **maf_bins**: Bins for Minor Allele Frequency (MAF) stratification.
- **stratified**: Boolean to generate a MAF stratified QQ plot.
- **n_processes**: Number of parallel processes for execution.

## Script Functionality
The script reads GWAS summary statistics files listed in `gwas_files_list_dir` and processes each file in parallel. For each file, it generates a QQ plot and/or Manhattan plot based on the configurations set in the JSON file. The output plots are saved in the directory specified by `save_path`.

## Note on automatic truncation in Manhattan Plots
The "cut" value is used to set the level at which Manhattan plots should be truncated along the y-axis (to allow for clear visualization in the case of outlier pQTLs with extremely high MLOG10P values) 

**The script calculates the "cut" value dynamically for each GWAS summary statistics file it processes.**
- It uses the significance level (`sig_level`) from the configuration file and a distance value (`cut_dist`), also specified in the configuration.
- The calculation is as follows: cut = round(-np.log10(`sig_level`)) + `cut_dist`.
This calculation means that the "cut" value is set to a certain number of orders of magnitude above the significance threshold (defined by `sig_level`), where the number of orders is determined by `cut_dist`.


## Usage
To run the script, use the following command in the terminal:

```bash
python bulk_qqmanh_maker.py plots_config.json
```

## Required Modules:
1. **gwaslab**: The core library used for GWAS analysis.
2. **numpy**: Provides support for large, multi-dimensional arrays and matrices, along with a collection of mathematical functions to operate on these arrays.
3. **argparse**: Used for parsing command line arguments.
4. **json**: Required for parsing and handling JSON files.
5. **multiprocessing**: Provides support for concurrent execution using processes.

### Installation Instructions:
You can install these modules using `pip`, which is the package installer for Python. Open your terminal or command prompt and run the following command:
```bash
pip install gwaslab numpy argparse json multiprocessing 
```
Note that some of these modules, like `argparse`, `json`, and `multiprocessing`, are part of the standard Python library and should be available in any standard Python installation. Therefore, they might not need to be installed separately.

## Additional Notes
- For more detailed information on GWASlab's capabilities, refer to the [GWASlab Documentation](https://cloufield.github.io/gwaslab/Visualization/).

