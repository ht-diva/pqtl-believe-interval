
# pQTL Cis/Trans Tagger

## Overview
This Python script processes pQTL (protein Quantitative Trait Loci) signals to categorize each SNP as being in cis or trans position relative to target proteins. It analyzes data from GWAS summary statistics' Clumping or Conditional Analysis, annotating the input tables with cis/trans tagging and the nearest gene to each SNP. 
To do that, it references a mapping file containing information about target proteins and their genomic locations, and the RefSeq `glist-hg19.txt` gene lists ([source](https://www.cog-genomics.org/plink/1.9/resources), generated from UCSC Table Browser RefSeq track data in May 2014) to annotate trans signals. The use of other reference gene lists has not been implemented yet. Future versions may include it, if needed. 

-----------------

## Requirements
- **Python 3.x**: Required for executing the script.
- **Pandas**: For data manipulation and analysis.
- **mpmath**: Used for arbitrary-precision arithmetic in numerical operations, to make sure the pvalues are reported w/o unwanted numerical approximations.
- **Pybedtools**: A Python tool used for genomic intervals manipulation and analysis.
- **JSON Configuration File**: Defines parameters and file paths for the script. This folder contains two versions, one prefilled with column names from plink Clumping pipeline, the other with column names from the Conditional Analysis pipeline (with COJO).

## Configuration File
The JSON configuration file must include:
- **`map_file_path`**: Path to a CSV file mapping proteins to genomic locations.
- **`gwas_files_dir`**: Directory containing SNP files from GWAS.
- **`gene_pos_file_dir`**: the directory to the `glist-hg19.txt` file.
- **`id_column_name`**: Column in the mapping file that uniquely identifies each protein.
- **Column Names**: Specific names of columns in SNP and mapping files required for processing (e.g., `snp_file_chrom_col`, `pval_column`).  
- **`region_buffer_value`**: Numerical value extending the genomic region of interest around each SNP to flag cis signals.
- **`pvalue_threshold`**: Threshold for filtering significant pQTL signals.
- **`columns_to_drop`**: List of columns from input files to exclude from the final output (optional, an empty list can be supplied).
- **`save_path`**: Directory for output and log files.
- **`output_filename`**: Name for the output file in pickle format.

----------------------

## Script Workflow

1. **Config File Parsing**: Reads parameters from the JSON configuration file.
2. **Logging Setup**: Initiates logging for detailed execution information.
3. **SeqID Processing**: Processes each unique SeqID:
   - Identifies the target protein's genomic location.
   - Finds corresponding SNP files, handling multiple files per SeqID.
4. **SNP File Analysis**: Processes each SNP file:
   - Filters SNPs using the p-value threshold.
   - Tags SNPs as 'cis' or 'trans' relative to the target protein.
   - Calculates distances for proximity flagging.
   - Adds the nearest gene information from the mapping file (in this case, the protein coding gene).
5. **Annotating Trans SNPs**: Determines the nearest gene for trans SNPs using genomic interval analysis with `pybedtools`, referencing the gene list txt file.
6. **Combining Results**: Aggregates processed data into a DataFrame.
7. **Output Saving**: Stores the final DataFrame as a pickle file.

**Important Note**: GWAS p-values should be in scientific notation (e.g., '1.2345e-678').

## Running the Script
Execute the script with the following command:
```
python cis_trans_tagger.py <path_to_config_file>
```

