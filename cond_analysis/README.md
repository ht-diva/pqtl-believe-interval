# Conditional Analysis

This scripts are meant to run a conditional analysis based on the `pygwas` packages.
Customization for slurm cluster and within process parallelization exploiting file separation in 
chromosomes. 

Script `run_cond_analysis.sh` will be the submitter, while the core code is `cond_analysis.py`.


## Arguments

### Process Submitter

- `RES_PATH`: location of the summary statistic normally where the `<pheno>.regenie.gz` files are located.
- `OUT_PATH`: where to store the conditional analysis results. If the directory does not exist it will be created.
- `LOG_PATH`: slurm log files for the script.
- `GENDATA` : location of a `json` file containing the genetic individual data to use for conditional analysis. 
- `IMPPANEL`: imputation panel key defined in the `${GENDATA}` json file
