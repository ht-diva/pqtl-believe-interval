# GWAS analysis by CHRIS

This analysis has been done at Institute of Biomedicine (Eurac Research). 
In the following document there will be the description and the configuration files
used to run the GWAS analysis.

## Software and packages

GWAS has run with the `regenie` pipeline (see more info [here](https://rgcgithub.github.io/regenie/overview/)) 
with the `nextflow` implementation (see [here](https://genepi.github.io/nf-gwas))
which allows automatic and standardize quality control based on configuration 
parameters.

Installation of `nextflow` on the slurm cluster has been done by out sysAdmin.
Once the framework is installed and available on the frontend and on all nodes
it's possible to launch the pipeline.

### Configuration 

This configuration parameters are CHRIS specific and should be adjusted for your
particular locations and set of data. A complete list of parameters to set will be available [here](https://genepi.github.io/nf-gwas/)


#### Parameters regarding models, phenotype and covariates

File: [`chris_somalogic_pGWAS.conf`](chris_somalogic_pGWAS.conf)

```{java}
params.project = 'somalogic'
params.phenotypes_binary_trait = false

// Phenotype file
// It should start with two columns
// FID\tIID such as the fam file in plink format.
params.phenotypes_filename = "<phenotypes_file>"

// Phenotype to run in the same batch
params.phenotypes_columns = 'x0so0034,x0so0035'
params.phenotypes_delete_missings = true

// Transform using Rank Inverse Normal
params.phenotypes_apply_rint = false

// Define here where all the results will be output
params.outdir="<output_directory>"

// maf/mac parameters
params.qc_maf=0.01
params.qc_mac=100
params.regenie_min_mac=10

// Include additional parameters
// See below
includeConfig "./regenie_base_conf_HRC13K.config"
```


#### Parameters for genotype locations

File: `./regenie_base_conf_HRC13K.config`

```{java}
// Define additional parameters for the processes
params {
    genotypes_build               = 'hg19'
    genotypes_prediction          = '<genotyped_data_plink>{bim,bed,fam}'

    // Files in our case are splitted by chromosome
    // that's why we have the * before the extension.
    genotypes_association         = '<imputed_vcf_files>*.dose.vcf.gz'
    genotypes_association_format  = 'vcf'
    regenie_test                  = 'additive'
    annotation_min_log10p         = 5
    prune_enabled                 = false 
    prune_maf                     = 0.01
    prune_window_kbsize           = 1000
    prune_step_size               = 100
    prune_r2_threshold            = 0.9
}
```

#### Parameters for cluster specific job execution

The following parameters are recommended for imputed dataset. Some of them may
change based on the specific cluster configuration.

```
process{

    executor='slurm'
    // queue='slow'
    // clusterOptions = '--time=48:00:00'

    // Set up CPUs and Memory for plink process
    withLabel: 'process_plink2' {
        cpus   = 8
        memory = 64.GB
    }
    // Set up memory for creating report
    withLabel: 'required_memory_report' {
        memory = 32.GB
    }

    // Set up memory for merging results
    withName: 'MERGE_RESULTS' {
        memory = 16.GB
    }

    // Set up cpus and memory for the whole process
    withName: 'REGENIE_STEP1|REGENIE_STEP2' {
        cpus   = 8
        memory = 16.GB
    }
}
```

### Run the pipeline

From the front-end
```{bash}
nextflow run genepi/gwas-regenie -r v1.0.0 -profile slurm -c chris_somalogic_pGWAS.conf
```

Submitting to a node

```{bash}
sbatch --job-name=gwas --wrap `nextflow run genepi/gwas-regenie -r v1.0.0 -profile slurm -c chris_somalogic_pGWAS.conf`
```

# References

1. [https://rgcgithub.github.io/regenie/overview/](https://rgcgithub.github.io/regenie/overview/) 
2. [https://genepi.github.io/nf-gwas/](https://genepi.github.io/nf-gwas/) 
