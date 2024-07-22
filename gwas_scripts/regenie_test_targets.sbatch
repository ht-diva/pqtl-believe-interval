#!/bin/bash
#SBATCH --job-name test_%A
#SBATCH --output test_%A.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --time 5:00:00

module load nextflow/22.10.1 singularity/3.8.5
###  For a single project ###
export NXF_OPTS="-Xms1G -Xmx8G"
nextflow run /exchange/healthds/pQTL/regenie/nf-pipeline-regenie \
  -profile singularity,ht_cluster -c single_project_targets_1_2.conf

