#!/bin/bash
#SBATCH --job-name=fastgwas_believe
#SBATCH --partition=cpuq
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --output=fastgwas_believe-%A.out
#SBATCH --error=fastgwas_believe-%A.err

source /center/healthds/singularity_functions

OUT_DIR=/exchange/healthds/pQTL/results/BELIEVE/fastGWA
GENO_DIR=/center/healthds/pQTL/BELIEVE/Genetic_QC_files/
GRM_DIR=/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/GRM 
PHENO=/exchange/healthds/pQTL/BELIEVE/Proteomics_QC_files/BELIEVEL_NonImp_residuals_new.txt ## here the path for the new residuals file 

OUT_PHENO_TEMP=${OUT_DIR}/pheno
mkdir -p $OUT_PHENO_TEMP

colNum=3
# you just need the FID, IID, and Phenotype.

nm=$(head -n 1 $PHENO | cut -f $colNum)
echo $colNum
echo $nm
echo $SLURM_JOBID

cut -f1,2,$colNum $PHENO > ${OUT_PHENO_TEMP}/phen_${nm}

gcta --pfile $GENO_DIR/HDS_BELIEVE_final_HDS --grm-sparse $GRM_DIR/newsparsegrm_0.025 --fastGWA-mlm --pheno ${OUT_PHENO_TEMP}/phen_${nm} --thread-num 10 --out ${OUT_DIR}/geno_assoc_phen_${nm}

# print memory usage at the end of the job to stdout:
#seff $SLURM_JOBID

clean=false
if [ "$clean" = true ] ; then
	rm -r ${OUT_PHENO_TEMP}
fi

# check if all finished ok?
# cd /home/c.giambartolomei/trials/scripts/mixed_models/temp_scripts_selectedTargts_models
# if everything went well, these should be the same:
# grep "Analysis finished" *.out | wc -l
# ls *.out | wc -l
