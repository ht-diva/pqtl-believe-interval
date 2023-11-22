#!/bin/bash

PROJPATH="/scratch/mfilosi/pQTL_somalogic"

# Define imputation panel
IMPPANEL="HRC13K"

# Create the result path based on the imputation panel
RES_PATH=${PROJPATH}/results_somalogic_${IMPPANEL}_int/results
OUT_PATH=${RES_PATH}/conditional_analysis
OUT_PATH=${PROJPATH}/conditional_analysis

if [ ! -e ${OUT_PATH} ]
then
    mkdir -p ${OUT_PATH}
fi

for cc in `ls ${RES_PATH} | grep .regenie.gz.tbi | grep x0so5582`
do
    p=${cc%%.regenie.gz.tbi}
    ofile=${OUT_PATH}/${p}_cond_analysis_merged.cma.cojo
    echo ${p} # rdsfile=${proj_path}/prot_pept_analysis/RDS/${p}_peptides_protclump.rds
    if [ ! -e ${ofile} ]
    then
        echo "Running ${p}"
        sbatch --nice=100 \
            --mem-per-cpu=8G \
            --partition=fast \
            -c 23 \
            --job-name=conditional_analysis \
            -o ${PROJPATH}/cond_log/cond_analysis_%j_${p}.log \
            ${PROJPATH}/bin/cond_analysis.py ${p} ${RES_PATH} -g ${IMPPANEL} -c 22  --conf ${PROJPATH}/genetic_data.json --outpath ${OUT_PATH}
    fi
done

