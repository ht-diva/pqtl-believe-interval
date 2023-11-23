#!/bin/bash

# Base path of the project
PROJPATH="/scratch/mfilosi/pQTL_somalogic"

# Create the result path based on the imputation panel
RES_PATH=${PROJPATH}/results_somalogic_${IMPPANEL}_int/results
OUT_PATH=${RES_PATH}/conditional_analysis
LOG_PATH=${PROJPATH}/cond_log

# Genetic data json configuration
GENDATA=${PROJPATH}/genetic_data.json

# A key present in the ${GENDATA} json file
IMPPANEL="HRC13K"

# Create directories
if [ ! -e ${OUT_PATH} ]
then
    mkdir -p ${OUT_PATH}
fi

if [ ! -e ${LOG_PATH} ]
then
    mkdir -p ${LOG_PATH}
fi

# Start submitting the scripts, one for each phenotype
for cc in `ls ${RES_PATH} | grep .regenie.gz.tbi | grep x0so5582`
do
    # Get phenotype
    p=${cc%%.regenie.gz.tbi}

    # Outfile
    ofile=${OUT_PATH}/${p}_cond_analysis_merged.cma.cojo
    echo ${p}

    # NB if the Conditional Analysis has been already run
    # do not resubmit the analysis
    if [ ! -e ${ofile} ]
    then
        echo "Running ${p}"
        sbatch --nice=100 \
            --mem-per-cpu=8G \
            --partition=fast \
            -c 23 \
            --job-name=conditional_analysis \
            -o ${LOG_PATH}/cond_analysis_%j_${p}.log \
            ${PROJPATH}/bin/cond_analysis.py ${p} ${RES_PATH} -g ${IMPPANEL} -c 22  --conf ${GENDATA} --outpath ${OUT_PATH}
    fi
done

