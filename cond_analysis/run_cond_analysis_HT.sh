#!/bin/bash

# Base path of the project
PROJPATH="/scratch/${USER}/pQTL_somalogic"

# Create the result path based on the imputation panel
RES_PATH=${PROJPATH}/results_somalogic_${IMPPANEL}_int/results
OUT_PATH=${RES_PATH}/conditional_analysis
LOG_PATH=${PROJPATH}/cond_log
SCRIPTS_PATH=${PROJPATH}/scripts

# Genetic data json configuration
GENDATA=${PROJPATH}/genetic_data.json

# A key present in the ${GENDATA} json file
IMPPANEL="HRC13K"

# Create directories
mkdir -p ${PROJPATH}
mkdir -p ${OUT_PATH}
mkdir -p ${LOG_PATH}
mkdir -p ${SCRIPTS_PATH}


# Start submitting the scripts, one for each phenotype
for cc in `ls ${RES_PATH} | grep .regenie.gz`
do
    # Get phenotype
    p=${cc%%.regenie.gz}

    # Outfile
    ofile=${OUT_PATH}/${p}_cond_analysis_merged.cma.cojo
    echo ${p}

    # NB if the Conditional Analysis has been already run
    # do not resubmit the analysis
    if [ ! -e ${ofile} ]; then
         cat << EOF > ./scripts/${p}.sbatch
#!/bin/bash
#SBATCH --job-name=cond_analysis
#SBATCH --mail-type=NONE
#SBATCH --mail-user=${USER}@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --output=${LOG_PATH}/logs/%j_tdb_store.log
#SBATCH --mem-per-cpu=10G
#SBATCH --time=08:00:00

source /exchange/healthds/singularity_functions
time cond_analysis ${p} ${RES_PATH} -g ${IMPPANEL} -c 22  --conf ${GENDATA} --outpath ${OUT_PATH}
EOF

      echo "Processing ${p}"
      sbatch ./scripts/${p}.sbatch;
    fi
    sleep 2
done
