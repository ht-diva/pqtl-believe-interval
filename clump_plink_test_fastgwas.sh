#!/bin/bash
#SBATCH --job-name=clump_fastgwas
#SBATCH --partition=cpuq
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=100GB
#SBATCH --output=clump_fastgwas-%A.out
#SBATCH --error=clump_fastgwas-%A.err

source /center/healthds/singularity_functions
# define output dirs
OUT_DIR=/exchange/healthds/pQTL/results/BELIEVE/fastGWA/post_imputation/
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

fastgwas_paths=/exchange/healthds/pQTL/pQTL_workplace/fastgwas_paths.txt

write_pathsFile=false
if [ "$write_pathsFile" = true ] ; then
	ASSOC_DIR=/group/diangelantonio/users/claudia/Claudia_Believe_test/fastGWAS_believe_selectedTargets_models2/
	find $ASSOC_DIR/ -maxdepth 1 -type f -not -path '*/\.*' -name "*fastGWA" | sort > $fastgwas_paths
fi

cat $fastgwas_paths | while read assocFile; do
echo $assocFile
seqid=$(basename $assocFile | cut -f2,3 -d'.' | cut -f1 -d'_')
model=$(basename $assocFile | cut -f2,3 -d'.' | cut -f2- -d'_')
echo $seqid
echo $model


convert_PGEN_BED=false
if [ "$convert_PGEN_BED" = true ] ; then
	GENO_DIR=/center/healthds/pQTL/BELIEVE/Genetic_QC_files/
	sample=/center/healthds/pQTL/BELIEVE/Genetic_QC_files/bgen_from_HDS_BELIEVE_final_HDS_erased_phase.sample
        plink2 --pgen $GENO_DIR/HDS_BELIEVE_final_HDS.pgen --pvar $GENO_DIR/HDS_BELIEVE_final_HDS.pvar --psam $GENO_DIR/HDS_BELIEVE_final_HDS.psam --keep $sample --extract /center/healthds/pQTL/BELIEVE/Working_shared/HDS_BELIEVE_final_HDS.snplist --id-delim ':' --threads 8  --memory 32768 --make-bed --out /group/diangelantonio/users/claudia/Claudia_Believe_test/fastGWAS_believe_selectedTargets_models2/post_imputation/HDS_BELIEVE_final_HDS_no_rel_snps
fi

GENO_DIR=/group/diangelantonio/users/claudia/Claudia_Believe_test/fastGWAS_believe_selectedTargets_models2/post_imputation/
# identified LD-independent loci using PLINK clumping function (parameters: -clump-p1 = 5*10^-8, -clump-p2 = 0.05, -clump-r2 = 0.4, -clump-kb = 500). 
snp_field=SNP
field=P

# these should be set and decided
p1=5e-08
p2=1 #5e-08
r2=0.01
kb=1000 #500

clumpedFile=${OUT_DIR}/clumped_seq.${seqid}_${model}_p1${p1}_p2${p2}_r2${r2}_kb${kb}.tab

/exchange/healthds/.bin/plink2 \
--bfile $GENO_DIR/HDS_BELIEVE_final_HDS_no_rel_snps \
--clump  $assocFile \
--clump-id-field ${snp_field} \
--clump-p-field ${field} \
--clump-p1 ${p1} \
--clump-p2 ${p2} \
--clump-r2 ${r2} \
--clump-kb ${kb} \
--out $clumpedFile
#--clump-p1 is the p-value threshold below which to consider SNPs for inclusion as the reported SNP from the clump
#--clump-p2 is the p-value threshold below which to consider SNPs for inclusion in the clump
#--clump-r2 is the LD R2 threshold above which SNPs must be to be included in the same clump
#--clump-kb is the maximum distance a clump SNP can be from the reported SNP
#fi

done

