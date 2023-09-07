#source /center/healthds/singularity_functions

SRC_DIR=/center/healthds/pQTL/BELIEVE/Genetic_QC_files/
seqid=11178.21
model=basic_covar

# define input dirs
GENO_DIR=/center/healthds/pQTL/BELIEVE/Genetic_QC_files/
ASSOC_DIR=/group/diangelantonio/users/claudia/Claudia_Believe_test/fastGWAS_believe_selectedTargets_models2/
assocFile=${ASSOC_DIR}/geno_assoc_phen_seq.${seqid}_${model}.fastGWA

# define output dirs
OUT_DIR=/group/diangelantonio/users/claudia/Claudia_Believe_test/fastGWAS_believe_selectedTargets_models2/post_imputation/
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

################################################ DO THIS ONCE OR USE MY FILES
convert_PGEN_BED=false
if [ "$convert_PGEN_BED" = true ] ; then

sample=/center/healthds/pQTL/BELIEVE/Genetic_QC_files/bgen_from_HDS_BELIEVE_final_HDS_erased_phase.sample

plink2  \
--pgen $SRC_DIR/HDS_BELIEVE_final_HDS.pgen  \
--pvar $SRC_DIR/HDS_BELIEVE_final_HDS.pvar  \
--psam $SRC_DIR/HDS_BELIEVE_final_HDS.psam  \
--keep $sample  \
--extract /center/healthds/pQTL/BELIEVE/Working_shared/HDS_BELIEVE_final_HDS.snplist  \
--id-delim ':'  \
--threads 8  --memory 32768  \
--make-bed  \
--out $OUT_DIR/HDS_BELIEVE_final_HDS_no_rel_snps

fi

################################################
# define fields (different for each software)
snp_field=SNP # this could be different in Regenie results
field=P # this could be different in Regenie results

# define fields (these should be set and decided after testing)
p1=5e-08
p2=1
r2=0.1
kb=1000 #500

clumpedFile=${OUT_DIR}/clumped_seq.${seqid}_${model}_p1${p1}_p2${p2}_r2${r2}_kb${kb}.tab

plink \
--bfile $OUT_DIR/HDS_BELIEVE_final_HDS_no_rel_snps \
--clump  $assocFile \
--clump-snp-field ${snp_field} \
--clump-field ${field} \
--clump-p1 ${p1} \
--clump-p2 ${p2} \
--clump-r2 ${r2} \
--clump-kb ${kb} \
--out $clumpedFile

#--clump-p1 is the p-value threshold below which to consider SNPs for inclusion as the reported SNP from the clump
#--clump-p2 is the p-value threshold below which to consider SNPs for inclusion in the clump
#--clump-r2 is the LD R2 threshold above which SNPs must be to be included in the same clump
#--clump-kb is the maximum distance a clump SNP can be from the reported SNP
