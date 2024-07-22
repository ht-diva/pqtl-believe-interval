#!/bin/bash
#SBATCH --job-name=grm_part
### #SBATCH --partition=cpuq
### #SBATCH --array=0-100
#SBATCH --cpus-per-task=1
#SBATCH --time=300:00:00
#SBATCH --mem=5GB
#SBATCH --output=grm_part.out
#SBATCH --error=grm_part.err

##part_num=${SLURM_ARRAY_TASK_ID}
##echo $part_num

source /center/healthds/singularity_functions

# To use GCTA to estimate the heritability accounted for by all autosomal genome-wide SNPs, you need to first estimate the GRM, and then use the GRM to estimate the (SNP) heritability.
# Partition the GRM into 100 parts and allocate 8GB memory to each job
# The total memory required is approximately [n * (n + 1) / 2 * 12] / 10243 GB + 0.5GB, where n is the sample size, so it should be 14646.36 
# https://humanprod.service-now.com/sp?id=kb_article&sys_id=4c6557ca871bf010cb700d86cebb351c
# OUT_DIR=/home/c.giambartolomei/trials/believe

OUT_DIR=/group/diangelantonio/users/claudia/Claudia_Interval_test/

by_chrom=false
merge_all=false
make_sparse=true
make_single_pheno_files=false
run_analyses=false
clump=false

if [ "$by_chrom" = true ] ; then
# mem=3GB
cd ${OUT_DIR}/grm/
test=/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/cleaned_genotype_INTERVAL

# We first calculate the GRMs for all additional chromosomes:
gcta --bfile $test --chr 1 --make-grm-bin --out GCTAgrmchr1
gcta --bfile $test --chr 2 --make-grm-bin --out GCTAgrmchr2
gcta --bfile $test --chr 3 --make-grm-bin --out GCTAgrmchr3
gcta --bfile $test --chr 4 --make-grm-bin --out GCTAgrmchr4
gcta --bfile $test --chr 5 --make-grm-bin --out GCTAgrmchr5
gcta --bfile $test --chr 6 --make-grm-bin --out GCTAgrmchr6
gcta --bfile $test --chr 7 --make-grm-bin --out GCTAgrmchr7
gcta --bfile $test --chr 8 --make-grm-bin --out GCTAgrmchr8
gcta --bfile $test --chr 9 --make-grm-bin --out GCTAgrmchr9
gcta --bfile $test --chr 10 --make-grm-bin --out GCTAgrmchr10
gcta --bfile $test --chr 11 --make-grm-bin --out GCTAgrmchr11
gcta --bfile $test --chr 12 --make-grm-bin --out GCTAgrmchr12
gcta --bfile $test --chr 13 --make-grm-bin --out GCTAgrmchr13
gcta --bfile $test --chr 14 --make-grm-bin --out GCTAgrmchr14
gcta --bfile $test --chr 15 --make-grm-bin --out GCTAgrmchr15
gcta --bfile $test --chr 16 --make-grm-bin --out GCTAgrmchr16
gcta --bfile $test --chr 17 --make-grm-bin --out GCTAgrmchr17
gcta --bfile $test --chr 18 --make-grm-bin --out GCTAgrmchr18
gcta --bfile $test --chr 19 --make-grm-bin --out GCTAgrmchr19
gcta --bfile $test --chr 20 --make-grm-bin --out GCTAgrmchr20
gcta --bfile $test --chr 21 --make-grm-bin --out GCTAgrmchr21
gcta --bfile $test  --chr 22 --make-grm-bin --out GCTAgrmchr22

#for chr in {1..22}; do
#gcta --bfile $test --make-grm-part 100 $part_num --thread-num 5 --out GCTAgrm_v2_; done
fi

# GCTA LOCO 
# Mixed linear models were computed using GCTA (26) (v1.92.3 with the—mlma option using GRM as computed above)
if [ "$merge_all" = true ] ; then
# mem=3GB
cd ${OUT_DIR}/grm/
# cd /group/diangelantonio/users/claudia/Claudia_Believe_test/grm/
# Merge all the parts together:
# cat GCTAgrm_v2_*.grm.id > GCTAgrm_v2.grm.id
# cat GCTAgrm_v2_*.grm.bin > GCTAgrm_v2.grm.bin
# cat GCTAgrm_v2_*.grm.N.bin > GCTAgrm_v2.grm.N.bin
#cat GCTAgrmchr*.grm.id > GCTAgrmALL.grm.id
#cat GCTAgrmchr*.grm.bin > GCTAgrmALL.grm.bin
#cat GCTAgrmchr*.grm.N.bin > GCTAgrmALL.grm.N.bin
# To estimate the GRMs from the SNPs on each chromosome, then merge them by the command
gcta  --mgrm multipleGRMs.txt  --make-grm  --out GCTAgrmALL
fi

if [ "$make_sparse" = true ] ; then
# cd /group/diangelantonio/users/claudia/Claudia_Believe_test/grm/
cd ${OUT_DIR}/grm/
#Now we can create a sparse GRM from this re-calculated version of the full- dense GRM:
gcta --grm GCTAgrmALL --make-bK-sparse 0.025 --out newsparsegrm_0.025

#Check that the top few lines of newsparsegrm.grm.sp seem to match GCTAsparsegrm.grm.sp:
#head *.sp
#--grm-cutoff 0.025
#Remove one of a pair of individuals with estimated relatedness larger than the specified cut-off value (e.g. 0.025). GCTA selectively removes individuals to maximize the remaining sample size rather than doing it at random. NOTE: When merging multiple GRMs, this option does not apply to each single GRM but to the final merged GRM.

fi

if [ "$make_single_pheno_files" = true ] ; then
cd /group/diangelantonio/users/claudia/Claudia_Believe_test/phenos
# /Users/c.giambartolomei/Library/CloudStorage/OneDrive-Htechnopole/Cambridge_dropbox/Test_targets/targets_to_test_imputed_IVN_basic_covar.Rds
#*# pheno_file="/group/diangelantonio/users/claudia/Claudia_Believe_test/targets_to_test_imputed_IVN_basic_covar.Rds"
#*# pheno = readRDS(pheno_file)
# So, you just need the FID, IID, and Phenotype.
#*# p1 = pheno[,c(1,2,3)]
#*# write.table(p1,file="phenos.txt",col.names=F,row.names=F,quote=F) 

# cut -f1,2,6 test.fam > test.phen
fi

if [ "$run_analyses" = true ] ; then
#base_dir=/group/diangelantonio/users/claudia/Claudia_Believe_test
out_dir=${OUT_DIR}/out_test
pheno_dir=${base_dir}/phenos
grm_dir=${base_dir}/grm
test=/center/healthds/pQTL/BELIEVE/Genetic_QC_files/BELIEVE_genotype_final_bed_123456_withoutPCoutliers_mac
cd $grm_dir

# gcta64 --reml --mgrm-bin multipleGRMs.txt --pheno phenos.txt --out GCTAherit22GRMs
# Note this command makes use of a file multipleGRMs.txt which we created for you in advance, listing the stem names of the individual GRM files.
# We then run the analysis:
# simple LMM: gcta --mlma --bfile $test --pheno phenos/phenos.txt --thread-num 10 --out out_test/GCTAresults_p1
# gcta --reml --mgrm-bin multipleGRMs.txt --pheno ${pheno_dir}/phenos.txt --out ${out_dir}/GCTAherit22GRMs
#*# gcta --mlma --bfile $test --mgrm multipleGRMs.txt --pheno ${pheno_dir}/phenos.txt --out ${out_dir}/GCTAmlmaGRMs --thread-num 10 
# MLMA-LOCO analysis
for chr in {1..22}; do
gcta --mlma --grm GCTAgrmALL --mlma-subtract-grm GCTAgrmchr${chr} --bfile $test --chr ${chr} --pheno ${pheno_dir}/phenos.txt --out ${out_dir}/GCTAmlmaGRMs_loco_chr${chr} --thread-num 10
done


fi



if [ "$clump" = true ] ; then
# https://github.com/JoniColeman/gwas_scripts/blob/master/README.md
#SNP Clumping to identify independent hits
#Limit associations to lowest p-value in each region of linkage disequilibrium

$plink \
--bfile $test \
--clump $base_dir.post_imputation_final_analysis.assoc.logistic \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.25 \
--clump-kb 250 \
--out $root.post_imputation_final_analysis_clumped
#--clump-p1 is the p-value threshold below which to consider SNPs for inclusion as the reported SNP from the clump
#--clump-p2 is the p-value threshold below which to consider SNPs for inclusion in the clump
#--clump-r2 is the LD R2 threshold above which SNPs must be to be included in the same clump
#--clump-kb is the maximum distance a clump SNP can be from the reported SNP
fi
